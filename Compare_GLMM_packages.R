# Comparison of computational time for generalized linear mixed effects models
# Author: Fabio Sigrist, May 2021
# Comparison has been done using 'gpboost' version 0.6.3, 'lme4' version 1.1-27, 'statsmodels' version 0.12.2

library(gpboost)
library(lme4)
library(reticulate)
# py_install("statsmodels", pip=TRUE)
path <- "C:/GLMM_comparison/"

sim_GLMM_data <- function(n, m, randef_type, likelihood, 
                          num_covariates = 10, sigma2 = 1){
  ## Purpose: Simulate data from (generalized) linear mixed effects models
  ## -----------------------------------------------------------
  ## Arguments: n = number of data points
  ## m = number of random groups / categories for (first) random effect
  ## randef_type = type of random effects
  ## likelihood = likelihood / data distribution
  ## num_covariates = number of covariates in linear predictor
  ## sigma2 = variance of random effects
  ## --------------------------------------------------------------
  ## Value: a list with 3 entries: y, X, group_data
  ## Author: Fabio Sigrist, Date: May 28, 2021
  ## --------------------------------------------------------------
  
  # Simulate random effects
  group <- rep(1,n) 
  for(i in 1:m) group[((i-1)*n/m+1):(i*n/m)] <- i
  b1 <- sqrt(sigma2) * rnorm(m)
  if(randef_type == "One_random_effect"){
    eps <- b1[group]
    group_data <- group
  } else if(randef_type == "Two_completely_crossed_random_effects"){
    n_obs_gr <- n/m # number of samples per group
    group2 <- rep(1,n) # grouping variable for second crossed random effect
    for(i in 1:m) group2[(1:n_obs_gr)+n_obs_gr*(i-1)] <- 1:n_obs_gr
  } else if(randef_type == "Two_randomly_crossed_random_effects"){
    group2 <-  group[sample.int(n=n, size=n, replace=FALSE)]
  } else if(randef_type == "Two_nested_random_effects"){
    m_nested <- m*2 # number of categories / levels for the second nested grouping variable
    group2 <- rep(1,n)  # grouping variable for nested lower level random effects
    for(i in 1:m_nested) group2[((i-1)*n/m_nested+1):(i*n/m_nested)] <- i
  }
  if(randef_type != "One_random_effect"){
    b2 <- sqrt(sigma2) * rnorm(length(unique(group2)))
    if (length(unique(group2)) != max(group2)) stop("not all levels samples -> gives index problems")
    eps <- b1[group] + b2[group2]
    group_data <- cbind(group,group2)
  }
  # Simulate fixed effects
  beta <- c(0,rep(1,num_covariates))
  X <- cbind(rep(1,n),matrix(rnorm(num_covariates*n),nrow=n))
  lp <- as.vector(X%*%beta)
  delta_var_lp <- sqrt(sigma2/var(lp))
  X[,-1] <- X[,-1] * delta_var_lp
  lp <- as.vector(X%*%beta)
  colnames(X) <- c("1",paste0("Cov_",1:num_covariates))
  # Simulate response variable
  if (likelihood == "bernoulli_probit") {
    probs <- pnorm(lp + eps)
    y <- as.numeric(runif(n) < probs)
  } else if (likelihood == "bernoulli_logit") {
    probs <- 1/(1+exp(-(lp + eps)))
    y <- as.numeric(runif(n) < probs)
  } else if (likelihood == "poisson") {
    mu <- exp(lp + eps)
    y <- qpois(runif(n), lambda = mu)
  } else if (likelihood == "gamma") {
    mu <- exp(f + eps)
    y <- qgamma(runif(n), scale = mu, shape = 1)
  } else if (likelihood == "gaussian") {
    mu <- lp + eps
    y <- sqrt(sigma2/5) * rnorm(n) + mu
  }
  group_data <- as.data.frame(group_data)
  colnames(group_data) <- paste0("RE_",1:dim(group_data)[2])
  return(list(y=y, X=X, group_data=group_data))
}

run_simulation_experiment <- function(nsim, sigma2, likelihood, path, 
                                      ndata_try, n_obs_per_group_try, 
                                      randef_type_try, num_covariates_try) {
  ## Purpose: Perform simulation study / experiments for generalized linear mixed effects models
  ## -----------------------------------------------------------
  ## Arguments: nsim = number of simulation runs
  ## sigma2 = variance of random effects
  ## likelihood = likelihood / data distribution
  ## path = Path where to save data for runnig Python code via reticulate
  ## ndata_try = vector with number of samples
  ## n_obs_per_group_try = vector with number of samples per group
  ## randef_type_try = vector with types of random effects
  ## num_covariates_try = vector with number of covariates
  ## --------------------------------------------------------------
  ## Value: a data.frame with the results from the simulated experiments
  ## Author: Fabio Sigrist, Date: May 28, 2021
  ## --------------------------------------------------------------
  
  results <- data.frame()
  for (n_obs_per_group in n_obs_per_group_try) {
    for (randef_type in randef_type_try) {
      if (randef_type == "One_random_effect") {
        vcs_true <- sigma2
      } else {
        vcs_true <- rep(sigma2,2)
      }
      for (n in ndata_try) {
        for (num_covariates in num_covariates_try) {
          coef_true <- c(0,rep(1,num_covariates))
          for (iter in 1:nsim) {
            
            m <- n / n_obs_per_group
            if (as.integer(m) != m) stop("'n' is not divisible by 'n_obs_per_group'")
            print(paste0("********* randef_type = ",randef_type,", n = ", n,", m = ", m,", num_covariates = ",
                         num_covariates,", iter = ", iter, " *********"))
            
            # Simulate data
            sim_data <- sim_GLMM_data(n = n, m = m, randef_type = randef_type, likelihood = likelihood,
                                      num_covariates = num_covariates, sigma2 = sigma2)
            y <- sim_data$y
            X <- sim_data$X
            group_data <- sim_data$group_data
            # Export for python
            write.csv(y,file = paste0(path,"y.csv"), row.names = FALSE)
            write.csv(X,file = paste0(path,"X.csv"), row.names = FALSE)
            write.csv(group_data,file = paste0(path,"group_data.csv"), row.names = FALSE)
            write.csv(likelihood,file = paste0(path,"likelihood.csv"), row.names = FALSE)
            
            #-----------------gpboost----------------- 
            t1 <- Sys.time()
            gp_model <- fitGPModel(group_data = group_data, y = y, X = X, likelihood = likelihood)
            # Note: gpboost uses Nesterov accelerated gradient descent by default.
            #  This can be changed to "Nelder-Mead" using e.g. 'params=list(optimizer_cov="nelder_mead",maxit=10000)'
            t2 <- Sys.time()
            ctime <- signif(difftime(t2, t1, units = "secs"), 3)
            # summary(gp_model)
            mse_coef <- mean((coef_true - gp_model$get_coef())^2)
            mse_vcs <- mean((vcs_true - gp_model$get_cov_pars())^2)
            results <- rbind(results,c("gpboost",randef_type,n,m,num_covariates,iter,
                                       ctime,mse_coef,mse_vcs))
            print("gpboost finished")
            
            #-----------------lme4----------------- 
            formula = as.formula(paste0("y ~ ",paste(colnames(X), collapse = ' + ')," + ",paste("(1|",colnames(group_data),")", collapse = ' + ')))
            if (likelihood=="bernoulli_probit") {
              family <- binomial(link = "probit")
            } else if (likelihood=="poisson") {
              family <- poisson(link = "log")
            }
            t1 <- Sys.time()
            mod <- glmer(formula, data=data.frame(y=y,X,group_data),family=family)
            t2 <- Sys.time()
            ctime <- signif(difftime(t2, t1, units = "secs"), 3)
            # summary(mod)
            vcs_fit <- c()
            for (ii in 1:length(vcs_true)) vcs_fit[ii] <- summary(mod)$varcor[[colnames(group_data)[ii]]][1]
            mse_coef <- mean((coef_true - summary(mod)$coefficients[,1])^2)
            mse_vcs <- mean((vcs_true - vcs_fit)^2)
            results <- rbind(results,c("lme4",randef_type,n,m,num_covariates,iter,
                                       ctime,mse_coef,mse_vcs))
            print("lme4 finished")
            
            #-----------------statsmodels-----------------
            source_python(paste0("GLMM_statsmodels.py"))
            results <- rbind(results,c("statsmodels",randef_type,n,m,num_covariates,iter,
                                       time_statsmodels,mse_coefs_statsmodels,mse_vcs_statsmodels))
            print("statsmodels finished")
          } # end loop over nsim
        } # end loop over num_covariates_try
      } # end loop over ndata_try
    } # end loop over randef_type_try
  } # end loop over n_obs_per_group_try
  colnames(results) <- c("package", "randef_type", "n", "m" ,"num_covariates","iter", "time", "mse_coefs", "mse_vcs")
  for (ic in 3:dim(results)[2]) results[,ic] <- as.numeric(results[,ic])
  return(results)
}

nsim <- 100
sigma2 <- 1

###############################
## Varying number of samples
###############################
ndata_try <- c(100,200,500,1000,2000)
randef_type_try <- c("One_random_effect")
num_covariates_try <- c(10)
n_obs_per_group_try <- c(10)
likelihood <- "bernoulli_probit"

set.seed(1)
results_sample_size <- run_simulation_experiment(nsim = nsim, sigma2 = sigma2, likelihood = likelihood,
                                                 n_obs_per_group_try = n_obs_per_group_try,
                                                 randef_type_try = randef_type_try, ndata_try = ndata_try, 
                                                 num_covariates_try = num_covariates_try, path = path)
write.csv(results_sample_size,file = paste0(path,"results/results_sample_size.csv"), row.names = FALSE)

###############################
## Varying number of covariates
###############################
ndata_try <- c(1000)
randef_type_try <- c("One_random_effect")
num_covariates_try <- c(1,2,5,10,20)
n_obs_per_group_try <- c(10)
likelihood <- "bernoulli_probit"

set.seed(1)
results_num_covariates <- run_simulation_experiment(nsim = nsim, sigma2 = sigma2, likelihood = likelihood,
                                                    n_obs_per_group_try = n_obs_per_group_try,
                                                    randef_type_try = randef_type_try, ndata_try = ndata_try, 
                                                    num_covariates_try = num_covariates_try, path = path)
write.csv(results_num_covariates,file = paste0(path,"results/results_num_covariates.csv"), row.names = FALSE)

###############################
## Varying number of groups
###############################
ndata_try <- c(1000)
randef_type_try <- c("One_random_effect")
num_covariates_try <- c(10)
n_obs_per_group_try <- c(2,5,10,20,50,100,200,500)
likelihood <- "bernoulli_probit"

set.seed(1)
results_number_groups <- run_simulation_experiment(nsim = nsim, sigma2 = sigma2, likelihood = likelihood,
                                                   n_obs_per_group_try = n_obs_per_group_try,
                                                   randef_type_try = randef_type_try, ndata_try = ndata_try, 
                                                   num_covariates_try = num_covariates_try, path = path)
write.csv(results_number_groups,file = paste0(path,"results/results_number_groups.csv"), row.names = FALSE)

###############################
## Different types of random effects
###############################
ndata_try <- c(1000)
randef_type_try <- c("One_random_effect","Two_completely_crossed_random_effects",
                     "Two_randomly_crossed_random_effects","Two_nested_random_effects")
num_covariates_try <- c(10)
n_obs_per_group_try <- c(10)
likelihood <- "bernoulli_probit"

set.seed(1)
results_randef_type <- run_simulation_experiment(nsim = nsim, sigma2 = sigma2, likelihood = likelihood,
                                                 n_obs_per_group_try = n_obs_per_group_try,
                                                 randef_type_try = randef_type_try, ndata_try = ndata_try, 
                                                 num_covariates_try = num_covariates_try, path = path)
write.csv(results_randef_type,file = paste0(path,"results/results_randef_type.csv"), row.names = FALSE)

###############################
## Other likelihood: varying number of samples
###############################
ndata_try <- c(100,200,500,1000,2000)
randef_type_try <- c("One_random_effect")
num_covariates_try <- c(10)
n_obs_per_group_try <- c(10)
likelihood <- "poisson"

set.seed(1)
results_poisson_sample_size <- run_simulation_experiment(nsim = nsim, sigma2 = sigma2, likelihood = likelihood,
                                                         n_obs_per_group_try = n_obs_per_group_try,
                                                         randef_type_try = randef_type_try, ndata_try = ndata_try, 
                                                         num_covariates_try = num_covariates_try, path = path)
write.csv(results_poisson_sample_size,file = paste0(path,"results/results_poisson_sample_size.csv"), row.names = FALSE)



###############################
## Plot results
###############################
results_sample_size <- read.csv(file = paste0(path,"results/results_sample_size.csv"))
results_num_covariates <- read.csv(file = paste0(path,"results/results_num_covariates.csv"))
results_number_groups <- read.csv(file = paste0(path,"results/results_number_groups.csv"))
results_randef_type <- read.csv(file = paste0(path,"results/results_randef_type.csv"))
results_poisson_sample_size <- read.csv(file = paste0(path,"results/results_poisson_sample_size.csv"))

library(ggplot2)
library(gridExtra)
height <- 8
height_all <- 4

## Varying number of samples
p1 <- ggplot(data=results_sample_size, aes(x=n,y=time,color=package,shape=package)) + geom_jitter(alpha=0.1,width=0.01,height=0) +
  stat_summary(fun=mean, geom="point",size=3,stroke=2.5,show.legend=FALSE) +
  stat_summary(fun=mean, geom="line",size=1,show.legend=FALSE) + 
  ylab("Time (sec)") + xlab("") +
  scale_y_log10() + scale_x_log10() + theme(plot.title=element_text(hjust=0.5),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=16)) +
  theme(legend.position = "none")

p2 <- ggplot(data=results_sample_size, aes(x=n,y=mse_coefs,color=package,shape=package)) + geom_jitter(alpha=0.1,width=0.01,height=0) +
  stat_summary(fun=mean, geom="point",size=3,stroke=2.5,show.legend=FALSE) +
  stat_summary(fun=mean, geom="line",size=1,show.legend=FALSE) +
  ylab("MSE coefficients") + xlab("Number samples") +
  scale_y_log10() + scale_x_log10() + theme(plot.title=element_text(hjust=0.5),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=16)) +
  theme(legend.position = "none")

p3 <- ggplot(data=results_sample_size, aes(x=n,y=mse_coefs,color=package,shape=package)) + geom_jitter(alpha=0.1,width=0.01,height=0) +
  stat_summary(fun=mean, geom="point",size=3,stroke=2.5,show.legend=FALSE) +
  stat_summary(fun=mean, geom="line",size=1,show.legend=FALSE) +
  ylab("MSE variance components") + xlab("") +
  scale_y_log10() + scale_x_log10() + theme(plot.title=element_text(hjust=0.5),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=16),
                                            legend.text=element_text(size=14),
                                            legend.title=element_text(size=16)) +
  guides(colour=guide_legend(override.aes=list(alpha=1,size=3,stroke=2.5),title="Package"),
         shape=guide_legend(title="Package"))

pall <- grid.arrange(p1, p2, p3, ncol=3, widths=c(1,1,1.45))
ggsave(pall,file=paste0(path,"plots/","Results_sample_size.jpeg"),height=height_all,width=4*height_all)


## Varying number of covariates
p1 <- ggplot(data=results_num_covariates, aes(x=num_covariates,y=time,color=package,shape=package)) + geom_jitter(alpha=0.1,width=0.01,height=0) +
  stat_summary(fun=mean, geom="point",size=3,stroke=2.5,show.legend=FALSE) +
  stat_summary(fun=mean, geom="line",size=1,show.legend=FALSE) +
  ylab("Time (sec)") + xlab("") +
  scale_y_log10() + scale_x_log10() + theme(plot.title=element_text(hjust=0.5),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=16)) +
  theme(legend.position = "none")

p2 <- ggplot(data=results_num_covariates, aes(x=num_covariates,y=mse_coefs,color=package,shape=package)) + geom_jitter(alpha=0.1,width=0.01,height=0) +
  stat_summary(fun=mean, geom="point",size=3,stroke=2.5,show.legend=FALSE) +
  stat_summary(fun=mean, geom="line",size=1,show.legend=FALSE) +
  ylab("MSE coefficients") + xlab("Number of covariates") +
  scale_y_log10() + scale_x_log10() + theme(plot.title=element_text(hjust=0.5),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=16)) +
  theme(legend.position = "none")

p3 <- ggplot(data=results_num_covariates, aes(x=num_covariates,y=mse_vcs,color=package,shape=package)) + geom_jitter(alpha=0.1,width=0.01,height=0) +
  stat_summary(fun=mean, geom="point",size=3,stroke=2.5,show.legend=FALSE) +
  stat_summary(fun=mean, geom="line",size=1,show.legend=FALSE) +
  ylab("MSE variance components") + xlab("") +
  scale_y_log10() + scale_x_log10() + theme(plot.title=element_text(hjust=0.5),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=16),
                                            legend.text=element_text(size=14),
                                            legend.title=element_text(size=16)) +
  guides(colour=guide_legend(override.aes=list(alpha=1,size=3,stroke=2.5),title="Package"),
         shape=guide_legend(title="Package"))

pall <- grid.arrange(p1, p2, p3, ncol=3, widths=c(1,1,1.45))
ggsave(pall,file=paste0(path,"plots/","Results_num_covariates.jpeg"),height=height_all,width=4*height_all)


## Varying number of groups
p1 <- ggplot(data=results_number_groups, aes(x=m,y=time,color=package,shape=package)) + geom_jitter(alpha=0.1,width=0.01,height=0) +
  stat_summary(fun=mean, geom="point",size=3,stroke=2.5,show.legend=FALSE) +
  stat_summary(fun=mean, geom="line",size=1,show.legend=FALSE) +
  ylab("Time (sec)") + xlab("") +
  scale_y_log10() + scale_x_log10() + theme(plot.title=element_text(hjust=0.5),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=16)) +
  theme(legend.position = "none")

p2 <- ggplot(data=results_number_groups, aes(x=m,y=mse_coefs,color=package,shape=package)) + geom_jitter(alpha=0.1,width=0.01,height=0) +
  stat_summary(fun=mean, geom="point",size=3,stroke=2.5,show.legend=FALSE) +
  stat_summary(fun=mean, geom="line",size=1,show.legend=FALSE) +
  ylab("MSE coefficients") + xlab("Number of groups") +
  scale_y_log10() + scale_x_log10() + theme(plot.title=element_text(hjust=0.5),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=16)) +
  theme(legend.position = "none")

p3 <- ggplot(data=results_number_groups, aes(x=m,y=mse_vcs,color=package,shape=package)) + geom_jitter(alpha=0.1,width=0.01,height=0) +
  stat_summary(fun=mean, geom="point",size=3,stroke=2.5,show.legend=FALSE) +
  stat_summary(fun=mean, geom="line",size=1,show.legend=FALSE) +
  ylab("MSE variance components") + xlab("") + 
  scale_y_log10() + scale_x_log10() + theme(plot.title=element_text(hjust=0.5),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=16),
                                            legend.text=element_text(size=14),
                                            legend.title=element_text(size=16)) +
  guides(colour=guide_legend(override.aes=list(alpha=1,size=3,stroke=2.5),title="Package"),
         shape=guide_legend(title="Package"))

pall <- grid.arrange(p1, p2, p3, ncol=3, widths=c(1,1,1.45))
ggsave(pall,file=paste0(path,"plots/","Results_num_groups.jpeg"),height=height_all,width=4*height_all)


## Varying types of random effects
randef_type_try <- c("One_random_effect","Two_completely_crossed_random_effects",
                     "Two_randomly_crossed_random_effects","Two_nested_random_effects")
new_names <- c("One RE","Two crossed REs","Two crossed REs v.2","Two nested REs")
for (i in 1:length(randef_type_try)){
  results_randef_type$randef_type[results_randef_type$randef_type==randef_type_try[i]] = new_names[i]
}

p1 <- ggplot(data=results_randef_type, aes(x=randef_type,y=time,color=package,shape=package)) + geom_jitter(alpha=0.1,width=0.01,height=0) +
  stat_summary(fun=mean, geom="point",size=3,stroke=2.5,show.legend=FALSE) +
  ylab("Time (sec)") + xlab("") +
  scale_y_log10() + theme(plot.title=element_text(hjust=0.5),
                          axis.text=element_text(size=12),
                          axis.title=element_text(size=16),
                          axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  theme(legend.position = "none")

p2 <- ggplot(data=results_randef_type, aes(x=randef_type,y=mse_coefs,color=package,shape=package)) + geom_jitter(alpha=0.1,width=0.01,height=0) +
  stat_summary(fun=mean, geom="point",size=3,stroke=2.5,show.legend=FALSE) +
  ylab("MSE coefficients") + xlab("Random effect type") + 
  scale_y_log10() + theme(plot.title=element_text(hjust=0.5),
                          axis.text=element_text(size=12),
                          axis.title=element_text(size=16),
                          axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  theme(legend.position = "none")

p3 <- ggplot(data=results_randef_type, aes(x=randef_type,y=mse_vcs,color=package,shape=package)) + geom_jitter(alpha=0.1,width=0.01,height=0) +
  stat_summary(fun=mean, geom="point",size=3,stroke=2.5,show.legend=FALSE) +
  ylab("MSE variance components") + xlab("") +
  scale_y_log10() + theme(plot.title=element_text(hjust=0.5),
                          legend.text=element_text(size=14),
                          legend.title=element_text(size=16),
                          axis.text=element_text(size=12),
                          axis.title=element_text(size=16),
                          axis.text.x=element_text(angle=45, hjust=1,vjust=1)) +
  guides(colour=guide_legend(override.aes=list(alpha=1,size=3,stroke=2.5),title="Package"),
         shape=guide_legend(title="Package"))

pall <- grid.arrange(p1, p2, p3, ncol=3, widths=c(1,1,1.45))
ggsave(pall,file=paste0(path,"plots/","Results_randef_type.jpeg"),height=height_all,width=3.5*height_all)


## Poisson likelihood
p1 <- ggplot(data=results_poisson_sample_size, aes(x=n,y=time,color=package,shape=package)) + geom_jitter(alpha=0.1,width=0.01,height=0) +
  stat_summary(fun=mean, geom="point",size=3,stroke=2.5,show.legend=FALSE) +
  stat_summary(fun=mean, geom="line",size=1,show.legend=FALSE) +
  ylab("Time (sec)") + xlab("") + 
  scale_y_log10() + scale_x_log10() + theme(plot.title=element_text(hjust=0.5),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=16)) +
  theme(legend.position = "none")

p2 <- ggplot(data=results_poisson_sample_size, aes(x=n,y=mse_coefs,color=package,shape=package)) + geom_jitter(alpha=0.1,width=0.01,height=0) +
  stat_summary(fun=mean, geom="point",size=3,stroke=2.5,show.legend=FALSE) +
  stat_summary(fun=mean, geom="line",size=1,show.legend=FALSE) +
  ylab("MSE coefficients") + xlab("Number samples") +
  scale_y_log10() + scale_x_log10() + theme(plot.title=element_text(hjust=0.5),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=16)) +
  theme(legend.position = "none")

p3 <- ggplot(data=results_poisson_sample_size, aes(x=n,y=mse_coefs,color=package,shape=package)) + geom_jitter(alpha=0.1,width=0.01,height=0) +
  stat_summary(fun=mean, geom="point",size=3,stroke=2.5,show.legend=FALSE) +
  stat_summary(fun=mean, geom="line",size=1,show.legend=FALSE) +
  ylab("MSE variance components") + xlab("") +
  scale_y_log10() + scale_x_log10() + theme(plot.title=element_text(hjust=0.5),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=16),
                                            legend.text=element_text(size=14),
                                            legend.title=element_text(size=16)) +
  guides(colour=guide_legend(override.aes=list(alpha=1,size=3,stroke=2.5),title="Package"),
         shape=guide_legend(title="Package"))

pall <- grid.arrange(p1, p2, p3, ncol=3, widths=c(1,1,1.45))
ggsave(pall,file=paste0(path,"plots/","Results_poisson_sample_size.jpeg"),height=height_all,width=4*height_all)


