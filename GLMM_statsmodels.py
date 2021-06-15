# -*- coding: utf-8 -*-
"""
Comparison of computational time for generalized linear mixed effects models
Author: Fabio Sigrist, May 2021
"""

import pandas as pd
import numpy as np
import os
import time
import statsmodels.genmod.bayes_mixed_glm as glm

path_data = "C:\\GLMM_comparison\\"

# load data
group_data_P = pd.read_csv(os.path.join(path_data, 'group_data.csv'),squeeze=True)
X_P = pd.read_csv(os.path.join(path_data, 'X.csv'),squeeze=True)
y_P = pd.read_csv(os.path.join(path_data, 'y.csv'),squeeze=True)
likelihood_P = pd.read_csv(os.path.join(path_data, 'likelihood.csv'),squeeze=True)[0]

if(len(group_data_P.shape)==1):
    num_randeff = 1
else:
    num_randeff = group_data_P.shape[1]
num_coef = X_P.shape[1]
# Fit model
if likelihood_P == 'bernoulli_probit':
    model = glm.BinomialBayesMixedGLM(endog=y_P, exog=X_P, exog_vc=group_data_P,
                                      ident=np.arange(0,num_randeff))
elif likelihood_P == 'poisson':
    model = glm.PoissonBayesMixedGLM(endog=y_P, exog=X_P, exog_vc=group_data_P,
                      ident=np.arange(0,num_randeff))
start = time.time()
model = model.fit_map()
end = time.time()
time_statsmodels = (end - start)
# Extract fitted values
summary = model.summary()
coefs = summary.tables[0]['Post. Mean'][0:num_coef]
vcs = np.exp(summary.tables[0]['Post. Mean'][num_coef:]) ** 2
coefs = coefs.astype(float)
vcs = vcs.astype(float)
# Calculate MSEs
coefs_true = pd.Series([0]).append(pd.Series(np.ones(num_coef-1)))
vcs_true = pd.Series(np.ones(num_randeff))
mse_coefs_statsmodels = np.mean((coefs.values - coefs_true.values)**2)
mse_vcs_statsmodels = np.mean((vcs.values - vcs_true.values)**2)
