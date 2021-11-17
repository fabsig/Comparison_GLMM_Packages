# Comparing Software Packages for Generalized Linear Mixed Effects Models (GLMMs)

This repository contains the code to reproduce the full simulation study described in [this blog post](https://towardsdatascience.com/generalized-linear-mixed-effects-models-in-r-and-python-with-gpboost-89297622820c). 

The simulation is done by the file [Compare_GLMM_packages.R](https://github.com/fabsig/Comparison_GLMM_Packages/blob/master/Compare_GLMM_packages.R). The `reticulate` package is used to call the Python code in [GLMM_statsmodels.py](https://github.com/fabsig/Comparison_GLMM_Packages/blob/master/GLMM_statsmodels.py), which runs the `statsmodels` GLMMs, from R.
