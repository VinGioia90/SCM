## SCM: Smoothing for Covariance matrix Modelling ##

The `SCM` R package offers methods for fitting multivariate Gaussian additive models, where the unconstrained entries of a covariance matrix parametrisation are allowed to vary with smooth functions of the covariates. 
The models are fitted using the `gam_scm()` wrapper,  which uses the model fitting routines of the `mgcv` package. The diagnostic and model visualising tools of the `mgcViz` package can be used. See the vignettes for an introduction to the main functionality of the  `SCM` package.
