# Train a PrestoGP model.

This method fits a PrestoGP model given a set of locations and predictor
and outcome variables.

## Usage

``` r
# S4 method for class 'PrestoGPModel'
prestogp_fit(
  model,
  Y,
  X,
  locs,
  Y.names = NULL,
  X.names = NULL,
  scaling = NULL,
  common_scale = NULL,
  covparams = NULL,
  beta.hat = NULL,
  tol = 0.999999,
  max.iters = 100,
  center.y = NULL,
  impute.y = FALSE,
  lod.upper = NULL,
  lod.lower = NULL,
  n.impute = 10,
  eps.impute = 0.01,
  maxit.impute = 0,
  quiet = FALSE,
  verbose = FALSE,
  optim.method = "Nelder-Mead",
  optim.control = list(trace = 0, reltol = 0.001, maxit = 5000),
  penalty = c("lasso", "relaxed", "MCP", "SCAD"),
  alpha = 1,
  family = c("gaussian", "binomial"),
  nfolds = 10,
  foldid = NULL,
  parallel = FALSE,
  cluster = NULL,
  adaptive = FALSE
)
```

## Arguments

- model:

  The PrestoGP model object being fit.

- Y:

  The values of the response variable(s). Should be a matrix or vector
  for univariate models or a list for multivariate models.

- X:

  The values of the predictor variable(s). Should be a matrix for
  univariate models or a list for multivariate models.

- locs:

  The values of the locations. Should be a matrix for univariate models
  or a list for multivariate models.

- Y.names:

  The name(s) for the response variable(s). Should be a vector with
  length equal to the number of response variables. Defaults to NULL.

- X.names:

  The names for the predictor variables. Should be a vector for
  univariate models or a list for multivariate models.

- scaling:

  A vector of consecutive positive integers that is used to specify
  which columns of locs should have the same scaling parameter. For
  example, in a spatiotemporal model with two spatial measures and a
  time measure, the value of scaling would be c(1, 1, 2). The length of
  scaling must match the number of columns of locs. If it is not
  specified, all columns of locs will have a common scale parameter.

- common_scale:

  Do all columsn of locs have a common scales parameter? See Details for
  the effects of this parameter. Defaults to TRUE if there is only one
  scale parameter for each outcome and FALSE otherwise.

- covparams:

  The initial covariance parameters estimate (optional).

- beta.hat:

  The initial beta parameters estimates (optional).

- tol:

  The model is considered converged when error is not less than
  tol\*previous_error (optional). Defaults to 0.999999.

- max.iters:

  Maximum number of iterations for the model fitting procedure. Defaults
  to 100.

- center.y:

  Should the Y's be mean centered before fitting the model? Defaults to
  TRUE for gaussian models and FALSE for binomial models.

- impute.y:

  Should missing Y's be imputed? Defaults to FALSE.

- lod.upper:

  Upper limit of detection value(s). Any missing Y is assumed to be less
  than lod.upper when performing missing data imputation. Should be
  numeric for univariate models and a list for multivariate models,
  where each element of the list corresponds to an outcome. The *i*th
  element of the lod is the limit of detection for observation i.
  Alternatively, one can specify a single limit of detection that is
  assumed to be the same for all observations. If not specified, it is
  assumed that no upper limit of detection exists. Ignored if impute.y
  is FALSE.

- lod.lower:

  Same as lod.upper, except it specifies the lower limit of detection:
  missing Y values are assumed to be greater than lod.lower. This is
  often equal to 0.

- n.impute:

  How many multiply imputed values of missing y's should be generated?
  Defaults to 10.

- eps.impute:

  The multiple imputation procedure will stop if successive coefficient
  estimates are within eps.impute of one another. Defaults to 0.01.

- maxit.impute:

  Maximum number of iterations for the multiple imputation procedure .
  Defaults to 0 (meaning that the multiple imputation procedure is not
  executed). See Details.

- quiet:

  If FALSE, the penalized log likelihood and the model MSE will be
  printed for each iteration of the model fitting procedure. No
  intermediate output will be printed if TRUE. Defaults to FALSE.

- verbose:

  If TRUE, the estimated theta/beta parameters will be printed for each
  iteration of the model fitting procedure. Defaults to FALSE. Ignored
  if quiet is TRUE.

- optim.method:

  Optimization method to be used for the maximum likelihood estimation
  that is passed to optim. Defaults to "Nelder-Mead". See
  [`optim`](https://rdrr.io/r/stats/optim.html).

- optim.control:

  Control parameter that is passed to optim. See
  [`optim`](https://rdrr.io/r/stats/optim.html).

- penalty:

  The type of penalized regression to be used. Should be one of "lasso",
  "relaxed", "MCP", or "SCAD". Note that "lasso" and "relaxed" will fit
  the model using glmnet and "MCP" and "SCAD" will fit the model using
  ncvreg. See
  [`glmnet`](https://glmnet.stanford.edu/reference/glmnet.html) or
  [`ncvreg`](https://pbreheny.github.io/ncvreg/reference/ncvreg.html).
  Defaults to "lasso".

- alpha:

  The elastic net mixing parameter. 'alpha=1' corresponds to lasso (or
  SCAD/MCP) penalty; 'alpha=0' corresponds to ridge regression. See
  [`glmnet`](https://glmnet.stanford.edu/reference/glmnet.html) or
  [`ncvreg`](https://pbreheny.github.io/ncvreg/reference/ncvreg.html).
  Defaults to 1.

- family:

  Family parameter for the glmnet or ncvreg model. Currently only
  "gaussian" and "binomial" are supported. Defaults to "gaussian". See
  [`glmnet`](https://glmnet.stanford.edu/reference/glmnet.html) or
  [`ncvreg`](https://pbreheny.github.io/ncvreg/reference/ncvreg.html).

- nfolds:

  Number of cross-validation folds for cv.glmnet. Defaults to 10. See
  [`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html).

- foldid:

  Optional vector of values between 1 and "nfolds" specifying what fold
  each observation should be assigned to in the cv.glmnet
  cross-validation. See
  [`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html)
  and
  [`cv.ncvreg`](https://pbreheny.github.io/ncvreg/reference/cv.ncvreg.html).

- parallel:

  Should parallel "foreach" be used to speed up the model fitting
  procedure where possible? Defaults to FALSE. Specifically,
  parallelization will be used for imputation and fitting the cv.glmnet
  object. See
  [`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html).
  Note that this only applies to glmnet models (where penalty="lasso" or
  penalty="relaxed"). Models using ncvreg (where penalty="MCP" or
  penalty="SCAD") require a cluster argument for parallelization (see
  below).

- cluster:

  A cluster for running cv.ncvreg in parallel. See
  [`cv.ncvreg`](https://pbreheny.github.io/ncvreg/reference/cv.ncvreg.html)
  and [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html). This
  must be specified to run cv.ncvreg in parallel. It is ignored for
  glmnet models.

- adaptive:

  Should adaptive lasso be used? Defaults to FALSE. It is ignored for
  SCAD and MCP models.

## Value

A PrestoGPModel object with slots updated based on the results of the
model fitting procedure. See
[`PrestoGPModel-class`](https://niehs.github.io/PrestoGP/reference/PrestoGPModel-class.md)
for details.

## Details

If common_scale is TRUE, multivariate models will use the Matern
cross-covariance function described in Apanasovich et al. (2012). This
model can only be used if each outcome has only a single scale
parameter. If common_scale is FALSE, each column of locs will be divded
by the corresponding element of scaling, and the cross-covariance will
be computed under the assumption that all scale parameters are equal
to 1. For univariate models, the Vecchia approximation will not be
recomputed in each iteration of the model fitting procedure if
common_scale is TRUE, but the parameter is otherwise ignored.

If impute.y is TRUE, missing values of y will be imputed. If there is no
limit of detection specified, then missing y's will be imputed using an
EM algorithm. If there are limits of detection, then missing y's will be
imputed using a multiple imputation algorithm based on the current
estimates of the Matern (theta) parameters. Unfortunately, this
procedure is not robust to inaccurate estimates of theta and produce
inaccurate results if theta is misspecified. Thus, the maximum number of
iterations is set to 0 by default, meaning that this step is skipped.
Results can be improved by increasing the number of iterations if one
has confidence in the theta estimates. If impute.y is FALSE, then any
missing y's will produce an error.

## References

- Apanasovich, T.V., Genton, M.G. and Sun, Y. "A valid Matérn class of
  cross-covariance functions for multivariate random fields with any
  number of components", Journal of the American Statistical
  Association (2012) 107(497):180-193.

- Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
  land-use regression for ground-level nitrogen dioxide", The Annals of
  Applied Statistics (2021) 15(2):688-710.

- Zou, H. "The adaptive lasso and its oracle properties", Journal of the
  American Statistical Association (2006) 101(476):1418-1429.

## See also

[`PrestoGPModel-class`](https://niehs.github.io/PrestoGP/reference/PrestoGPModel-class.md),
[`glmnet`](https://glmnet.stanford.edu/reference/glmnet.html)

## Examples

``` r
data(soil)
soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
y <- soil[,4]                   # predict moisture content
X <- as.matrix(soil[,5:9])
locs <- as.matrix(soil[,1:2])

# Vecchia model
soil.vm <- new("VecchiaModel", n_neighbors = 10)
soil.vm <- prestogp_fit(soil.vm, y, X, locs)
#> 
#> Estimating initial beta... 
#> Estimation of initial beta complete 
#> 
#> Beginning iteration 1 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 1 complete 
#> Current penalized negative log likelihood: 487.6744 
#> Current MSE: 9.104869 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 482.5802 
#> Current MSE: 9.038949 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 482.3463 
#> Current MSE: 9.040137 
#> Beginning iteration 4 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 4 complete 
#> Current penalized negative log likelihood: 482.3463 
#> Current MSE: 9.023641 

# Impute missing y's
miss <- sample(1:nrow(soil), 20)
y.miss <- y
y.miss[miss] <- NA
soil.vm2 <- new("VecchiaModel", n_neighbors = 10)
soil.vm2 <- prestogp_fit(soil.vm, y, X, locs, impute.y = TRUE)
#> 
#> Estimating initial beta... 
#> Estimation of initial beta complete 
#> 
#> Beginning iteration 1 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 1 complete 
#> Current penalized negative log likelihood: 487.3579 
#> Current MSE: 9.104869 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 482.0535 
#> Current MSE: 9.054379 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 481.7182 
#> Current MSE: 9.049372 
#> Beginning iteration 4 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 4 complete 
#> Current penalized negative log likelihood: 481.6796 
#> Current MSE: 9.051876 
#> Beginning iteration 5 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 5 complete 
#> Current penalized negative log likelihood: 481.6796 
#> Current MSE: 9.054146 

# Impute y's missing due to limit of detection
soil.lod <- quantile(y, 0.1)
y.lod <- y
y.lod[y.lod <= soil.lod] <- NA
soil.vm3 <- new("VecchiaModel", n_neighbors = 10)
soil.vm3 <- prestogp_fit(soil.vm, y, X, locs, impute.y = TRUE,
lod.upper = soil.lod)
#> 
#> Estimating initial beta... 
#> Estimation of initial beta complete 
#> 
#> Beginning iteration 1 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 1 complete 
#> Current penalized negative log likelihood: 487.4757 
#> Current MSE: 9.104869 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 483.0882 
#> Current MSE: 9.020892 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 482.2932 
#> Current MSE: 9.029623 
#> Beginning iteration 4 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 4 complete 
#> Current penalized negative log likelihood: 482.2406 
#> Current MSE: 9.033159 
#> Beginning iteration 5 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 5 complete 
#> Current penalized negative log likelihood: 482.0968 
#> Current MSE: 9.046215 
#> Beginning iteration 6 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 6 complete 
#> Current penalized negative log likelihood: 482.0968 
#> Current MSE: 9.043857 

# Full model
soil.fm <- new("FullModel")
soil.fm <- prestogp_fit(soil.fm, y, X, locs)
#> 
#> Estimating initial beta... 
#> Estimation of initial beta complete 
#> 
#> Beginning iteration 1 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 1 complete 
#> Current penalized negative log likelihood: 487.4422 
#> Current MSE: 9.104869 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 486.7817 
#> Current MSE: 9.105443 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 482.0513 
#> Current MSE: 9.099646 
#> Beginning iteration 4 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 4 complete 
#> Current penalized negative log likelihood: 482.0384 
#> Current MSE: 9.094437 
#> Beginning iteration 5 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 5 complete 
#> Current penalized negative log likelihood: 481.9948 
#> Current MSE: 9.094437 
#> Beginning iteration 6 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 6 complete 
#> Current penalized negative log likelihood: 481.9948 
#> Current MSE: 9.061991 

# Multivariate model
ym <- list()
ym[[1]] <- soil[,5]             # predict two nitrogen concentration levels
ym[[2]] <- soil[,7]
Xm <- list()
Xm[[1]] <- Xm[[2]] <- as.matrix(soil[,c(4,6,8,9)])
locsm <- list()
locsm[[1]] <- locsm[[2]] <- locs

soil.mvm <-  new("MultivariateVecchiaModel", n_neighbors = 10)
soil.mvm <- prestogp_fit(soil.mvm, ym, Xm, locsm)
#> 
#> Estimating initial beta... 
#> Estimation of initial beta complete 
#> 
#> Beginning iteration 1 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 1 complete 
#> Current penalized negative log likelihood: 1306.597 
#> Current MSE: 46.92982 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 1292.651 
#> Current MSE: 161.7518 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 1292.651 
#> Current MSE: 170.9654 

# Space/elevation model
data(soil250, package="geoR")
y2 <- soil250[,7]               # predict pH level
X2 <- as.matrix(soil250[,c(4:6,8:22)])
# columns 1+2 are location coordinates; column 3 is elevation
locs2 <- as.matrix(soil250[,1:3])

soil.vm2 <- new("VecchiaModel", n_neighbors = 10)
# fit separate scale parameters for location and elevation
soil.vm2 <- prestogp_fit(soil.vm2, y2, X2, locs2, scaling = c(1, 1, 2))
#> 
#> Estimating initial beta... 
#> Estimation of initial beta complete 
#> 
#> Beginning iteration 1 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 1 complete 
#> Current penalized negative log likelihood: -258.4624 
#> Current MSE: 0.008200979 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: -258.4624 
#> Current MSE: 0.009064528 
```
