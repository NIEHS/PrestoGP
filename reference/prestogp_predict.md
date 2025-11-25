# Prediction for PrestoGP models

After fitting a PrestoGP model, this method can be used to make
predictions on an independent test set (consisting of both new locations
and new values of the predictor variables).

## Usage

``` r
prestogp_predict(
  model,
  X = "matrix",
  locs = "matrix",
  m = "numeric",
  ordering.pred = c("obspred", "general"),
  pred.cond = c("independent", "general"),
  return.values = c("mean", "meanvar")
)

# S4 method for class 'VecchiaModel'
prestogp_predict(
  model,
  X = "matrix",
  locs = "matrix",
  m = NULL,
  ordering.pred = c("obspred", "general"),
  pred.cond = c("independent", "general"),
  return.values = c("mean", "meanvar")
)

# S4 method for class 'FullModel'
prestogp_predict(
  model,
  X = "matrix",
  locs = "matrix",
  m = NULL,
  ordering.pred = c("obspred", "general"),
  pred.cond = c("independent", "general"),
  return.values = c("mean", "meanvar")
)

# S4 method for class 'MultivariateVecchiaModel'
prestogp_predict(
  model,
  X = "matrix",
  locs = "matrix",
  m = NULL,
  ordering.pred = c("obspred", "general"),
  pred.cond = c("independent", "general"),
  return.values = c("mean", "meanvar")
)
```

## Arguments

- model:

  A PrestoGP model object obtained after running `link{prestogp_fit}`.

- X:

  The values of the predictor variable(s) for which prediction will be
  performed. Should be a matrix for univariate models or a list for
  multivariate models.

- locs:

  The locations where prediction will be performed. Should be a matrix
  for univariate models or a list for multivariate models.

- m:

  The number of neighbors to condition on. If not specified, it will
  default to the value of m used to fit the model.

- ordering.pred:

  Should "obspred" or "general" ordering be used for prediction? See
  [`vecchia_specify`](https://rdrr.io/pkg/GPvecchia/man/vecchia_specify.html)
  or
  [`vecchia_Mspecify`](https://niehs.github.io/PrestoGP/reference/vecchia_Mspecify.md).
  Defaults to "obspred".

- pred.cond:

  Should prediction conditioning be "general" or "indepedent"? See
  [`vecchia_specify`](https://rdrr.io/pkg/GPvecchia/man/vecchia_specify.html)
  or
  [`vecchia_Mspecify`](https://niehs.github.io/PrestoGP/reference/vecchia_Mspecify.md).
  Defaults to "independent".

- return.values:

  Values that should be returned. Possible values include "mean" and
  "meanvar". See
  [`vecchia_prediction`](https://rdrr.io/pkg/GPvecchia/man/vecchia_prediction.html)
  or
  [`vecchia_Mprediction`](https://niehs.github.io/PrestoGP/reference/vecchia_Mprediction.md).
  Defaults to "mean".

## Value

A list containing the estimated mean values and (if requested) standard
deviations for the predictions.

## Details

It is important to note that the variance estimates produced by this
function assume that the Matern covariance parameters are fixed and
known. It does not consider any variance resulting from estimating the
covariance parameters or the regression coefficients in PrestoGP models.
These variance estimates will be anticonservative in such cases (and a
warning will be returned when these estimates are calculated).

Prediction is currently not implemented for full models. This function
will return an error if it is applied to a full model.

## References

- Katzfuss, M., and Guinness, J. "A general framework for Vecchia
  approximations of Gaussian processes", Statistical Science (2021)
  36(1):124-141.

- Katzfuss, M., Guinness, J., Gong, W. and Zilber, D. "Vecchia
  approximations of Gaussian-process predictions", Journal of
  Agricultural, Biological and Environmental Statistics (2020)
  25:383-414.

- Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
  land-use regression for ground-level nitrogen dioxide", The Annals of
  Applied Statistics (2021) 15(2):688-710.

## See also

[`PrestoGPModel-class`](https://niehs.github.io/PrestoGP/reference/PrestoGPModel-class.md),[`prestogp_fit`](https://niehs.github.io/PrestoGP/reference/prestogp_fit-PrestoGPModel-method.md),
[`vecchia_prediction`](https://rdrr.io/pkg/GPvecchia/man/vecchia_prediction.html),
[`vecchia_Mprediction`](https://niehs.github.io/PrestoGP/reference/vecchia_Mprediction.md)

## Examples

``` r
data(soil)
soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
y <- soil[,4]                   # predict moisture content
X <- as.matrix(soil[,5:9])
locs <- as.matrix(soil[,1:2])

# Create training and test sets
n <- length(y)
otr <- rep(FALSE, n)
otr[sample(1:n, size=floor(n/2))] <- TRUE
otst <- !otr

# Fit the model on the training set
soil.vm <- new("VecchiaModel", n_neighbors = 10)
soil.vm <- prestogp_fit(soil.vm, y[otr], X[otr,], locs[otr,])
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
#> Current penalized negative log likelihood: 253.7902 
#> Current MSE: 9.779005 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 253.0128 
#> Current MSE: 9.689043 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 252.2248 
#> Current MSE: 9.683512 
#> Beginning iteration 4 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 4 complete 
#> Current penalized negative log likelihood: 252.2248 
#> Current MSE: 9.692544 

# Perform predictions on the test set
soil.yhat <- prestogp_predict(soil.vm, X[otst,], locs[otst,])
```
