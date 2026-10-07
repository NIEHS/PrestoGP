# Extract the Matern (theta) parameters from a PrestoGP model.

This method extracts a list containing the Matern variance/covariance
parameters from a PrestoGP model.

## Usage

``` r
# S4 method for class 'PrestoGPModel'
get_theta(model)
```

## Arguments

- model:

  The PrestoGP model object

## Value

A list containing the following elements:

- sigma::

  A vector containing the estimates of the variance parameter sigma for
  each outcome.

- scale::

  A vector containing the estimates of the the scale parameter(s) for
  each outcome.

- smoothness::

  A vector containing the estimates of the smoothness parameter for each
  outcome.

- nuggets::

  A vector containing the estimates of the nugget variance for each
  outcome.

- correlation::

  A matrix containing the estimates of the correlations between
  outcomes. Omitted for univariate models.

## References

- Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
  land-use regression for ground-level nitrogen dioxide", The Annals of
  Applied Statistics (2021) 15(2):688-710.

- Porcu, E., Bevilacqua, M., Schaback, R., and Oates RJ. "The Matérn
  Model: A Journey Through Statistics, Numerical Analysis and Machine
  Learning", Statistical Science (2024) 39(3):469-492.

## See also

[`PrestoGPModel-class`](https://niehs.github.io/PrestoGP/reference/PrestoGPModel-class.md),
[`prestogp_fit`](https://niehs.github.io/PrestoGP/reference/prestogp_fit-PrestoGPModel-method.md)

## Examples

``` r
data(soil)
soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
y <- soil[,4]                   # predict moisture content
X <- as.matrix(soil[,5:9])
locs <- as.matrix(soil[,1:2])

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
#> Current penalized negative log likelihood: 487.5317 
#> Current MSE: 9.104869 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 482.4013 
#> Current MSE: 9.05704 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 482.3368 
#> Current MSE: 9.057425 
#> Beginning iteration 4 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 4 complete 
#> Current penalized negative log likelihood: 482.3368 
#> Current MSE: 9.055174 
get_theta(soil.vm)
#> $sigma
#> [1] 10.26467
#> 
#> $scale
#> [1] 17.71528
#> 
#> $smoothness
#> [1] 0.7058993
#> 
#> $nuggets
#> [1] 0.7606272
#> 
```
