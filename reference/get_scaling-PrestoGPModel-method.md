# Extract the scaling parameter for a PrestoGP model.

This method extracts the scaling parameter for a PrestoGP model.

## Usage

``` r
# S4 method for class 'PrestoGPModel'
get_scaling(model)
```

## Arguments

- model:

  The PrestoGP model object

## Value

A vector containing the scaling parameter. See
[`prestogp_fit`](https://niehs.github.io/PrestoGP/reference/prestogp_fit-PrestoGPModel-method.md)
for an explanation of this parameter.

## References

- Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
  land-use regression for ground-level nitrogen dioxide", The Annals of
  Applied Statistics (2021) 15(2):688-710.

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
#> Current penalized negative log likelihood: 487.4296 
#> Current MSE: 9.104869 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 483.5328 
#> Current MSE: 9.033187 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 482.8815 
#> Current MSE: 9.056389 
#> Beginning iteration 4 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 4 complete 
#> Current penalized negative log likelihood: 482.7846 
#> Current MSE: 9.044995 
#> Beginning iteration 5 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 5 complete 
#> Current penalized negative log likelihood: 482.3403 
#> Current MSE: 9.058457 
#> Beginning iteration 6 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 6 complete 
#> Current penalized negative log likelihood: 482.3403 
#> Current MSE: 9.048252 
get_scaling(soil.vm)
#> [1] 1 1
```
