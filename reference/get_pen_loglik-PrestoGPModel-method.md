# Extract the penalized log likelihood for a PrestoGP model.

This method extracts the (negative) penalized log likelihood for a
PrestoGP model.

## Usage

``` r
# S4 method for class 'PrestoGPModel'
get_pen_loglik(model)
```

## Arguments

- model:

  The PrestoGP model object

## Value

The value of the (negative) penalized log likelihood. See References for
details.

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
#> Current penalized negative log likelihood: 487.6736 
#> Current MSE: 9.104869 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 482.6845 
#> Current MSE: 9.039046 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 482.5234 
#> Current MSE: 9.05557 
#> Beginning iteration 4 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 4 complete 
#> Current penalized negative log likelihood: 482.5095 
#> Current MSE: 9.059903 
#> Beginning iteration 5 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 5 complete 
#> Current penalized negative log likelihood: 482.5095 
#> Current MSE: 9.053168 
get_pen_loglik(soil.vm)
#> [1] 482.5095
```
