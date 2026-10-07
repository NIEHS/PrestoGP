# Extract the regression coefficients (beta) for a PrestoGP model.

This method extracts a list containing the regression coefficients for a
PrestoGP model.

## Usage

``` r
# S4 method for class 'PrestoGPModel'
get_beta(model)
```

## Arguments

- model:

  The PrestoGP model object

## Value

A list containing the regression coefficients. Each element of the list
corresponds to an outcome variable. The final element of the list
contains the intercept(s).

## Details

It is important to note that the intercepts are estimated on the
transformed data. They are not useful for making predictions on the
original (untransformed) data. Use
[`prestogp_predict`](https://niehs.github.io/PrestoGP/reference/prestogp_predict.md)
to make predictions based on a PrestoGP model.

## References

- Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
  land-use regression for ground-level nitrogen dioxide", The Annals of
  Applied Statistics (2021) 15(2):688-710.

## See also

[`PrestoGPModel-class`](https://niehs.github.io/PrestoGP/reference/PrestoGPModel-class.md),
[`prestogp_fit`](https://niehs.github.io/PrestoGP/reference/prestogp_fit-PrestoGPModel-method.md),
[`prestogp_predict`](https://niehs.github.io/PrestoGP/reference/prestogp_predict.md)

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
#> Current penalized negative log likelihood: 487.387 
#> Current MSE: 9.104869 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 481.3886 
#> Current MSE: 9.048082 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 481.3435 
#> Current MSE: 9.038025 
#> Beginning iteration 4 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 4 complete 
#> Current penalized negative log likelihood: 481.2817 
#> Current MSE: 9.048029 
#> Beginning iteration 5 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 5 complete 
#> Current penalized negative log likelihood: 481.2817 
#> Current MSE: 9.043287 
get_beta(soil.vm)
#> $Y
#>       NO3.N     Total.N       NH4.N         DOC        N20N 
#> -0.03702712  0.00000000  0.02874770  0.00240729 34.35233338 
#> 
#> $`(Intercept)`
#> (Intercept) 
#>    11.39463 
#> 
```
