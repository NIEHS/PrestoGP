# show method for PrestoGP models

This method prints a summary of a PrestoGP model and its parameters.

## Usage

``` r
# S4 method for class 'PrestoGPModel'
show(object)
```

## Arguments

- object:

  The PrestoGP model object

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
#> Current penalized negative log likelihood: 487.4952 
#> Current MSE: 9.104869 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 482.95 
#> Current MSE: 9.019932 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 482.3916 
#> Current MSE: 9.052201 
#> Beginning iteration 4 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 4 complete 
#> Current penalized negative log likelihood: 482.1548 
#> Current MSE: 9.042599 
#> Beginning iteration 5 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 5 complete 
#> Current penalized negative log likelihood: 482.0701 
#> Current MSE: 9.048612 
#> Beginning iteration 6 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 6 complete 
#> Current penalized negative log likelihood: 482.0552 
#> Current MSE: 9.0508 
#> Beginning iteration 7 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 7 complete 
#> Current penalized negative log likelihood: 482.0552 
#> Current MSE: 9.0508 
show(soil.vm)
#> Matern covariance parameters (theta): 
#> $sigma
#> [1] 10.51871
#> 
#> $scale
#> [1] 15.58912
#> 
#> $smoothness
#> [1] 0.7962241
#> 
#> $nuggets
#> [1] 0.8344075
#> 
#> Regression coefficients (beta): 
#> $Y
#>        NO3.N        NH4.N          DOC         N20N 
#> -0.039103489  0.029995180  0.002587154 28.703457422 
#> 
#> $`(Intercept)`
#> (Intercept) 
#>    11.39708 
#> 
#> Model type: VecchiaModel 
#> Nearest neighbors: 10 
#> Scaling: 1 1 
#> Penalized likelihood: 482.0552 
#> MSE: 9.0508 
```
