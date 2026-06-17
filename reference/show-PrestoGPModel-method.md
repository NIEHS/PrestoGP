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
#> Current penalized negative log likelihood: 487.6106 
#> Current MSE: 9.104869 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 482.1937 
#> Current MSE: 9.057586 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 482.1937 
#> Current MSE: 9.046636 
show(soil.vm)
#> Matern covariance parameters (theta): 
#> $sigma
#> [1] 10.63738
#> 
#> $scale
#> [1] 14.44042
#> 
#> $smoothness
#> [1] 0.796426
#> 
#> $nuggets
#> [1] 0.7330837
#> 
#> Regression coefficients (beta): 
#> $Y
#>        NO3.N        NH4.N          DOC         N20N 
#> -0.037611577  0.028946406  0.002497739 34.120803768 
#> 
#> $`(Intercept)`
#> (Intercept) 
#>    11.39514 
#> 
#> Model type: VecchiaModel 
#> Nearest neighbors: 10 
#> Scaling: 1 1 
#> Penalized likelihood: 482.1937 
#> MSE: 9.046636 
```
