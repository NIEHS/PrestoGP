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
#> Current penalized negative log likelihood: 487.6794 
#> Current MSE: 9.104869 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 482.4042 
#> Current MSE: 9.041971 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 481.9999 
#> Current MSE: 9.034115 
#> Beginning iteration 4 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 4 complete 
#> Current penalized negative log likelihood: 481.8401 
#> Current MSE: 9.047439 
#> Beginning iteration 5 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 5 complete 
#> Current penalized negative log likelihood: 481.8401 
#> Current MSE: 9.049394 
show(soil.vm)
#> Matern covariance parameters (theta): 
#> $sigma
#> [1] 10.28376
#> 
#> $scale
#> [1] 13.52572
#> 
#> $smoothness
#> [1] 0.9088149
#> 
#> $nuggets
#> [1] 0.7629514
#> 
#> Regression coefficients (beta): 
#> $Y
#>        NO3.N        NH4.N          DOC         N20N 
#> -0.039404924  0.030027624  0.002530095 35.497666810 
#> 
#> $`(Intercept)`
#> (Intercept) 
#>    11.39525 
#> 
#> Model type: VecchiaModel 
#> Nearest neighbors: 10 
#> Scaling: 1 1 
#> Penalized likelihood: 481.8401 
#> MSE: 9.049394 
```
