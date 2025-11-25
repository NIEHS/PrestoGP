# Extract the Y values for a PrestoGP model.

This method extracts the values of Y for a PrestoGP model. If imputation
was used when fitting the model, the missing Y's will be replaced with
their imputed values.

## Usage

``` r
get_Y(model)

# S4 method for class 'VecchiaModel'
get_Y(model)

# S4 method for class 'MultivariateVecchiaModel'
get_Y(model)
```

## Arguments

- model:

  The PrestoGP model object

## Value

A vector or list containing the values of Y. Any missing Y's will be
replaced with their imputed values.

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
#> Current penalized negative log likelihood: 487.4966 
#> Current MSE: 9.104869 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 482.5612 
#> Current MSE: 9.028533 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 482.2526 
#> Current MSE: 9.046949 
#> Beginning iteration 4 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 4 complete 
#> Current penalized negative log likelihood: 482.2526 
#> Current MSE: 9.037678 
get_Y(soil.vm)
#>   [1]  7.233333 14.633333 12.666667 13.733333  7.400000 15.966667 10.266667
#>   [8]  8.033333  9.633333  7.433333  8.966667  7.933333  8.200000  7.700000
#>  [15]  8.933333  7.066667 15.000000  9.766667 10.500000 13.366667 13.700000
#>  [22] 12.866667 22.300000 14.600000 15.566667 17.066667 10.700000  7.233333
#>  [29] 10.700000 12.133333 14.066667  9.700000  8.866667  7.833333  8.633333
#>  [36]  7.566667 13.233333 13.633333 11.233333 12.700000  9.400000  6.966667
#>  [43]  7.533333 12.533333  9.200000  8.833333  9.466667 10.133333  6.933333
#>  [50]  8.533333  8.300000  4.566667  8.450000 11.100000  6.050000 10.950000
#>  [57]  8.700000 14.050000 14.500000 10.250000 11.100000  7.300000 16.200000
#>  [64] 14.950000 13.400000 11.200000  8.750000 13.500000 11.350000 11.050000
#>  [71]  9.500000 10.300000 10.100000  9.200000  8.550000  9.600000 10.600000
#>  [78]  7.233333  9.100000  8.633333  9.366667  8.333333  9.300000 10.766667
#>  [85] 13.266667 15.433333 10.733333 11.400000  9.733333  9.866667 13.166667
#>  [92] 13.633333 11.100000 10.633333 12.400000 11.200000  9.433333 14.000000
#>  [99] 12.300000 15.300000 14.333333 10.400000  9.666667 14.450000 14.833333
#> [106] 13.666667  9.233333 11.850000 12.900000 12.466667 15.566667 13.266667
#> [113] 10.700000 10.666667 12.533333 12.566667 14.366667 12.333333 12.566667
#> [120] 12.600000 11.800000 12.866667  8.225000 12.600000  8.200000  8.850000
#> [127] 12.250000  8.500000  9.450000  8.300000 13.550000  9.650000 11.150000
#> [134] 13.650000 12.800000  8.500000 11.700000 11.150000  5.550000 10.000000
#> [141] 10.100000 10.100000  9.400000  7.550000  8.700000  8.900000  9.700000
#> [148]  6.950000  8.650000 11.400000 12.350000 12.050000 11.950000 11.850000
#> [155] 10.600000 12.450000 13.500000 13.100000  9.300000 11.600000 10.650000
#> [162]  9.300000  9.150000  9.350000 10.733333 16.950000 10.900000 12.766667
#> [169] 15.433333  9.833333 12.466667 13.450000 14.266667 11.566667  6.266667
#> [176] 15.466667 12.500000 13.233333 10.566667 12.166667  8.866667  8.933333
#> [183]  9.066667  8.466667 17.966667 14.800000 19.366667 19.300000 13.733333
#> [190] 13.600000 20.800000  7.966667  8.833333 20.466667 14.700000 15.000000
#> [197] 14.533333 12.066667 18.633333 19.100000
```
