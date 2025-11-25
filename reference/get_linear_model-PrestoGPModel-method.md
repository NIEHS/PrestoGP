# Extract the fitted linear model for a PrestoGP model

This method return the fitted linear model (of class cv.glmnet or
cv.ncvreg) for a PrestoGP model.

## Usage

``` r
# S4 method for class 'PrestoGPModel'
get_linear_model(model)
```

## Arguments

- model:

  The PrestoGP model object

## Value

The fitted linear model (of class cv.glmnet or cv.ncvreg).

## Details

It is important to note that the model is fit to the transformed data.
The CV error rate and predicted values of Y will not be correct for the
original (untransformed) data. This method should be used primarily for
examining the coefficient paths and generating plots.

## References

- Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
  land-use regression for ground-level nitrogen dioxide", The Annals of
  Applied Statistics (2021) 15(2):688-710.

## See also

[`PrestoGPModel-class`](https://niehs.github.io/PrestoGP/reference/PrestoGPModel-class.md),
[`prestogp_fit`](https://niehs.github.io/PrestoGP/reference/prestogp_fit-PrestoGPModel-method.md),
[`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html),
[`cv.ncvreg`](https://pbreheny.github.io/ncvreg/reference/cv.ncvreg.html)

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
#> Current penalized negative log likelihood: 487.5444 
#> Current MSE: 9.104869 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 483.4359 
#> Current MSE: 9.056759 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 482.989 
#> Current MSE: 9.074136 
#> Beginning iteration 4 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 4 complete 
#> Current penalized negative log likelihood: 482.989 
#> Current MSE: 9.061437 
get_linear_model(soil.vm)
#> 
#> Call:  cv.glmnet(x = as.matrix(model@X_tilde), y = as.matrix(model@y_tilde),      nfolds = nfolds, foldid = foldid, parallel = parallel, relax = penalty ==          "relaxed", alpha = model@alpha, family = family, penalty.factor = pen.factor) 
#> 
#> Measure: Mean-Squared Error 
#> 
#>      Lambda Index Measure      SE Nonzero
#> min 0.00634    26  0.4752 0.05124       4
#> 1se 0.06484     1  0.4888 0.05022       0
```
