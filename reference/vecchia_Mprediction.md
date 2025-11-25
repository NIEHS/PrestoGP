# Multivariate Vecchia prediction

This function is used to make predictions based on multivariate Vecchia
models. It is a multivariate version of
[`vecchia_prediction`](https://rdrr.io/pkg/GPvecchia/man/vecchia_prediction.html).

## Usage

``` r
vecchia_Mprediction(
  z,
  vecchia.approx,
  covparms,
  var.exact = NULL,
  return.values = "mean"
)
```

## Arguments

- z:

  The observed data.

- vecchia.approx:

  A Vecchia object returned by
  [`vecchia_Mspecify`](https://niehs.github.io/PrestoGP/reference/vecchia_Mspecify.md).

- covparms:

  Vector of covariance parameters. See
  [`create_param_sequence`](https://niehs.github.io/PrestoGP/reference/create_param_sequence.md)
  or the examples below for details about the format of this vector.

- var.exact:

  Should prediction variances by computed exactly, or is a (faster)
  approximation acceptable? See
  [`vecchia_prediction`](https://rdrr.io/pkg/GPvecchia/man/vecchia_prediction.html).

- return.values:

  Values that should be returned. Possible values include "mean",
  "meanvar", "meanmat", and "all". See
  [`vecchia_prediction`](https://rdrr.io/pkg/GPvecchia/man/vecchia_prediction.html).
  Defaults to "mean".

## Value

The posterior means/variances/V matrices at the observed and unobserved
locations. See
[`vecchia_prediction`](https://rdrr.io/pkg/GPvecchia/man/vecchia_prediction.html).

## References

- Katzfuss, M., and Guinness, J. "A general framework for Vecchia
  approximations of Gaussian processes", Statistical Science (2021)
  36(1):124-141.

- Katzfuss, M., Guinness, J., Gong, W. and Zilber, D. "Vecchia
  approximations of Gaussian-process predictions", Journal of
  Agricultural, Biological and Environmental Statistics (2020)
  25:383-414.

## See also

[`vecchia_prediction`](https://rdrr.io/pkg/GPvecchia/man/vecchia_prediction.html),
[`vecchia_Mspecify`](https://niehs.github.io/PrestoGP/reference/vecchia_Mspecify.md),
[`create_param_sequence`](https://niehs.github.io/PrestoGP/reference/create_param_sequence.md)

## Examples

``` r
data(soil)
soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
locs <- as.matrix(soil[,1:2])
locsm <- list()
locsm[[1]] <- locsm[[2]] <- locs
locsp <- locsm
locsp[[1]] <- locsp[[1]] + 0.5
locsp[[2]] <- locsp[[2]] - 0.5
soil.vap <- vecchia_Mspecify(locsm, m=10, locs.list.pred=locsp)

pseq <- create_param_sequence(2)
# Initialize the vector of covariance parameters
params <- rep(NA, pseq[5,2])
# Sigma parameters:
params[pseq[1,1]:pseq[1,2]] <- c(100, 80)
# Scale parameters:
params[pseq[2,1]:pseq[2,2]] <- c(60, 50)
# Smoothness parameters:
params[pseq[3,1]:pseq[3,2]] <- c(0.5, 0.5)
# Nuggets:
params[pseq[4,1]:pseq[4,2]] <- c(30, 30)
# Correlation:
params[pseq[5,1]:pseq[5,2]] <- -0.9

soil.yhat <- vecchia_Mprediction(rnorm(nrow(locs)), soil.vap, params)
```
