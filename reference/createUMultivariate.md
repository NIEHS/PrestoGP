# Create the sparse triangular matrix U for multivariate Vecchia models

This creates the sparse triangular matrix U for multivariate Vecchia
models. This matrix can be used to estimate the likelihood or transform
the data to be iid. This function is a multivariate version of
[`createU`](https://rdrr.io/pkg/GPvecchia/man/createU.html).

## Usage

``` r
createUMultivariate(vec.approx, params, cov_func = NULL)
```

## Arguments

- vec.approx:

  Object returned by
  [`vecchia_Mspecify`](https://niehs.github.io/PrestoGP/reference/vecchia_Mspecify.md).

- params:

  Vector of covariance parameters. See
  [`create_param_sequence`](https://niehs.github.io/PrestoGP/reference/create_param_sequence.md)
  or the examples below for details about the format of this vector.

- cov_func:

  The function used to compute the covariance between two observations.
  Defaults to a Matern model.

## Value

A list containing the sparse upper trianguler U, plus additional objects
required for other functions.

## Details

This function will be much slower if a non-default cov_func is
specified. More importantly, there is no guarantee that the resulting
covariance matrices will be positive definite. We recommend using the
default (Matern) covariance function unless you know exactly what you
are doing. See Apanasovich et al. (2012) for a description of how the
cross-covariances are computed.

## References

- Apanasovich, T.V., Genton, M.G. and Sun, Y. "A valid Matérn class of
  cross-covariance functions for multivariate random fields with any
  number of components", Journal of the American Statistical
  Association (2012) 107(497):180-193.

- Katzfuss, M., and Guinness, J. "A general framework for Vecchia
  approximations of Gaussian processes", Statistical Science (2021)
  36(1):124-141.

## See also

[`createU`](https://rdrr.io/pkg/GPvecchia/man/createU.html),
[`vecchia_Mspecify`](https://niehs.github.io/PrestoGP/reference/vecchia_Mspecify.md),
[`create_param_sequence`](https://niehs.github.io/PrestoGP/reference/create_param_sequence.md)

## Examples

``` r
data(soil)
soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
locs <- as.matrix(soil[,1:2])
locsm <- list()
locsm[[1]] <- locsm[[2]] <- locs
soil.va <- vecchia_Mspecify(locsm, m=10)

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

soil.u <- createUMultivariate(soil.va, params)
```
