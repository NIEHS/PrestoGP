# Evaluation of the multivariate Vecchia likelihood

This function is used to evaluate the multivariate Vecchia likelihood.

## Usage

``` r
vecchia_Mlikelihood(z, vecchia.approx, covparams)
```

## Arguments

- z:

  The observed data.

- vecchia.approx:

  A Vecchia object returned by
  [`vecchia_Mspecify`](https://niehs.github.io/PrestoGP/reference/vecchia_Mspecify.md).

- covparams:

  Vector of covariance parameters. See
  [`create_param_sequence`](https://niehs.github.io/PrestoGP/reference/create_param_sequence.md)
  or the examples below for details about the format of this vector.

## Value

The log likelihood implied by the multivariate Vecchia approximation.

## References

- Katzfuss, M., and Guinness, J. "A general framework for Vecchia
  approximations of Gaussian processes", Statistical Science (2021)
  36(1):124-141.

## See also

[`vecchia_likelihood`](https://rdrr.io/pkg/GPvecchia/man/vecchia_likelihood.html),
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

vecchia_Mlikelihood(rnorm(nrow(locs)), soil.va, params)
#> [1] NA
```
