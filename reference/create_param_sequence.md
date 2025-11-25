# Extract specific Matern parameters from a parameter sequence

This function is used to obtain specific Matern parameters (e.g., range
or smoothness) from the covparams slot of a PrestoGPModel object.

## Usage

``` r
create_param_sequence(P, ns = 1)
```

## Arguments

- P:

  Number of outcome variables

- ns:

  Number of scale parameters

## Value

A matrix with five rows and two columns as described below:

- Row 1::

  Starting and ending indices for the sigma parameter(s)

- Row 2::

  Starting and ending indices for the scale parameter(s)

- Row 3::

  Starting and ending indices for the smoothness parameter(s)

- Row 4::

  Starting and ending indices for the nugget(s)

- Row 5::

  Starting and ending indices for the correlation parameter(s)

## Details

This function is intended for advanced users who want to specify the
input Matern parameters for functions such as
[`vecchia_Mlikelihood`](https://niehs.github.io/PrestoGP/reference/vecchia_Mlikelihood.md)
or
[`createUMultivariate`](https://niehs.github.io/PrestoGP/reference/createUMultivariate.md).
To extract the Matern parameters from a fitted PrestoGP model, it is
strongly recommended to use `link{get_theta}` instead.

## References

- Apanasovich, T.V., Genton, M.G. and Sun, Y. "A valid Matérn class of
  cross-covariance functions for multivariate random fields with any
  number of components", Journal of the American Statistical
  Association (2012) 107(497):180-193.

- Genton, M.G. "Classes of kernels for machine learning: a statistics
  perspective", The Journal of Machine Learning Research (2001)
  2:299-312.

## See also

[`PrestoGPModel-class`](https://niehs.github.io/PrestoGP/reference/PrestoGPModel-class.md)

## Examples

``` r
# Space/elevation model
data(soil250, package="geoR")
y2 <- soil250[,7]               # predict pH level
X2 <- as.matrix(soil250[,c(4:6,8:22)])
# columns 1+2 are location coordinates; column 3 is elevation
locs2 <- as.matrix(soil250[,1:3])

soil.vm2 <- new("VecchiaModel", n_neighbors = 10)
# fit separate scale parameters for location and elevation
soil.vm2 <- prestogp_fit(soil.vm2, y2, X2, locs2, scaling = c(1, 1, 2))
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
#> Current penalized negative log likelihood: -258.9123 
#> Current MSE: 0.008209381 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: -258.9123 
#> Current MSE: 0.009980023 

pseq <- create_param_sequence(1, 2)
soil2.params <- soil.vm2@covparams
# sigma
soil2.params[pseq[1,1]:pseq[1,2]]
#> [1] 0.004592186
# scale parameters
soil2.params[pseq[2,1]:pseq[2,2]]
#> [1] 4.09758059 0.03708606
# smoothness parameter
soil2.params[pseq[3,1]:pseq[3,2]]
#> [1] 1.320382
# nugget
soil2.params[pseq[4,1]:pseq[4,2]]
#> [1] 0.003849531

# Multivariate model
ym <- list()
ym[[1]] <- soil250[,4] # predict sand/silt portion of the sample
ym[[2]] <- soil250[,5]
ym[[3]] <- soil250[,6]
Xm <- list()
Xm[[1]] <- Xm[[2]] <- Xm[[3]] <- as.matrix(soil250[,7:22])
locsm <- list()
locsm[[1]] <- locsm[[2]] <- locsm[[3]] <- as.matrix(soil250[,1:3])

soil.mvm <-  new("MultivariateVecchiaModel", n_neighbors = 10)
soil.mvm <- prestogp_fit(soil.mvm, ym, Xm, locsm)
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
#> Current penalized negative log likelihood: 1052.6 
#> Current MSE: 2.033446 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 1006.342 
#> Current MSE: 3.285929 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 1004.776 
#> Current MSE: 3.333343 
#> Beginning iteration 4 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 4 complete 
#> Current penalized negative log likelihood: 1004.776 
#> Current MSE: 3.966228 

pseq <- create_param_sequence(3, 2)
soil.params <- soil.mvm@covparams
# sigmas
soil.params[pseq[1,1]:pseq[1,2]]
#> [1] 0.8121059 4.4694598 7.6915030
# scale parameters
scale.seq <- pseq[2,1]:pseq[2,2]
# scale parameter for location, outcome 1
soil.params[scale.seq[1]]
#> [1] 17.07644
# scale parameter for elevation, outcome 1
soil.params[scale.seq[2]]
#> [1] 19.80937
# scale parameter for location, outcome 2
soil.params[scale.seq[3]]
#> [1] 9.165145
# scale parameter for elevation, outcome 2
soil.params[scale.seq[4]]
#> [1] 0.5874919
# scale parameter for location, outcome 3
soil.params[scale.seq[5]]
#> [1] 0.3783888
# scale parameter for elevation, outcome 3
soil.params[scale.seq[6]]
#> [1] 0.935525
# smoothness parameters
soil.params[pseq[3,1]:pseq[3,2]]
#> [1] 0.1433850 0.1598858 0.5561018
# nuggets
soil.params[pseq[4,1]:pseq[4,2]]
#> [1] -0.2868219 -0.3519231 -0.8842306
# correlation
soil.corr <- diag(2) / 2
soil.corr[upper.tri(soil.corr)] <- soil.params[pseq[5,1]:pseq[5,2]]
#> Warning: number of items to replace is not a multiple of replacement length
soil.corr <- soil.corr + t(soil.corr)
```
