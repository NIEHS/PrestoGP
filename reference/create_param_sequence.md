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
# Columns 1+2 are location coordinates; column 3 is elevation
locs2 <- as.matrix(soil250[,1:3])

soil.vm2 <- new("VecchiaModel", n_neighbors = 10)
# Fit separate scale parameters for location and elevation
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
#> Current penalized negative log likelihood: -259.1528 
#> Current MSE: 0.008211436 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: -259.1528 
#> Current MSE: 0.009062582 

pseq <- create_param_sequence(1, 2)
soil2.params <- soil.vm2@covparams
# Extract the sigmas
soil2.params[pseq[1,1]:pseq[1,2]]
#> [1] 0.005046449
# Extract the scale parameters
soil2.params[pseq[2,1]:pseq[2,2]]
#> [1] 3.14302152 0.02809024
# Extract the smoothness parameters
soil2.params[pseq[3,1]:pseq[3,2]]
#> [1] 1.580789
# Extract the nuggets
soil2.params[pseq[4,1]:pseq[4,2]]
#> [1] 0.003012928

# Multivariate model with user-specified initial Matern parameter
# estimates
ym <- list()
ym[[1]] <- soil250[,4] # predict sand/silt portion of the sample
ym[[2]] <- soil250[,5]
ym[[3]] <- soil250[,6]
Xm <- list()
Xm[[1]] <- Xm[[2]] <- Xm[[3]] <- as.matrix(soil250[,7:22])
locsm <- list()
locsm[[1]] <- locsm[[2]] <- locsm[[3]] <- as.matrix(soil250[,1:3])

# Initialize the vector of initial Matern parameters estimates
pseq2 <- create_param_sequence(3, 2)
soil.params0 <- rep(NA, pseq2[5, 2])

# Specify the initial sigma estimates
soil.params0[pseq2[1, 1]:pseq2[1, 2]] <- c(1, 5, 8)
# Scale parameters
scale.seq <- pseq2[2,1]:pseq2[2,2]
# Specify the scale parameter for location, outcome 1
soil.params0[scale.seq[1]] <- 12.8
# Specify the scale parameter for elevation, outcome 1
soil.params0[scale.seq[2]] <- 12.8
# Specify the scale parameter for location, outcome 2
soil.params0[scale.seq[3]] <- 21.5
# Specify the scale parameter for elevation, outcome 2
soil.params0[scale.seq[4]] <- 21.5
# Specify the scale parameter for location, outcome 3
soil.params0[scale.seq[5]] <- 17.8
# Specify the scale parameter for elevation, outcome 3
soil.params0[scale.seq[6]] <- 17.8
# Specify the initial smoothness parameter estimates
soil.params0[pseq2[3, 1]:pseq2[3, 2]] <- c(0.5, 0.5, 0.5)
# Specify the initial nugget estimates
soil.params0[pseq2[4, 1]:pseq2[4, 2]] <- c(0.25, 0.5, 0.5)
# Specify the initial correlation estimates
soil.params0[pseq2[5, 1]:pseq2[5, 2]] <- c(0, 0, 0)

soil.mvm <-  new("MultivariateVecchiaModel", n_neighbors = 25)
soil.mvm <- prestogp_fit(soil.mvm, ym, Xm, locsm, scaling= c(1, 1, 2),
covparams = soil.params0, Y.names = names(soil250)[4:6])
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
#> Current penalized negative log likelihood: 1174.223 
#> Current MSE: 2.020575 
#> Beginning iteration 2 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 2 complete 
#> Current penalized negative log likelihood: 1137.111 
#> Current MSE: 2.767673 
#> Beginning iteration 3 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 3 complete 
#> Current penalized negative log likelihood: 1136.835 
#> Current MSE: 2.829191 
#> Beginning iteration 4 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 4 complete 
#> Current penalized negative log likelihood: 1136.452 
#> Current MSE: 2.882692 
#> Beginning iteration 5 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 5 complete 
#> Current penalized negative log likelihood: 1135.68 
#> Current MSE: 2.870947 
#> Beginning iteration 6 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 6 complete 
#> Current penalized negative log likelihood: 1135.045 
#> Current MSE: 2.867495 
#> Beginning iteration 7 
#> Estimating theta... 
#> Estimation of theta complete 
#> Estimating beta... 
#> Estimation of beta complete 
#> Iteration 7 complete 
#> Current penalized negative log likelihood: 1135.045 
#> Current MSE: 2.870512 

# Extract the estimated Matern paramaters
soil.params <- soil.mvm@covparams
# Extract the sigmas
soil.params[pseq2[1,1]:pseq2[1,2]]
#> [1] 0.4038428 1.9177732 4.7231172
# Extract the scale parameter for location, outcome 1
soil.params[scale.seq[1]]
#> [1] 34.20424
# Extract the scale parameter for elevation, outcome 1
soil.params[scale.seq[2]]
#> [1] 18.93287
# Extract the scale parameter for location, outcome 2
soil.params[scale.seq[3]]
#> [1] 44.7757
# Extract the scale parameter for elevation, outcome 2
soil.params[scale.seq[4]]
#> [1] 140.7954
# Extract the scale parameter for location, outcome 3
soil.params[scale.seq[5]]
#> [1] 21.36912
# Extract the scale parameter for elevation, outcome 3
soil.params[scale.seq[6]]
#> [1] 15.78034
# Extract the smoothness parameters
soil.params[pseq2[3,1]:pseq2[3,2]]
#> [1] 0.3725119 0.5690135 0.5023576
# Extract the nuggets
soil.params[pseq2[4,1]:pseq2[4,2]]
#> [1] 0.2432345 1.3513640 0.9633006
# Extract the correlation parameters
soil.corr <- diag(2) / 2
soil.corr[upper.tri(soil.corr)] <- soil.params[pseq2[5,1]:pseq2[5,2]]
#> Warning: number of items to replace is not a multiple of replacement length
soil.corr <- soil.corr + t(soil.corr)
```
