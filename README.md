# **PrestoGP**: **P**enalized **Re**gression on **S**patio-**T**emporal **O**utcomes using **G**aussian **P**rocess <a href="https://niehs.github.io/PrestoGP/index.html/"><img src="man/figures/prestoGP-logo.png" align="right" height="139" alt="penalized regression and Gaussian Process plots" /></a>


[![R-CMD-check](https://github.com/NIEHS/PrestoGP/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/NIEHS/PrestoGP/actions/workflows/check-standard.yaml)
[![test-coverage-local](https://github.com/NIEHS/PrestoGP/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/NIEHS/PrestoGP/actions/workflows/test-coverage.yaml)
[![cov](https://NIEHS.github.io/PrestoGP/badges/coverage.svg)](https://github.com/NIEHS/PrestoGP/actions)
[![lint](https://github.com/NIEHS/PrestoGP/actions/workflows/lint.yaml/badge.svg)](https://github.com/NIEHS/PrestoGP/actions/workflows/lint.yaml)
[![pkgdown](https://github.com/NIEHS/PrestoGP/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/NIEHS/PrestoGP/actions/workflows/pkgdown.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)




# Overview

PrestoGP is an R package for scalable penalized regression on spatio-temporal outcomes using Gaussian processes. The package includes the methods described in the paper "Scalable penalized spatiotemporal land-use regression for ground-level nitrogen dioxide" by Messier and Katzfuss (2020). The package is designed to handle big data and is useful for large-scale geospatial exposure assessment and geophysical modeling. PrestoGP expands the methods in Messier and Katzfuss (2020) to include the following important feature:

1. Multivariate spatiotemporal outcomes using multivariate Matern Gaussian process. The method is described in the paper "A valid Matérn class of cross-covariance functions for multivariate random fields with any number of components" by Apanasovich, Genton, and Sun (2012).

2. Simultaneous variable selection and estimation of the fixed effects (i.e. land-use regression variables) and the spatiotemporal random effects. The method is described in the paper "Scalable penalized spatiotemporal land-use regression for ground-level nitrogen dioxide" by Messier and Katzfuss (2020).

3. Imputation of missing outcome data, including outcome data that is missing due to limit of detection effects. 

4. And as always, scalability of the Gaussian Process via the General Vecchia approximation as described in the paper "A general framework for Vecchia approximations of Gaussian processes" by Katzfuss and Guinness (2021).


# Installation

You can install the development version of PrestoGP from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("NIEHS/PrestoGP")
```

# Usage and Examples

## Data Preprocessing

``` r
data(soil)
soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
y <- soil[,4]                   # predict moisture content
X <- as.matrix(soil[,5:9])
locs <- as.matrix(soil[,1:2]) 
```

## Vecchia Model

``` r
soil.vm <- new("VecchiaModel", n_neighbors = 10)
soil.vm <- prestogp_fit(soil.vm, y, X, locs)
```

## Impute Missing y's

``` r
miss <- sample(1:nrow(soil), 20)
y.miss <- y
y.miss[miss] <- NA
soil.vm2 <- new("VecchiaModel", n_neighbors = 10)
soil.vm2 <- prestogp_fit(soil.vm, y, X, locs, impute.y = TRUE)
```

## Impute y's Missing Due to Limit of Detection

``` r
soil.lod <- quantile(y, 0.1)
y.lod <- y
y.lod[y.lod <= soil.lod] <- NA
soil.vm3 <- new("VecchiaModel", n_neighbors = 10)
soil.vm3 <- prestogp_fit(soil.vm, y, X, locs, impute.y = TRUE, lod = soil.lod)
```

## Full Model

``` r
soil.fm <- new("FullModel")
soil.fm <- prestogp_fit(soil.fm, y, X, locs)
```

## Multivariate Model

``` r
ym <- list()
ym[[1]] <- soil[,5]             # predict two nitrogen concentration levels
ym[[2]] <- soil[,7]
Xm <- list()
Xm[[1]] <- Xm[[2]] <- as.matrix(soil[,c(4,6,8,9)])
locsm <- list()
locsm[[1]] <- locsm[[2]] <- locs

soil.mvm <-  new("MultivariateVecchiaModel", n_neighbors = 10)
soil.mvm <- prestogp_fit(soil.mvm, ym, Xm, locsm)
```

## Space/Elevation Model

``` r
data(soil250, package="geoR")
y2 <- soil250[,7]               # predict pH level
X2 <- as.matrix(soil250[,c(4:6,8:22)])
# columns 1+2 are location coordinates; column 3 is elevation
locs2 <- as.matrix(soil250[,1:3])

soil.vm2 <- new("VecchiaModel", n_neighbors = 10)
# fit separate scale parameters for location and elevation
soil.vm2 <- prestogp_fit(soil.vm2, y2, X2, locs2, scaling = c(1, 1, 2))
```

## Test Set Prediction

``` r
# Create training and test sets
n <- length(y)
otr <- rep(FALSE, n)
otr[sample(1:n, size=floor(n/2))] <- TRUE
otst <- !otr
     
# Fit the model on the training set
soil.vmp <- new("VecchiaModel", n_neighbors = 10)
soil.vmp <- prestogp_fit(soil.vm, y[otr], X[otr,], locs[otr,])
     
# Perform predictions on the test set
soil.yhat <- prestogp_predict(soil.vmp, X[otst,], locs[otst,])
```

# Work In Progress

This package is currently under development. Please report any issues on the [Issues page](https://github.com/NIEHS/PrestoGP/issues)

### Multivariate and Censored Geospatial Models For External Exposomics of Data-Sparse Chemicals: United States Chlorotriazine Groundwater Distributions from 1990-2022

The manuscript is in progress and expected to be submitted soon. Please check back for updates or contact Kyle Messier for more information.

### Authors
- Kyle P Messier<sup>1,2</sup>
- Insang Song<sup>1,2</sup>
- Matt W Wheeler<sup>2</sup>
- Myeongjon Kang<sup>3</sup>
- Matthias Katzfuss<sup>3</sup>
- Ruchir R Shah<sup>4</sup>
- Deepak Mav<sup>4</sup>
- Brian Kidd<sup>4</sup>
- Eric Bair<sup>4</sup>

### Affiliations
1. National Institute of Environmental Health Sciences, Division of the National Toxicology Program, Predictive Toxicology Branch
2. National Institute of Environmental Health Sciences, Division of Intramural Research, Biostatistics and Computational Biology Branch
3. University of Wisconsin, Department of Statistics
4. Sciome, LLC
