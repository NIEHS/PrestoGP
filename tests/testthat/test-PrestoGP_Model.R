test_that("beta.hat not numeric", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs, beta.hat = "foo"),
    "beta.hat parameter must be a numeric vector"
  )
})

test_that("beta.hat not a vector", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs, beta.hat = matrix(1:3)),
    "beta.hat parameter must be a numeric vector"
  )
})

test_that("beta.hat incorrect length", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs, beta.hat = 1:3),
    "Length of beta.hat must match the number of predictors"
  )
})

test_that("tol not numeric", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs, tol = "foo"),
    "tol must be numeric"
  )
})

test_that("tol not a scalar", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs, tol = 1:2),
    "tol must be a scalar"
  )
})

test_that("tol out of range", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs, tol = 1.1),
    "tol must satisfy 0<tol<=1"
  )
})

test_that("length(scaling) != ncol(locs)", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs, scaling = 1:5),
    "Length of scaling must equal ncol of locs"
  )
})

test_that("scaling not sequential", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs, scaling = c(1, 3)),
    "scaling must consist of sequential integers starting at 1"
  )
})

test_that("Invalid Apanasovich model", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs,
      scaling = 1:2,
      apanasovich = TRUE
    ),
    "Apanasovich models require a common scale parameter"
  )
})

test_that("covparams not numeric", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs, covparams = "foo"),
    "covparams must be a numeric vector"
  )
})

test_that("covparams not numeric", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs, covparams = "foo"),
    "covparams must be a numeric vector"
  )
})

test_that("covparams not a vector", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs, covparams = matrix(1:3)),
    "covparams must be a numeric vector"
  )
})

test_that("covparams incorrect length", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs, covparams = 1:3),
    "Incorrect number of parameters in covparams"
  )
})

test_that("covparams incorrect length - multivariate", {
  source("sim_multivariate_small.R")
  pgp.mmodel1 <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(pgp.mmodel1, y.list, X.st, locs.list,
      covparams = c(
        0.9, 0.9, 0.5, 0.5, 0.5, 0.5,
        0.1, 0.1
      )
    ),
    "Incorrect number of parameters in covparams"
  )
})

test_that("Invalid initial variance estimate", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs,
      covparams = c(-0.9, 0.5, 0.5, 0.1)
    ),
    "Initial variance estimates must be positive"
  )
})

test_that("Invalid initial range estimate", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs,
      covparams = c(0.9, -0.5, 0.5, 0.1)
    ),
    "Initial range estimates must be positive"
  )
})

test_that("Invalid initial nugget estimate", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs,
      covparams = c(0.9, 0.5, 0.5, -0.1)
    ),
    "Initial nugget estimates must be positive"
  )
})

test_that("Invalid initial smoothness estimate", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel")
  expect_error(
    prestogp_fit(pgp.model1, y, X, locs,
      covparams = c(0.9, 0.5, -0.5, 0.1)
    ),
    "Initial smoothness estimates must be between 0 and 2.5"
  )
})

test_that("Invalid initial correlation estimate", {
  source("sim_multivariate_small.R")
  pgp.mmodel1 <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(pgp.mmodel1, y.list, X.st, locs.list,
      covparams = c(
        0.9, 0.9, 0.5, 0.5, 0.5, 0.5,
        0.1, 0.1, -1.2
      )
    ),
    "Initial correlation estimates must be between -1 and 1"
  )
})

test_that("m too small", {
  source("sim_multivariate_small.R")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 2)
  expect_error(
    prestogp_fit(pgp.mmodel1, y.list, X.st, locs.list),
    "m must be at least 3"
  )
})

test_that("m too large", {
  source("sim_vecchia_small.R")
  pgp.model1 <- new("VecchiaModel", n_neighbors = 101)
  expect_warning(
    prestogp_fit(pgp.model1, y, X, locs, quiet = TRUE),
    "Conditioning set size m chosen to be >=n. Changing to m=n-1"
  )
})
