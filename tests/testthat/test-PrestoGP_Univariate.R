test_that("Invalid locs input", {
  model <- new("VecchiaModel")
  expect_error(
    prestogp_fit(
      model, as.matrix(1:3), as.matrix(1:3),
      1:3
    ),
    "locs must be a matrix"
  )
})

test_that("Invalid Y input", {
  model <- new("VecchiaModel")
  expect_error(
    prestogp_fit(
      model, "y", as.matrix(1:3),
      as.matrix(1:3)
    ),
    "Y must be a numeric vector or matrix"
  )
})

test_that("Invalid X input", {
  model <- new("VecchiaModel")
  expect_error(
    prestogp_fit(
      model, as.matrix(1:3), 1:3,
      as.matrix(1:3)
    ),
    "X must be a matrix"
  )
})

test_that("Y has multiple columns", {
  model <- new("VecchiaModel")
  expect_error(
    prestogp_fit(
      model, matrix(1:4, nrow = 2), as.matrix(1:4),
      as.matrix(1:4)
    ),
    "Y must have only 1 column"
  )
})

test_that("nrow(Y) != nrow(locs)", {
  model <- new("VecchiaModel")
  expect_error(
    prestogp_fit(
      model, as.matrix(1:4), as.matrix(1:4),
      as.matrix(1:3)
    ),
    "Y must have the same number of rows as locs"
  )
})

test_that("nrow(Y) != nrow(X)", {
  model <- new("VecchiaModel")
  expect_error(
    prestogp_fit(
      model, as.matrix(1:4), as.matrix(1:3),
      as.matrix(1:4)
    ),
    "Y must have the same number of rows as X"
  )
})

test_that("Simulated dataset spatial", {
  load("sim_vecchia.RData")
  pgp.model1 <- new("VecchiaModel", n_neighbors = 25)
  pgp.model1 <- prestogp_fit(pgp.model1, y, X, locs,
    scaling = c(1, 1), apanasovich = TRUE, verbose = FALSE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  beta.out <- as.vector(pgp.model1@beta)
  params.out <- pgp.model1@covparams

  expect_length(beta.out, 11)
  expect_length(params.out, 5)
  expect_equal(beta.out, c(0.01, 0.86, 0.98, 0.94, 0.9, rep(0, 6)),
    tolerance = 0.03
  )
  expect_equal(params.out[1], 1.6, tolerance = 0.5)
  expect_equal(params.out[2], 0.4, tolerance = 0.2)
  expect_equal(params.out[3], 0.59, tolerance = 0.2)
  expect_equal(params.out[4], 2.0, tolerance = 0.15)

  pgp.model2 <- new("FullModel")
  pgp.model2 <- prestogp_fit(pgp.model2, y, X, locs,
    scaling = c(1, 1), apanasovich = TRUE, verbose = FALSE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  beta.out2 <- as.vector(pgp.model2@beta)
  params.out2 <- pgp.model2@covparams

  expect_length(beta.out2, 11)
  expect_length(params.out2, 5)
  expect_equal(beta.out2, c(0.01, 0.86, 0.98, 0.95, 0.9, rep(0, 6)),
    tolerance = 0.03
  )
  expect_equal(params.out2[1], 1.5, tolerance = 0.6)
  expect_equal(params.out2[2], 0.4, tolerance = 0.15)
  expect_equal(params.out2[3], 0.62, tolerance = 0.2)
  expect_equal(params.out2[4], 2.0, tolerance = 0.15)

  # Vecchia and full models should be approximately equal
  expect_equal(beta.out[1], beta.out2[1], tolerance = 0.07)
  expect_equal(beta.out[-1], beta.out2[-1], tolerance = 0.04)
  expect_equal(params.out[1], params.out2[1], tolerance = 1)
  expect_equal(params.out[2] - params.out2[2], 0, tolerance = 0.2)
  expect_equal(params.out[3], params.out2[3], tolerance = 0.3)
  expect_equal(params.out[4], params.out2[4], tolerance = 0.2)
})

test_that("Simulated dataset spatiotemporal", {
  source("sim_vecchia_st.R")
  pgp.model1 <- new("VecchiaModel", n_neighbors = 25)
  pgp.model1 <- prestogp_fit(pgp.model1, y, X, locs,
    scaling = c(1, 1, 2),
    apanasovich = FALSE, verbose = FALSE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  beta.out <- as.vector(pgp.model1@beta)
  params.out <- pgp.model1@covparams

  expect_length(beta.out, 11)
  expect_length(params.out, 6)
  expect_equal(beta.out, c(0.01, 0.93, 1.01, 0.91, 0.99, rep(0, 6)),
    tolerance = 0.015
  )
  expect_equal(params.out[1], 1.7, tolerance = 0.55)
  expect_equal(params.out[2] - 0.19, 0, tolerance = 0.05)
  expect_equal(params.out[3] - 0.22, 0, tolerance = 0.05)
  expect_equal(params.out[4], 1.12, tolerance = 0.3)
  expect_equal(params.out[5], 0.62, tolerance = 0.05)

  pgp.model2 <- new("FullModel", n_neighbors = 25)
  pgp.model2 <- prestogp_fit(pgp.model2, y, X, locs,
    scaling = c(1, 1, 2),
    apanasovich = FALSE, verbose = FALSE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  beta.out2 <- as.vector(pgp.model2@beta)
  params.out2 <- pgp.model2@covparams

  expect_length(beta.out2, 11)
  expect_length(params.out2, 6)
  expect_equal(beta.out2, c(-0.03, 0.93, 1, 0.91, 0.98, rep(0, 6)),
    tolerance = 0.02
  )
  expect_equal(params.out2[1], 1.6, tolerance = 0.5)
  expect_equal(params.out2[2] - 0.19, 0, tolerance = 0.05)
  expect_equal(params.out2[3] - 0.22, 0, tolerance = 0.05)
  expect_equal(params.out2[4], 1.18, tolerance = 0.15)
  expect_equal(params.out2[5], 0.64, tolerance = 0.05)

  # Vecchia and full models should be approximately equal
  expect_equal(beta.out[1], beta.out2[1], tolerance = 0.08)
  expect_equal(beta.out[-1], beta.out2[-1], tolerance = 0.03)
  expect_equal(params.out[1], params.out2[1], tolerance = 0.6)
  expect_equal(params.out[2] - params.out2[2], 0, tolerance = 0.06)
  expect_equal(params.out[3] - params.out2[3], 0, tolerance = 0.06)
  expect_equal(params.out[4], params.out2[4], tolerance = 0.3)
  expect_equal(params.out[5], params.out2[5], tolerance = 0.1)
})

test_that("Invalid locs input for prediction", {
  source("sim_vecchia_small_pred.R")
  pgp.model1 <- new("VecchiaModel", n_neighbors = 5)
  pgp.model1 <- prestogp_fit(pgp.model1, y.otr, X.otr,
    locs.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, verbose = FALSE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  expect_error(
    prestogp_predict(pgp.model1, X.otst, 1:3),
    "locs must be a matrix"
  )
})

test_that("Invalid X input for prediction", {
  source("sim_vecchia_small_pred.R")
  pgp.model1 <- new("VecchiaModel", n_neighbors = 5)
  pgp.model1 <- prestogp_fit(pgp.model1, y.otr, X.otr,
    locs.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, verbose = FALSE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  expect_error(
    prestogp_predict(pgp.model1, 1:3, locs.otst),
    "X must be a matrix"
  )
})

test_that("ncol(locs) != ncol(locs_train)", {
  source("sim_vecchia_small_pred.R")
  pgp.model1 <- new("VecchiaModel", n_neighbors = 5)
  pgp.model1 <- prestogp_fit(pgp.model1, y.otr, X.otr,
    locs.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, verbose = FALSE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  locs.otst <- cbind(locs.otst, rep(1, nrow(locs.otst)))
  expect_error(
    prestogp_predict(pgp.model1, X.otst, locs.otst),
    "locs must have the same number of columns as locs_train"
  )
})

test_that("nrow(X) != nrow(locs) for prediction", {
  source("sim_vecchia_small_pred.R")
  pgp.model1 <- new("VecchiaModel", n_neighbors = 5)
  pgp.model1 <- prestogp_fit(pgp.model1, y.otr, X.otr,
    locs.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, verbose = FALSE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  X.otst <- rbind(X.otst, rep(1, ncol(X.otst)))
  expect_error(
    prestogp_predict(pgp.model1, X.otst, locs.otst),
    "X must have the same number of rows as locs"
  )
})

test_that("ncol(X) != ncol(X_train) for prediction", {
  source("sim_vecchia_small_pred.R")
  pgp.model1 <- new("VecchiaModel", n_neighbors = 5)
  pgp.model1 <- prestogp_fit(pgp.model1, y.otr, X.otr,
    locs.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, verbose = FALSE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  X.otst <- cbind(X.otst, rep(1, nrow(X.otst)))
  expect_error(
    prestogp_predict(pgp.model1, X.otst, locs.otst),
    "X and X_train must have the same number of predictors"
  )
})

test_that("m too small for prediction", {
  source("sim_vecchia_small_pred.R")
  pgp.model1 <- new("VecchiaModel", n_neighbors = 5)
  pgp.model1 <- prestogp_fit(pgp.model1, y.otr, X.otr,
    locs.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, verbose = FALSE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  expect_error(
    prestogp_predict(pgp.model1, X.otst, locs.otst, m = 2),
    "m must be at least 3"
  )
})

test_that("m too large for prediction", {
  source("sim_vecchia_small_pred.R")
  pgp.model1 <- new("VecchiaModel", n_neighbors = 5)
  pgp.model1 <- prestogp_fit(pgp.model1, y.otr, X.otr,
    locs.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, verbose = FALSE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  expect_warning(
    prestogp_predict(pgp.model1, X.otst, locs.otst, m = 201),
    "Conditioning set size m chosen to be >=n. Changing to m=n-1"
  )
})

test_that("Simulated spatial prediction", {
  source("sim_vecchia_pred.R")
  pgp.model1 <- new("VecchiaModel", n_neighbors = 25)
  pgp.model1 <- prestogp_fit(pgp.model1, y.otr, X.otr, locs.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, verbose = FALSE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )

  pgp.model1.pred <- prestogp_predict(pgp.model1, X.otst, locs.otst)

  mse <- mean((pgp.model1.pred$means - y.otst)^2)
  me <- mean(pgp.model1.pred$means - y.otst)

  expect_equal(mse, 2.24, tolerance = 0.07)
  expect_equal(me, 0.03, tolerance = 0.05)
})
