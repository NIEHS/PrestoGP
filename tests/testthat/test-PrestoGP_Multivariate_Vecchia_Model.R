test_that("Invalid locs input", {
  model <- new("MultivariateVecchiaModel")
  expect_error(prestogp_fit(model, list(as.matrix(1:3)), list(as.matrix(1:3)),
                            as.matrix(1:3)),
               "locs must be a list for multivariate models")
})

test_that("Invalid Y input", {
  model <- new("MultivariateVecchiaModel")
  expect_error(prestogp_fit(model, as.matrix(1:3), list(as.matrix(1:3)),
                            list(as.matrix(1:3))),
               "Y must be a list for multivariate models")
})

test_that("Invalid X input", {
  model <- new("MultivariateVecchiaModel")
  expect_error(prestogp_fit(model, list(as.matrix(1:3)), as.matrix(1:3),
                            list(as.matrix(1:3))),
               "X must be a list for multivariate models")
})

test_that("locs/Y length mismatch", {
  model <- new("MultivariateVecchiaModel")
  expect_error(prestogp_fit(model, list(1:3, 2:4), list(as.matrix(1:3)),
                            list(as.matrix(1:3))),
               "locs and Y must have the same length")
})

test_that("locs/X length mismatch", {
  model <- new("MultivariateVecchiaModel")
  expect_error(prestogp_fit(model, list(1:3), list(as.matrix(1:3),
                                                   as.matrix(2:4)),
                            list(as.matrix(1:3))),
               "locs and X must have the same length")
})

test_that("locs not a matrix", {
  model <- new("MultivariateVecchiaModel")
  expect_error(prestogp_fit(model, list(1:3), list(as.matrix(1:3)),
                            list(1:3)),
               "Each locs must be a matrix")
})

test_that("X not a matrix", {
  model <- new("MultivariateVecchiaModel")
  expect_error(prestogp_fit(model, list(1:3), list(1:3),
                            list(as.matrix(1:3))),
               "Each X must be a matrix")
})

test_that("Y not a numeric matrix/vector", {
  model <- new("MultivariateVecchiaModel")
  expect_error(prestogp_fit(model, list("foo"), list(as.matrix(1:3)),
                            list(as.matrix(1:3))),
               "Each Y must be a numeric vector or matrix")
})

test_that("locs with differing numbers of columns", {
  model <- new("MultivariateVecchiaModel")
  expect_error(prestogp_fit(model, list(1:4, 1:4),
                            list(as.matrix(1:4), as.matrix(1:4)),
                            list(as.matrix(1:4), matrix(1:4, nrow=2))),
               "All locs must have the same number of columns")
})

test_that("Y has multiple columns", {
  model <- new("MultivariateVecchiaModel")
  expect_error(prestogp_fit(model, list(matrix(1:4, nrow=2)),
                            list(as.matrix(1:4)), list(as.matrix(1:4))),
               "Each Y must have only 1 column")
})

test_that("nrow(Y) != nrow(locs)", {
  model <- new("MultivariateVecchiaModel")
  expect_error(prestogp_fit(model, list(as.matrix(1:4)),
                            list(as.matrix(1:4)), list(as.matrix(1:3))),
               "Each Y must have the same number of rows as locs")
})

test_that("nrow(Y) != nrow(X)", {
  model <- new("MultivariateVecchiaModel")
  expect_error(prestogp_fit(model, list(as.matrix(1:4)),
                            list(as.matrix(1:3)), list(as.matrix(1:4))),
               "Each Y must have the same number of rows as X")
})

test_that("Simulated dataset multivariate spatial", {
    load("sim_multivariate_big.RData")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors=25)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list, X.st, locs.list,
        scaling = c(1, 1), apanasovich = TRUE, verbose = FALSE,
        optim.control = list(
            trace = 0, maxit = 5000,
            reltol = 1e-3
        )
    )
    beta.out <- as.vector(pgp.mmodel1@beta)
    params.out <- pgp.mmodel1@covparams

    expect_length(beta.out, 31)
    expect_length(params.out, 15)
    expect_equal(beta.out, c(0, 0.96, 0.93, 0.92, 0.89, rep(0, 2), 0.03,
                             rep(0,3), 0.57, 0.72, 1.11, 1, 0, 0.06, rep(0,4),
                             0.93, 0.87, 1.03, 0.92, 0.05, rep(0, 2), 0.01,
                             0.05, 0), tolerance=0.05)
    expect_equal(params.out[1], 1.7, tolerance=2.5)
    expect_equal(params.out[2], 2.9, tolerance=3.3)
    expect_equal(params.out[3], 3.2, tolerance=5)
    expect_equal(params.out[4], 0.71, tolerance=10)
    expect_equal(params.out[5], 0.56, tolerance=2.8)
    expect_equal(params.out[6], 0.4, tolerance=1.7)
    expect_equal(params.out[7], 0.63, tolerance=1.3)
    expect_equal(params.out[8], 0.47, tolerance=1.2)
    expect_equal(params.out[9], 0.74, tolerance=1.1)
    expect_equal(params.out[10], 1.8, tolerance=1.4)
    expect_equal(params.out[11], 2.6, tolerance=2.4)
    expect_equal(params.out[12], 1.4, tolerance=0.8)
    expect_equal(params.out[13], 0.15, tolerance=0.6)
    expect_equal(params.out[14], 0.32, tolerance=0.4)
    expect_equal(params.out[15], 0.17, tolerance=0.3)
})

test_that("Simulated dataset multivariate spatiotemporal", {
    source("sim_multivariate_big_st.R")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 25)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list, X.st, locs.list,
        scaling = c(1, 1, 2), verbose = FALSE,
        optim.control = list(
            trace = 0, maxit = 5000,
            reltol = 1e-3
        )
    )
    beta.out <- as.vector(pgp.mmodel1@beta)
    params.out <- pgp.mmodel1@covparams

    expect_length(beta.out, 31)
    expect_length(params.out, 18)
    expect_equal(beta.out, c(
        0, 0.91, 0.86, 0.82, 0.97, rep(0, 6), 0.95,
        0.97, 0.92, 0.78, rep(0, 6), 0.8, 0.97,
        1.04, 0.81, rep(0, 6)
    ), tolerance = 1.1)
})

test_that("Invalid locs input for prediction", {
    source("sim_multivariate_small_pred.R")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors=5)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
                                locs.list.otr, scaling=c(1,1),
                                apanasovich=TRUE, verbose=FALSE,
                                optim.control=list(trace=0,maxit=5000,
                                                   reltol=1e-3))
    expect_error(prestogp_predict(pgp.mmodel1, X.st.otst, as.matrix(1:3)),
                 "locs must be a list for multivariate models")
})

test_that("Invalid X input for prediction", {
    source("sim_multivariate_small_pred.R")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors=5)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
                                locs.list.otr, scaling=c(1,1),
                                apanasovich=TRUE, verbose=FALSE,
                                optim.control=list(trace=0,maxit=5000,
                                                   reltol=1e-3))
    expect_error(prestogp_predict(pgp.mmodel1, as.matrix(1:3), locs.list.otst),
                 "X must be a list for multivariate models")
})

test_that("locs/X length mismatch for prediction", {
    source("sim_multivariate_small_pred.R")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors=5)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
                                locs.list.otr, scaling=c(1,1),
                                apanasovich=TRUE, verbose=FALSE,
                                optim.control=list(trace=0,maxit=5000,
                                                   reltol=1e-3))
    expect_error(prestogp_predict(pgp.mmodel1, X.st.otst, list(as.matrix(1:3))),
                 "locs and X must have the same length")
})

test_that("locs/locs_train length mismatch", {
    source("sim_multivariate_small_pred.R")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors=5)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
                                locs.list.otr, scaling=c(1,1),
                                apanasovich=TRUE, verbose=FALSE,
                                optim.control=list(trace=0,maxit=5000,
                                                   reltol=1e-3))
    expect_error(prestogp_predict(pgp.mmodel1, list(1:3, 1:3, 1:3),
                                  list(1:3, 1:3, 1:3)),
                 "Training and test set locs must have the same length")
})

test_that("locs not a matrix for prediction", {
    source("sim_multivariate_small_pred.R")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors=5)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
                                locs.list.otr, scaling=c(1,1),
                                apanasovich=TRUE, verbose=FALSE,
                                optim.control=list(trace=0,maxit=5000,
                                                   reltol=1e-3))
    expect_error(prestogp_predict(pgp.mmodel1, list(as.matrix(1:3),
                                                    as.matrix(1:3)),
                                  list(1:3, as.matrix(1:3))),
                 "Each locs must be a matrix")
})

test_that("ncol(locs) != ncol(locs_train)", {
    source("sim_multivariate_small_pred.R")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors=5)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
                                locs.list.otr, scaling=c(1,1),
                                apanasovich=TRUE, verbose=FALSE,
                                optim.control=list(trace=0,maxit=5000,
                                                   reltol=1e-3))
    expect_error(prestogp_predict(pgp.mmodel1, list(as.matrix(1:3),
                                                    as.matrix(1:3)),
                                  list(as.matrix(1:3), as.matrix(1:3))),
                 "All locs must have the same number of columns as locs_train")
})

test_that("X not a matrix for prediction", {
    source("sim_multivariate_small_pred.R")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors=5)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
                                locs.list.otr, scaling=c(1,1),
                                apanasovich=TRUE, verbose=FALSE,
                                optim.control=list(trace=0,maxit=5000,
                                                   reltol=1e-3))
    X.st.otst[[1]] <- as.vector(X.st.otst[[1]])
    expect_error(prestogp_predict(pgp.mmodel1, X.st.otst, locs.list.otst),
                 "Each X must be a matrix")
})

test_that("nrow(X) != nrow(locs) for prediction", {
    source("sim_multivariate_small_pred.R")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors=5)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
                                locs.list.otr, scaling=c(1,1),
                                apanasovich=TRUE, verbose=FALSE,
                                optim.control=list(trace=0,maxit=5000,
                                                   reltol=1e-3))
    locs.list.otst[[1]] <- rbind(locs.list.otst[[1]], rep(0, 2))
    expect_error(prestogp_predict(pgp.mmodel1, X.st.otst, locs.list.otst),
                 "Each X must have the same number of rows as locs")
})

test_that("ncol(X) != ncol(X_train) for prediction", {
    source("sim_multivariate_small_pred.R")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors=5)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
                                locs.list.otr, scaling=c(1,1),
                                apanasovich=TRUE, verbose=FALSE,
                                optim.control=list(trace=0,maxit=5000,
                                                   reltol=1e-3))
    X.st.otst[[1]] <- cbind(X.st.otst[[1]], rep(0, nrow(X.st.otst[[1]])))
    expect_error(prestogp_predict(pgp.mmodel1, X.st.otst, locs.list.otst),
                 "X and X_train must have the same number of predictors")
})

test_that("m too small for prediction", {
    source("sim_multivariate_small_pred.R")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors=5)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
                                locs.list.otr, scaling=c(1,1),
                                apanasovich=TRUE, verbose=FALSE,
                                optim.control=list(trace=0,maxit=5000,
                                                   reltol=1e-3))
    expect_error(prestogp_predict(pgp.mmodel1, X.st.otst, locs.list.otst, m=2),
                 "m must be at least 3")
})

test_that("m too large for prediction", {
    source("sim_multivariate_small_pred.R")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors=5)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
                                locs.list.otr, scaling=c(1,1),
                                apanasovich=TRUE, verbose=FALSE,
                                optim.control=list(trace=0,maxit=5000,
                                                   reltol=1e-3))
    expect_warning(prestogp_predict(pgp.mmodel1, X.st.otst, locs.list.otst,
                                    m=51),
                   "Conditioning set size m chosen to be >=n. Changing to m=n-1")
})

test_that("Simulated dataset multivariate spatial prediction", {
    source("sim_multivariate_big_pred.R")
    pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 25)
    pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
        locs.list.otr,
        scaling = c(1, 1),
        apanasovich = TRUE, verbose = FALSE,
        optim.control = list(
            trace = 0, maxit = 5000,
            reltol = 1e-3
        )
    )

    pgp.mmodel1.pred <- prestogp_predict(pgp.mmodel1, X.st.otst, locs.list.otst)

    mse <- mean((pgp.mmodel1.pred$means - unlist(y.list.otst))^2)
    me <- mean(pgp.mmodel1.pred$means - unlist(y.list.otst))

    expect_equal(mse, 1.99, tolerance = 0.1)
    expect_equal(me, -0.04, tolerance = 0.03)
})
