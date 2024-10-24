test_that("Invalid locs input", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(1:3)), list(as.matrix(1:3)),
      as.matrix(1:3)
    ),
    "locs must be a list for multivariate models"
  )
})

test_that("Invalid Y input", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, as.matrix(1:3), list(as.matrix(1:3)),
      list(as.matrix(1:3))
    ),
    "Y must be a list for multivariate models"
  )
})

test_that("Invalid X input", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(1:3)), as.matrix(1:3),
      list(as.matrix(1:3))
    ),
    "X must be a list for multivariate models"
  )
})

test_that("locs/Y length mismatch", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(1:3, 2:4), list(as.matrix(1:3)),
      list(as.matrix(1:3))
    ),
    "locs and Y must have the same length"
  )
})

test_that("locs/X length mismatch", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(1:3), list(
        as.matrix(1:3),
        as.matrix(2:4)
      ),
      list(as.matrix(1:3))
    ),
    "locs and X must have the same length"
  )
})

test_that("Invalid lod input", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(1:3), list(as.matrix(1:3)),
      list(as.matrix(1:3)), lod = 1:2
    ),
    "lod must be a list for multivariate models"
  )
})

test_that("locs/lod length mismatch", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(1:3), list(as.matrix(1:3)),
      list(as.matrix(1:3)), lod = list(1, 2)
    ),
    "locs and lod must have the same length"
  )
})

test_that("locs not a matrix", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(1:3), list(as.matrix(1:3)),
      list(1:3)
    ),
    "Each locs must be a matrix"
  )
})

test_that("X not a matrix", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(1:3), list(1:3),
      list(as.matrix(1:3))
    ),
    "Each X must be a matrix"
  )
})

test_that("Y not a numeric matrix/vector", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list("foo"), list(as.matrix(1:3)),
      list(as.matrix(1:3))
    ),
    "Each Y must be a numeric vector or matrix"
  )
})

test_that("locs with differing numbers of columns", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(1:4, 1:4),
      list(as.matrix(1:4), as.matrix(1:4)),
      list(as.matrix(1:4), matrix(1:4, nrow = 2))
    ),
    "All locs must have the same number of columns"
  )
})

test_that("Y has multiple columns", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(matrix(1:4, nrow = 2)),
      list(as.matrix(1:4)), list(as.matrix(1:4))
    ),
    "Each Y must have only 1 column"
  )
})

test_that("nrow(Y) != nrow(locs)", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(1:4)),
      list(as.matrix(1:4)), list(as.matrix(1:3))
    ),
    "Each Y must have the same number of rows as locs"
  )
})

test_that("nrow(Y) != nrow(X)", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(1:4)),
      list(as.matrix(1:3)), list(as.matrix(1:4))
    ),
    "Each Y must have the same number of rows as X"
  )
})

test_that("NA's in X", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(1:4)), list(as.matrix(c(1:3, NA))),
      list(as.matrix(1:4))
    ),
    "X must not contain NA's"
  )
})

test_that("NA's in locs", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(1:4)), list(as.matrix(1:4)),
      list(as.matrix(c(1:3, NA)))
    ),
    "locs must not contain NA's"
  )
})

test_that("NA's in Y", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(c(1:3, NA))), list(as.matrix(1:4)),
      list(as.matrix(1:4))
    ),
    "Y contains NA's and impute.y is FALSE. Set impute.y=TRUE to impute missing Y's."
  )
})

test_that("lod not numeric", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(c(1:4))), list(as.matrix(1:4)),
      list(as.matrix(1:4)), lod = list("foo")
    ),
    "Each lod must be numeric"
  )
})

test_that("length(lod) != nrow(X)", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(c(1:4))), list(as.matrix(1:4)),
      list(as.matrix(1:4)), lod = list(1:2)
    ),
    "Length of each lod must equal the number of observations"
  )
})

test_that("length(Y.names) != length(locs)", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(c(1:4))), list(as.matrix(1:4)),
      list(as.matrix(1:4)), Y.names = c("test1", "test2"),
    ),
    "Length of Y.names must match the number of response variables"
  )
})

test_that("X.names is not a list", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(c(1:4))), list(as.matrix(1:4)),
      list(as.matrix(1:4)), X.names = c("test1", "test2"),
    ),
    "X.names must be a list for multivariate models"
  )
})

test_that("length(X.names) != length(locs)", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(c(1:4))), list(as.matrix(1:4)),
      list(as.matrix(1:4)), X.names = list(1, 2),
    ),
    "Length of X.names must match the number of response variables"
  )
})

test_that("length(X.names[[i]]) != length(locs[[i]])", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(c(1:4))), list(as.matrix(1:4)),
      list(as.matrix(1:4)), X.names = list(c("test1", "test2")),
    ),
    "Length of each X.names must match the number of predictors"
  )
})

test_that("Simulated dataset multivariate spatial", {
  load("sim_multivariate_big.RData")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 25)
  pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list, X.st, locs.list,
    Y.names = c("test1", "test2", "test3"), X.names = list(c("x1_1", "x1_2",
        "x1_3", "x1_4", "x1_5", "x1_6", "x1_7", "x1_8", "x1_9", "x1_10"),
      c("x2_1", "x2_2", "x2_3", "x2_4", "x2_5", "x2_6", "x2_7", "x2_8",
        "x2_9", "x2_10"), c("x3_1", "x3_2", "x3_3", "x3_4", "x3_5", "x3_6",
        "x3_7", "x3_8", "x3_9", "x3_10")),
    scaling = c(1, 1), apanasovich = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )

  expect_true(validObject(pgp.mmodel1))

  expect_equal(get_neighbors(pgp.mmodel1), pgp.mmodel1@n_neighbors)
  expect_equal(get_scaling(pgp.mmodel1), pgp.mmodel1@scaling)
  expect_equal(get_converged(pgp.mmodel1), pgp.mmodel1@converged)
  expect_equal(get_pen_loglik(pgp.mmodel1), pgp.mmodel1@error)

  beta.out <- get_beta(pgp.mmodel1)
  params.out <- pgp.mmodel1@covparams
  theta.out <- get_theta(pgp.mmodel1)

  expect_length(beta.out, 4)
  expect_length(beta.out[[1]], 10)
  expect_length(beta.out[[2]], 10)
  expect_length(beta.out[[3]], 10)
  expect_length(beta.out[[4]], 3)
  expect_length(params.out, 15)
  expect_length(theta.out[[1]], 3)
  expect_length(theta.out[[2]], 3)
  expect_length(theta.out[[3]], 3)
  expect_length(theta.out[[4]], 3)
  expect_equal(dim(theta.out[[5]]), c(3, 3))

  expect_named(theta.out, c("sigma", "scale", "smoothness", "nuggets",
      "correlation"))
  expect_named(theta.out[[1]], c("test1", "test2", "test3"))
  expect_named(theta.out[[2]], c("test1", "test2", "test3"))
  expect_named(theta.out[[3]], c("test1", "test2", "test3"))
  expect_named(theta.out[[4]], c("test1", "test2", "test3"))
  expect_equal(colnames(theta.out[[5]]), c("test1", "test2", "test3"))
  expect_equal(rownames(theta.out[[5]]), c("test1", "test2", "test3"))

  expect_named(beta.out, c("test1", "test2", "test3", "(Intercept)"))
  expect_named(beta.out[[1]], c("x1_1", "x1_2", "x1_3", "x1_4", "x1_5",
      "x1_6", "x1_7", "x1_8", "x1_9", "x1_10"))
  expect_named(beta.out[[2]], c("x2_1", "x2_2", "x2_3", "x2_4", "x2_5",
      "x2_6", "x2_7", "x2_8", "x2_9", "x2_10"))
  expect_named(beta.out[[3]], c("x3_1", "x3_2", "x3_3", "x3_4", "x3_5",
      "x3_6", "x3_7", "x3_8", "x3_9", "x3_10"))
  expect_named(beta.out[[4]], c("(Intercept)", "test2", "test3"))

  expect_identical(as.numeric(theta.out[[1]]), params.out[1:3])
  expect_identical(as.numeric(theta.out[[2]]), params.out[4:6])
  expect_identical(as.numeric(theta.out[[3]]), params.out[7:9])
  expect_identical(as.numeric(theta.out[[4]]), params.out[10:12])

  rho.mat <- matrix(0, nrow = 3, ncol = 3)
  rho.mat[upper.tri(rho.mat, diag = FALSE)] <- params.out[13:15]
  rho.mat <- rho.mat + t(rho.mat)
  diag(rho.mat) <- 1
  expect_equal(matrix(theta.out[[5]], nrow = nrow(theta.out[[5]])), rho.mat)

  expect_identical(as.vector(beta.out[[1]]), as.vector(pgp.mmodel1@beta[2:11]))
  expect_identical(as.vector(beta.out[[2]]), as.vector(pgp.mmodel1@beta[13:22]))
  expect_identical(as.vector(beta.out[[3]]), as.vector(pgp.mmodel1@beta[24:33]))
  expect_identical(as.numeric(beta.out[[4]][1]),
    as.numeric(pgp.mmodel1@beta[1]))
  expect_identical(as.numeric(beta.out[[4]][2]),
    as.numeric(pgp.mmodel1@beta[12]))
  expect_identical(as.numeric(beta.out[[4]][3]),
    as.numeric(pgp.mmodel1@beta[23]))

  expect_equal(as.numeric(beta.out[[1]]), c(0.98, 0.97, 0.92, 0.96, rep(0, 3),
      0.02, rep(0, 2)), tolerance = 0.04)
  expect_equal(as.numeric(beta.out[[2]]), c(0.86, 1.0, 0.8, 0.98, rep(0, 3),
      0.04, rep(0, 2)), tolerance = 0.04)
  expect_equal(as.numeric(beta.out[[3]]), c(1.0, 0.99, 0.96, 0.94, rep(0, 6)),
    tolerance = 0.04)
  expect_equal(as.numeric(beta.out[[4]]), c(-0.01, 0, 0), tolerance = 0.04)

  expect_equal(params.out[1], 2.3, tolerance = 5.1)
  expect_equal(params.out[2], 4.4, tolerance = 4.6)
  expect_equal(params.out[3], 1.9, tolerance = 4.3)
  expect_equal(params.out[4], 0.35, tolerance = 1.6)
  expect_equal(params.out[5], 0.46, tolerance = 2.6)
  expect_equal(params.out[6], 0.22, tolerance = 0.4)
  expect_equal(params.out[7], 0.66, tolerance = 0.8)
  expect_equal(params.out[8] - 0.63, 0, tolerance = 0.5)
  expect_equal(params.out[9], 1, tolerance = 0.7)
  expect_equal(params.out[10], 1, tolerance = 0.8)
  expect_equal(params.out[11], 1.5, tolerance = 1.4)
  expect_equal(params.out[12], 0.47, tolerance = 0.1)
  expect_equal(params.out[13], 0.39, tolerance = 0.4)
  expect_equal(params.out[14], 0.7, tolerance = 0.4)
  expect_equal(params.out[15], 0.54, tolerance = 0.5)

  # Missing data
  set.seed(1234)
  y.list.na <- y.list
  y.list.lod <- y.list
  lod.cut <- list()
  for (i in seq_along(y.list)) {
    y.list.na[[i]][sample(seq_along(y.list[[i]]),
        floor(0.1 * length(y.list[[i]])))] <- NA
    y.list.lod[[i]] <- y.list.lod[[i]] + 10
    lod.cut[[i]] <- quantile(y.list.lod[[i]], 0.1)
    y.list.lod[[i]][y.list.lod[[i]] <= lod.cut[[i]]] <- NA
  }

  pgp.mmodel2 <- new("MultivariateVecchiaModel", n_neighbors = 25)
  pgp.mmodel2 <- prestogp_fit(pgp.mmodel2, y.list.na, X.st, locs.list,
    scaling = c(1, 1), apanasovich = TRUE, quiet = TRUE,
    impute.y = TRUE, optim.control = list(trace = 0, maxit = 5000,
      reltol = 1e-3))
  beta.out2 <- get_beta(pgp.mmodel2)
  params.out2 <- pgp.mmodel2@covparams
  theta.out2 <- get_theta(pgp.mmodel2)

  expect_length(beta.out2, 4)
  expect_length(beta.out2[[1]], 10)
  expect_length(beta.out2[[2]], 10)
  expect_length(beta.out2[[3]], 10)
  expect_length(beta.out2[[4]], 3)
  expect_length(params.out2, 15)
  expect_length(theta.out2[[1]], 3)
  expect_length(theta.out2[[2]], 3)
  expect_length(theta.out2[[3]], 3)
  expect_length(theta.out2[[4]], 3)
  expect_equal(dim(theta.out2[[5]]), c(3, 3))

  expect_named(theta.out2[[1]], c("Y1", "Y2", "Y3"))
  expect_named(theta.out2[[2]], c("Y1", "Y2", "Y3"))
  expect_named(theta.out2[[3]], c("Y1", "Y2", "Y3"))
  expect_named(theta.out2[[4]], c("Y1", "Y2", "Y3"))
  expect_equal(colnames(theta.out2[[5]]), c("Y1", "Y2", "Y3"))
  expect_equal(rownames(theta.out2[[5]]), c("Y1", "Y2", "Y3"))

  expect_named(beta.out2, c("Y1", "Y2", "Y3", "(Intercept)"))
  expect_named(beta.out2[[4]], c("(Intercept)", "Y2", "Y3"))

  # Results should be the same after imputation
  expect_equal(as.numeric(beta.out[[1]]), as.numeric(beta.out2[[1]]),
    tolerance = 0.05)
  expect_equal(as.numeric(beta.out[[2]]), as.numeric(beta.out2[[2]]),
    tolerance = 0.05)
  expect_equal(as.numeric(beta.out[[3]]), as.numeric(beta.out2[[3]]),
    tolerance = 0.05)
  expect_equal(as.numeric(beta.out[[4]]), as.numeric(beta.out2[[4]]),
    tolerance = 0.05)
  expect_equal(params.out[1], params.out2[1], tolerance = 6.8)
  expect_equal(params.out[2], params.out2[2], tolerance = 10.0)
  expect_equal(params.out[3], params.out2[3], tolerance = 5.2)
  expect_equal(params.out[4], params.out2[4], tolerance = 2.5)
  expect_equal(params.out[5], params.out2[5], tolerance = 1.9)
  expect_equal(params.out[6], params.out2[6], tolerance = 0.6)
  expect_equal(params.out[7], params.out2[7], tolerance = 1.4)
  expect_equal(params.out[8], params.out2[8], tolerance = 1.1)
  expect_equal(params.out[9], params.out2[9], tolerance = 1.0)
  expect_equal(params.out[10], params.out2[10], tolerance = 0.9)
  expect_equal(params.out[11], params.out2[11], tolerance = 1.4)
  expect_equal(params.out[12], params.out2[12], tolerance = 0.2)
  expect_equal(params.out[13], params.out2[13], tolerance = 0.6)
  expect_equal(params.out[14] - params.out2[14], 0, tolerance = 0.4)
  expect_equal(params.out[15], params.out2[15], tolerance = 0.5)

  # Missing data with lod
  pgp.mmodel3 <- new("MultivariateVecchiaModel", n_neighbors = 25)
  pgp.mmodel3 <- prestogp_fit(pgp.mmodel3, y.list.lod, X.st, locs.list,
    scaling = c(1, 1), apanasovich = TRUE, verbose = TRUE,
    impute.y = TRUE, lod = lod.cut, optim.control = list(trace = 0,
      maxit = 5000, reltol = 1e-3))
  beta.out3 <- get_beta(pgp.mmodel3)
  params.out3 <- pgp.mmodel3@covparams

  # Results should be the same after imputation
  expect_equal(as.numeric(beta.out[[1]]), as.numeric(beta.out3[[1]]),
    tolerance = 0.05)
  expect_equal(as.numeric(beta.out[[2]]), as.numeric(beta.out3[[2]]),
    tolerance = 0.05)
  expect_equal(as.numeric(beta.out[[3]]), as.numeric(beta.out3[[3]]),
    tolerance = 0.05)
  expect_equal(as.numeric(beta.out[[4]]), as.numeric(beta.out3[[4]]),
    tolerance = 0.05)
  expect_equal(params.out[1], params.out3[1], tolerance = 7.5)
  expect_equal(params.out[2], params.out3[2], tolerance = 7.2)
  expect_equal(params.out[3], params.out3[3], tolerance = 7.1)
  expect_equal(params.out[4], params.out3[4], tolerance = 3.7)
  expect_equal(params.out[5], params.out3[5], tolerance = 1.7)
  expect_equal(params.out[6], params.out3[6], tolerance = 0.4)
  expect_equal(params.out[7], params.out3[7], tolerance = 1.0)
  expect_equal(params.out[8], params.out3[8], tolerance = 0.8)
  expect_equal(params.out[9], params.out3[9], tolerance = 1.0)
  expect_equal(params.out[10], params.out3[10], tolerance = 0.8)
  expect_equal(params.out[11], params.out3[11], tolerance = 1.6)
  expect_equal(params.out[12], params.out3[12], tolerance = 0.3)
  expect_equal(params.out[13], params.out3[13], tolerance = 0.5)
  expect_equal(params.out[14], params.out3[14], tolerance = 0.4)
  expect_equal(params.out[15], params.out3[15], tolerance = 0.6)
})

test_that("Simulated dataset multivariate spatiotemporal", {
  source("sim_multivariate_big_st.R")
  names(y.list) <- c("test1", "test2", "test3")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 25)
  pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list, X.st, locs.list,
    scaling = c(1, 1, 2), quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  beta.out <- as.vector(pgp.mmodel1@beta)
  params.out <- pgp.mmodel1@covparams
  theta.out <- get_theta(pgp.mmodel1)

  expect_length(beta.out, 33)
  expect_length(params.out, 18)
  expect_length(theta.out[[1]], 3)
  expect_length(theta.out[[2]], 6)
  expect_length(theta.out[[3]], 3)
  expect_length(theta.out[[4]], 3)
  expect_equal(dim(theta.out[[5]]), c(3, 3))

  expect_named(theta.out, c("sigma", "scale", "smoothness", "nuggets",
      "correlation"))
  expect_named(theta.out[[1]], c("test1", "test2", "test3"))
  expect_named(theta.out[[2]], c("test1_1", "test1_2", "test2_1", "test2_2",
      "test3_1", "test3_2"))
  expect_named(theta.out[[3]], c("test1", "test2", "test3"))
  expect_named(theta.out[[4]], c("test1", "test2", "test3"))
  expect_equal(colnames(theta.out[[5]]), c("test1", "test2", "test3"))
  expect_equal(rownames(theta.out[[5]]), c("test1", "test2", "test3"))

  expect_identical(as.numeric(theta.out[[1]]), params.out[1:3])
  expect_identical(as.numeric(theta.out[[2]]), params.out[4:9])
  expect_identical(as.numeric(theta.out[[3]]), params.out[10:12])
  expect_identical(as.numeric(theta.out[[4]]), params.out[13:15])

  rho.mat <- matrix(0, nrow = 3, ncol = 3)
  rho.mat[upper.tri(rho.mat, diag = FALSE)] <- params.out[16:18]
  rho.mat <- rho.mat + t(rho.mat)
  diag(rho.mat) <- 1
  expect_equal(matrix(theta.out[[5]], nrow = nrow(theta.out[[5]])), rho.mat)

  expect_equal(beta.out[1:24], c(
    0, 0.82, 0.95, 0.85, 1.02, rep(0, 7), 0.92,
    0.95, 1.01, 0.76, rep(0, 7), 0.81), tolerance = 0.2)
  expect_equal(beta.out[27] - 0.84, 0, tolerance = 0.2)
  expect_equal(beta.out[28:33], rep(0, 6), tolerance = 0.04)
  expect_gt(beta.out[25], 0.78)
  expect_lt(beta.out[25], 1.27)
  expect_gt(beta.out[26], 0.58)
  expect_lt(beta.out[26], 1.02)

  expect_gt(params.out[1], 1.4)
  expect_lt(params.out[1], 15.7)
  expect_gt(params.out[2], 3.2)
  expect_lt(params.out[2], 25.5)
  expect_gt(params.out[3], 3)
  expect_lt(params.out[3], 49.8)
  expect_gt(params.out[4], 0.04)
  expect_lt(params.out[4], 1.06)
  expect_gt(params.out[5], 5.3)
  expect_lt(params.out[5], 65.6)
  expect_gt(params.out[6], 0.04)
  expect_lt(params.out[6], 1.51)
  expect_gt(params.out[7], 5.9)
  expect_lt(params.out[7], 73.5)
  expect_gt(params.out[8], 0.02)
  expect_lt(params.out[8], 3.74)
  expect_gt(params.out[9], 0.04)
  expect_lt(params.out[9], 91.1)
  expect_gt(params.out[10], 0.4)
  expect_lt(params.out[10], 1.48)
  expect_gt(params.out[11], 0.3)
  expect_lt(params.out[11], 2.19)
  expect_gt(params.out[12], 0.48)
  expect_lt(params.out[12], 2.19)
  expect_gt(params.out[13], 0.01)
  expect_lt(params.out[13], 1.65)
  expect_gt(params.out[14], 0.04)
  expect_lt(params.out[14], 3.01)
  expect_gt(params.out[15], 0.4)
  expect_lt(params.out[15], 3.95)
  expect_gt(params.out[16], -0.5)
  expect_gt(params.out[17], -0.5)
  expect_gt(params.out[18], -0.5)
})

test_that("Invalid locs input for prediction", {
  source("sim_multivariate_small_pred.R")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 5)
  pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
    locs.list.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  expect_error(
    prestogp_predict(pgp.mmodel1, X.st.otst, as.matrix(1:3)),
    "locs must be a list for multivariate models"
  )
})

test_that("Invalid X input for prediction", {
  source("sim_multivariate_small_pred.R")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 5)
  pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
    locs.list.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  expect_error(
    prestogp_predict(pgp.mmodel1, as.matrix(1:3), locs.list.otst),
    "X must be a list for multivariate models"
  )
})

test_that("locs/X length mismatch for prediction", {
  source("sim_multivariate_small_pred.R")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 5)
  pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
    locs.list.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  expect_error(
    prestogp_predict(pgp.mmodel1, X.st.otst, list(as.matrix(1:3))),
    "locs and X must have the same length"
  )
})

test_that("locs/locs_train length mismatch", {
  source("sim_multivariate_small_pred.R")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 5)
  pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
    locs.list.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  expect_error(
    prestogp_predict(
      pgp.mmodel1, list(1:3, 1:3, 1:3),
      list(1:3, 1:3, 1:3)
    ),
    "Training and test set locs must have the same length"
  )
})

test_that("locs not a matrix for prediction", {
  source("sim_multivariate_small_pred.R")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 5)
  pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
    locs.list.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  expect_error(
    prestogp_predict(
      pgp.mmodel1, list(
        as.matrix(1:3),
        as.matrix(1:3)
      ),
      list(1:3, as.matrix(1:3))
    ),
    "Each locs must be a matrix"
  )
})

test_that("ncol(locs) != ncol(locs_train)", {
  source("sim_multivariate_small_pred.R")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 5)
  pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
    locs.list.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  expect_error(
    prestogp_predict(
      pgp.mmodel1, list(
        as.matrix(1:3),
        as.matrix(1:3)
      ),
      list(as.matrix(1:3), as.matrix(1:3))
    ),
    "All locs must have the same number of columns as locs_train"
  )
})

test_that("X not a matrix for prediction", {
  source("sim_multivariate_small_pred.R")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 5)
  pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
    locs.list.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  X.st.otst[[1]] <- as.vector(X.st.otst[[1]])
  expect_error(
    prestogp_predict(pgp.mmodel1, X.st.otst, locs.list.otst),
    "Each X must be a matrix"
  )
})

test_that("nrow(X) != nrow(locs) for prediction", {
  source("sim_multivariate_small_pred.R")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 5)
  pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
    locs.list.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  locs.list.otst[[1]] <- rbind(locs.list.otst[[1]], rep(0, 2))
  expect_error(
    prestogp_predict(pgp.mmodel1, X.st.otst, locs.list.otst),
    "Each X must have the same number of rows as locs"
  )
})

test_that("ncol(X) != ncol(X_train) for prediction", {
  source("sim_multivariate_small_pred.R")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 5)
  pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
    locs.list.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  X.st.otst[[1]] <- cbind(X.st.otst[[1]], rep(0, nrow(X.st.otst[[1]])))
  expect_error(
    prestogp_predict(pgp.mmodel1, X.st.otst, locs.list.otst),
    "X and X_train must have the same number of predictors"
  )
})

test_that("m too small for prediction", {
  source("sim_multivariate_small_pred.R")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 5)
  pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
    locs.list.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  expect_error(
    prestogp_predict(pgp.mmodel1, X.st.otst, locs.list.otst, m = 2),
    "m must be at least 3"
  )
})

test_that("m too large for prediction", {
  source("sim_multivariate_small_pred.R")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 5)
  pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
    locs.list.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  expect_warning(
    prestogp_predict(pgp.mmodel1, X.st.otst, locs.list.otst,
      m = 51
    ),
    "Conditioning set size m chosen to be >=n. Changing to m=n-1"
  )
})

test_that("Simulated dataset multivariate spatial prediction", {
  source("sim_multivariate_big_pred.R")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 25)
  pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
    locs.list.otr,
    scaling = c(1, 1),
    apanasovich = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )

  pgp.mmodel1.pred <- prestogp_predict(pgp.mmodel1, X.st.otst, locs.list.otst)

  mse <- mean((pgp.mmodel1.pred$means - unlist(y.list.otst))^2)
  me <- mean(pgp.mmodel1.pred$means - unlist(y.list.otst))

  expect_equal(mse, 1.99, tolerance = 0.1)
  expect_equal(me + 0.04, 0, tolerance = 0.03)
})
