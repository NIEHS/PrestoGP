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

test_that("Invalid lod.upper input", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(1:3), list(as.matrix(1:3)),
      list(as.matrix(1:3)), lod.upper = 1:2
    ),
    "lod.upper must be a list for multivariate models"
  )
})

test_that("Invalid lod.lower input", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(1:3), list(as.matrix(1:3)),
      list(as.matrix(1:3)), lod.lower = 1:2
    ),
    "lod.lower must be a list for multivariate models"
  )
})

test_that("locs/lod.upper length mismatch", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(1:3), list(as.matrix(1:3)),
      list(as.matrix(1:3)), lod.upper = list(1, 2)
    ),
    "locs and lod.upper must have the same length"
  )
})

test_that("locs/lod.lower length mismatch", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(1:3), list(as.matrix(1:3)),
      list(as.matrix(1:3)), lod.lower = list(1, 2)
    ),
    "locs and lod.lower must have the same length"
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

test_that("lod.upper not numeric", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(c(1:4))), list(as.matrix(1:4)),
      list(as.matrix(1:4)), lod.upper = list("foo")
    ),
    "Each lod.upper must be numeric"
  )
})

test_that("lod.lower not numeric", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(c(1:4))), list(as.matrix(1:4)),
      list(as.matrix(1:4)), lod.lower = list("foo")
    ),
    "Each lod.lower must be numeric"
  )
})

test_that("length(lod.upper) != nrow(X)", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(c(1:4))), list(as.matrix(1:4)),
      list(as.matrix(1:4)), lod.upper = list(1:2)
    ),
    "Length of each lod.upper must equal the number of observations"
  )
})

test_that("length(lod.lower) != nrow(X)", {
  model <- new("MultivariateVecchiaModel")
  expect_error(
    prestogp_fit(
      model, list(as.matrix(c(1:4))), list(as.matrix(1:4)),
      list(as.matrix(1:4)), lod.lower = list(1:2)
    ),
    "Length of each lod.lower must equal the number of observations"
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
    scaling = c(1, 1), common_scale = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )

  expect_silent(plot_beta(pgp.mmodel1))
  dev.off()

  expect_true(validObject(pgp.mmodel1))

  expect_equal(get_Y(pgp.mmodel1), y.list)
  expect_equal(get_linear_model(pgp.mmodel1), pgp.mmodel1@linear_model)
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
  expect_equal(as.numeric(beta.out[[4]][1] - mean(y.list[[1]])),
    as.numeric(pgp.mmodel1@beta[1]))
  expect_identical(as.numeric(beta.out[[4]][2] + mean(y.list[[1]]) -
        mean(y.list[[2]])), as.numeric(pgp.mmodel1@beta[12]))
  expect_identical(as.numeric(beta.out[[4]][3] + mean(y.list[[1]]) -
        mean(y.list[[3]])), as.numeric(pgp.mmodel1@beta[23]))

  expect_equal(as.numeric(beta.out[[1]]), c(1.05, 0.98, 0.97, 1.0, 0.02, 0,
      0, 0.01, rep(0, 2)), tolerance = 0.02)
  expect_equal(as.numeric(beta.out[[2]]), c(0.95, 1.06, 0.83, 1.03, rep(0, 3),
      0.06, rep(0, 2)), tolerance = 0.02)
  expect_equal(as.numeric(beta.out[[3]]), c(1.015, 1.04, 1.0, 0.95, 0.01,
      rep(0, 4), 0.01), tolerance = 0.02)
  expect_equal(as.numeric(beta.out[[4]]), c(-0.01 + mean(y.list[[1]]),
      mean(y.list[[2]]) - mean(y.list[[1]]),
      mean(y.list[[3]]) - mean(y.list[[1]])), tolerance = 0.04)

  expect_gt(params.out[1], 1.0)
  expect_lt(params.out[1], 3.4)
  expect_gt(params.out[2], 3.8)
  expect_lt(params.out[2], 10.9)
  expect_gt(params.out[3], 0.9)
  expect_lt(params.out[3], 1.6)
  expect_gt(params.out[4], 0.06)
  expect_lt(params.out[4], 0.22)
  expect_gt(params.out[5], 0.24)
  expect_lt(params.out[5], 0.53)
  expect_gt(params.out[6], 0.1)
  expect_lt(params.out[6], 0.24)
  expect_gt(params.out[7], 0.64)
  expect_lt(params.out[7], 1.88)
  expect_gt(params.out[8], 0.55)
  expect_lt(params.out[8], 0.99)
  expect_gt(params.out[9], 0.68)
  expect_lt(params.out[9], 1.38)
  expect_gt(params.out[10], 1.02)
  expect_lt(params.out[10], 1.21)
  expect_gt(params.out[11], 1.42)
  expect_lt(params.out[11], 1.84)
  expect_gt(params.out[12], 0.39)
  expect_lt(params.out[12], 0.52)
  expect_gt(params.out[13], 0.12)
  expect_lt(params.out[13], 0.6)
  expect_gt(params.out[14], 0.63)
  expect_lt(params.out[14], 0.83)
  expect_gt(params.out[15], 0.44)
  expect_lt(params.out[15], 0.76)

  # Missing data
  set.seed(1234)
  y.list.na <- y.list
  y.list.lod <- y.list
  lod.cut <- list()
  lodupper <- list()
  lodlower <- list()
  for (i in seq_along(y.list)) {
    y.list.na[[i]][sample(seq_along(y.list[[i]]),
        floor(0.1 * length(y.list[[i]])))] <- NA
    y.list.lod[[i]] <- y.list.lod[[i]] + 10
    lod.cut[[i]] <- quantile(y.list.lod[[i]], 0.1)
    y.list.lod[[i]][y.list.lod[[i]] <= lod.cut[[i]]] <- NA
    lodupper[[i]] <- rep(lod.cut[[i]], length(y.list.lod[[i]]))
    lodlower[[i]] <- rep(0, length(y.list.lod[[i]]))
  }

  pgp.mmodel2 <- new("MultivariateVecchiaModel", n_neighbors = 25)
  pgp.mmodel2 <- prestogp_fit(pgp.mmodel2, y.list.na, X.st, locs.list,
    scaling = c(1, 1), common_scale = TRUE, quiet = TRUE,
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
  expect_named(beta.out2[[1]], c("Y1_1", "Y1_2", "Y1_3", "Y1_4", "Y1_5",
      "Y1_6", "Y1_7", "Y1_8", "Y1_9", "Y1_10"))
  expect_named(beta.out2[[2]], c("Y2_1", "Y2_2", "Y2_3", "Y2_4", "Y2_5",
      "Y2_6", "Y2_7", "Y2_8", "Y2_9", "Y2_10"))
  expect_named(beta.out2[[3]], c("Y3_1", "Y3_2", "Y3_3", "Y3_4", "Y3_5",
      "Y3_6", "Y3_7", "Y3_8", "Y3_9", "Y3_10"))
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
  expect_equal(params.out[1], params.out2[1], tolerance = 3.5)
  expect_equal(params.out[2], params.out2[2], tolerance = 7.8)
  expect_equal(params.out[3], params.out2[3], tolerance = 1.5)
  expect_equal(params.out[4], params.out2[4], tolerance = 0.7)
  expect_equal(params.out[5], params.out2[5], tolerance = 0.6)
  expect_equal(params.out[6], params.out2[6], tolerance = 0.3)
  expect_equal(params.out[7], params.out2[7], tolerance = 1.3)
  expect_equal(params.out[8], params.out2[8], tolerance = 0.6)
  expect_equal(params.out[9], params.out2[9], tolerance = 0.7)
  expect_equal(params.out[10], params.out2[10], tolerance = 0.3)
  expect_equal(params.out[11], params.out2[11], tolerance = 0.4)
  expect_equal(params.out[12] - params.out2[12], 0, tolerance = 0.15)
  expect_equal(params.out[13] - params.out2[13], 0, tolerance = 0.5)
  expect_equal(params.out[14] - params.out2[14], 0, tolerance = 0.25)
  expect_equal(params.out[15], params.out2[15], tolerance = 0.25)

  # Missing data with lod
  doParallel::registerDoParallel(cores = 2)
  pgp.mmodel3 <- new("MultivariateVecchiaModel", n_neighbors = 25)
  pgp.mmodel3 <- prestogp_fit(pgp.mmodel3, y.list.lod, X.st, locs.list,
    scaling = c(1, 1), common_scale = TRUE, verbose = TRUE, parallel = TRUE,
    impute.y = TRUE, lod.upper = lodupper, lod.lower = lodlower,
    optim.control = list(trace = 0,
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
  expect_equal(as.numeric(beta.out[[4]][1] - pgp.mmodel1@Y_bar[1]),
    as.numeric(beta.out3[[4]][1] - pgp.mmodel3@Y_bar[1]), tolerance = 0.05)
  expect_equal(as.numeric(beta.out[[4]][2] - pgp.mmodel1@Y_bar[2] +
        pgp.mmodel1@Y_bar[1]), as.numeric(beta.out3[[4]][2] -
        pgp.mmodel3@Y_bar[2] + pgp.mmodel3@Y_bar[1]), tolerance = 0.05)
  expect_equal(as.numeric(beta.out[[4]][3] - pgp.mmodel1@Y_bar[3] +
        pgp.mmodel1@Y_bar[1]), as.numeric(beta.out3[[4]][3] -
        pgp.mmodel3@Y_bar[3] + pgp.mmodel3@Y_bar[1]), tolerance = 0.05)
  expect_equal(params.out[1], params.out3[1], tolerance = 2.2)
  expect_equal(params.out[2], params.out3[2], tolerance = 4.3)
  expect_equal(params.out[3], params.out3[3], tolerance = 0.6)
  expect_equal(params.out[4] - params.out3[4], 0, tolerance = 0.28)
  expect_equal(params.out[5], params.out3[5], tolerance = 0.42)
  expect_equal(params.out[6] - params.out3[6], 0, tolerance = 0.23)
  expect_equal(params.out[7], params.out3[7], tolerance = 1.0)
  expect_equal(params.out[8] - params.out3[8], 0, tolerance = 0.5)
  expect_equal(params.out[9], params.out3[9], tolerance = 0.8)
  expect_equal(params.out[10] - params.out3[10], 0, tolerance = 0.25)
  expect_equal(params.out[11], params.out3[11], tolerance = 0.52)
  expect_equal(params.out[12] - params.out3[12], 0, tolerance = 0.12)
  expect_equal(params.out[13] - params.out3[13], 0, tolerance = 0.4)
  expect_equal(params.out[14] - params.out3[14], 0, tolerance = 0.15)
  expect_equal(params.out[15], params.out3[15], tolerance = 0.3)

  # SCAD fit
  pgp.mmodel4 <- new("MultivariateVecchiaModel", n_neighbors = 25)
  pgp.mmodel4 <- prestogp_fit(pgp.mmodel4, y.list, X.st, locs.list,
    scaling = c(1, 1), common_scale = TRUE, quiet = TRUE, penalty = "SCAD",
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  beta.out4 <- get_beta(pgp.mmodel4)
  params.out4 <- pgp.mmodel4@covparams

  expect_equal(as.numeric(beta.out4[[1]]), c(1.08, 1.01, 1.01, 1.01, rep(0, 6)),
    tolerance = 0.01)
  expect_equal(as.numeric(beta.out4[[2]]), c(0.98, 1.07, 0.88, 1.06, rep(0, 6)),
    tolerance = 0.01)
  expect_equal(as.numeric(beta.out4[[3]]), c(1.03, 1.05, 1.01, 0.97,
      rep(0, 6)), tolerance = 0.01)
  expect_equal(as.numeric(beta.out4[[4]]), c(-0.01 + mean(y.list[[1]]),
      mean(y.list[[2]]) - mean(y.list[[1]]),
      mean(y.list[[3]]) - mean(y.list[[1]])), tolerance = 0.005)

  expect_gt(params.out4[1], 1.0)
  expect_lt(params.out4[1], 3.5)
  expect_gt(params.out4[2], 2.9)
  expect_lt(params.out4[2], 14.1)
  expect_gt(params.out4[3], 0.9)
  expect_lt(params.out4[3], 2.2)
  expect_gt(params.out4[4], 0.05)
  expect_lt(params.out4[4], 0.28)
  expect_gt(params.out4[5], 0.23)
  expect_lt(params.out4[5], 1.06)
  expect_gt(params.out4[6], 0.09)
  expect_lt(params.out4[6], 0.23)
  expect_gt(params.out4[7], 0.67)
  expect_lt(params.out4[7], 1.98)
  expect_gt(params.out4[8], 0.47)
  expect_lt(params.out4[8], 0.97)
  expect_gt(params.out4[9], 0.77)
  expect_lt(params.out4[9], 1.64)
  expect_gt(params.out4[10], 1.02)
  expect_lt(params.out4[10], 1.21)
  expect_gt(params.out4[11], 1.48)
  expect_lt(params.out4[11], 1.81)
  expect_gt(params.out4[12], 0.4)
  expect_lt(params.out4[12], 0.5)
  expect_gt(params.out4[13], 0.13)
  expect_lt(params.out4[13], 0.64)
  expect_gt(params.out4[14], 0.63)
  expect_lt(params.out4[14], 0.83)
  expect_gt(params.out4[15], 0.41)
  expect_lt(params.out4[15], 0.8)

  # Missing data
  pgp.mmodel5 <- new("MultivariateVecchiaModel", n_neighbors = 25)
  pgp.mmodel5 <- prestogp_fit(pgp.mmodel5, y.list.na, X.st, locs.list,
    scaling = c(1, 1), common_scale = TRUE, quiet = TRUE, penalty = "SCAD",
    impute.y = TRUE, optim.control = list(trace = 0, maxit = 5000,
      reltol = 1e-3))
  beta.out5 <- get_beta(pgp.mmodel5)
  params.out5 <- pgp.mmodel5@covparams

  # Results should be the same after imputation
  expect_equal(as.numeric(beta.out4[[1]]), as.numeric(beta.out5[[1]]),
    tolerance = 0.04)
  expect_equal(as.numeric(beta.out4[[2]]), as.numeric(beta.out5[[2]]),
    tolerance = 0.04)
  expect_equal(as.numeric(beta.out4[[3]]), as.numeric(beta.out5[[3]]),
    tolerance = 0.04)
  expect_equal(as.numeric(beta.out4[[4]]), as.numeric(beta.out5[[4]]),
    tolerance = 0.04)
  expect_equal(params.out4[1], params.out5[1], tolerance = 3.1)
  expect_equal(params.out4[2], params.out5[2], tolerance = 11.3)
  expect_equal(params.out4[3], params.out5[3], tolerance = 1.6)
  expect_equal(params.out4[4], params.out5[4], tolerance = 0.9)
  expect_equal(params.out4[5], params.out5[5], tolerance = 0.9)
  expect_equal(params.out4[6] - params.out5[6], 0, tolerance = 0.4)
  expect_equal(params.out4[7], params.out5[7], tolerance = 1.5)
  expect_equal(params.out4[8], params.out5[8], tolerance = 0.7)
  expect_equal(params.out4[9], params.out5[9], tolerance = 1.1)
  expect_equal(params.out4[10], params.out5[10], tolerance = 0.3)
  expect_equal(params.out4[11], params.out5[11], tolerance = 0.4)
  expect_equal(params.out4[12] - params.out5[12], 0, tolerance = 0.14)
  expect_equal(params.out4[13] - params.out5[13], 0, tolerance = 0.6)
  expect_equal(params.out4[14] - params.out5[14], 0, tolerance = 0.28)
  expect_equal(params.out4[15] - params.out5[15], 0, tolerance = 0.27)

  # Missing data with lod
  pgp.mmodel6 <- new("MultivariateVecchiaModel", n_neighbors = 25)
  pgp.mmodel6 <- prestogp_fit(pgp.mmodel6, y.list.lod, X.st, locs.list,
    scaling = c(1, 1), common_scale = TRUE, verbose = TRUE, penalty = "SCAD",
    impute.y = TRUE, lod.upper = lod.cut, lod.lower = as.list(rep(0, 3)),
    optim.control = list(trace = 0,
      maxit = 5000, reltol = 1e-3))
  beta.out6 <- get_beta(pgp.mmodel6)
  params.out6 <- pgp.mmodel6@covparams

  # Results should be the same after imputation
  expect_equal(as.numeric(beta.out4[[1]]), as.numeric(beta.out6[[1]]),
    tolerance = 0.02)
  expect_equal(as.numeric(beta.out4[[2]]), as.numeric(beta.out6[[2]]),
    tolerance = 0.02)
  expect_equal(as.numeric(beta.out4[[3]]), as.numeric(beta.out6[[3]]),
    tolerance = 0.02)
  expect_equal(as.numeric(beta.out4[[4]][1] - pgp.mmodel4@Y_bar[1]),
    as.numeric(beta.out6[[4]][1] - pgp.mmodel6@Y_bar[1]), tolerance = 0.02)
  expect_equal(as.numeric(beta.out4[[4]][2] - pgp.mmodel4@Y_bar[2] +
        pgp.mmodel4@Y_bar[1]), as.numeric(beta.out6[[4]][2] -
        pgp.mmodel6@Y_bar[2] + pgp.mmodel6@Y_bar[1]), tolerance = 0.02)
  expect_equal(as.numeric(beta.out4[[4]][3] - pgp.mmodel4@Y_bar[3] +
        pgp.mmodel4@Y_bar[1]), as.numeric(beta.out6[[4]][3] -
        pgp.mmodel6@Y_bar[3] + pgp.mmodel6@Y_bar[1]), tolerance = 0.02)
  expect_equal(params.out4[1], params.out6[1], tolerance = 2.1)
  expect_equal(params.out4[2], params.out6[2], tolerance = 9.4)
  expect_equal(params.out4[3], params.out6[3], tolerance = 0.6)
  expect_equal(params.out4[4] - params.out6[4], 0, tolerance = 0.2)
  expect_equal(params.out4[5], params.out6[5], tolerance = 0.7)
  expect_equal(params.out4[6] - params.out6[6], 0, tolerance = 0.29)
  expect_equal(params.out4[7] - params.out6[7], 0, tolerance = 1.1)
  expect_equal(params.out4[8] - params.out6[8], 0, tolerance = 0.5)
  expect_equal(params.out4[9] - params.out6[9], 0, tolerance = 0.83)
  expect_equal(params.out4[10] - params.out6[10], 0, tolerance = 0.27)
  expect_equal(params.out4[11] - params.out6[11], 0, tolerance = 0.83)
  expect_equal(params.out4[12] - params.out6[12], 0, tolerance = 0.12)
  expect_equal(params.out4[13] - params.out6[13], 0, tolerance = 0.39)
  expect_equal(params.out4[14] - params.out6[14], 0, tolerance = 0.2)
  expect_equal(params.out4[15] - params.out6[15], 0, tolerance = 0.42)
})

test_that("Simulated dataset multivariate spatiotemporal", {
  load("sim_multivariate_big_st2.RData")
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

  expect_gt(beta.out[1], -0.02)
  expect_lt(beta.out[1], 0.05)
  expect_gt(beta.out[2], 0.79)
  expect_lt(beta.out[2], 0.96)
  expect_gt(beta.out[3], 0.94)
  expect_lt(beta.out[3], 1.09)
  expect_gt(beta.out[4], 0.85)
  expect_lt(beta.out[4], 0.91)
  expect_gt(beta.out[5], 1.06)
  expect_lt(beta.out[5], 1.13)
  expect_gt(beta.out[6], -0.001)
  expect_lt(beta.out[6], 0.06)
  expect_gt(beta.out[7], -0.03)
  expect_lt(beta.out[7], 0.014)
  expect_gt(beta.out[8], -0.001)
  expect_lt(beta.out[8], 0.11)
  expect_gt(beta.out[9], -0.08)
  expect_lt(beta.out[9], 0.001)
  expect_gt(beta.out[10], -0.001)
  expect_lt(beta.out[10], 0.06)
  expect_gt(beta.out[11], -0.001)
  expect_lt(beta.out[11], 0.06)
  expect_gt(beta.out[12], -2.41)
  expect_lt(beta.out[12], 2.61)
  expect_gt(beta.out[13], 0.91)
  expect_lt(beta.out[13], 1.11)
  expect_gt(beta.out[14], 0.89)
  expect_lt(beta.out[14], 1.03)
  expect_gt(beta.out[15], 0.95)
  expect_lt(beta.out[15], 1.13)
  expect_gt(beta.out[16], 0.81)
  expect_lt(beta.out[16], 0.92)
  expect_gt(beta.out[17], -0.21)
  expect_lt(beta.out[17], 0.001)
  expect_gt(beta.out[18], -0.02)
  expect_lt(beta.out[18], 0.06)
  expect_gt(beta.out[19], -0.001)
  expect_lt(beta.out[19], 0.12)
  expect_gt(beta.out[20], -0.001)
  expect_lt(beta.out[20], 0.06)
  expect_gt(beta.out[21], -0.06)
  expect_lt(beta.out[21], 0.001)
  expect_gt(beta.out[22], -0.001)
  expect_lt(beta.out[22], 0.015)
  expect_gt(beta.out[23], -1.26)
  expect_lt(beta.out[23], 0.23)
  expect_gt(beta.out[24], 0.9)
  expect_lt(beta.out[24], 1.08)
  expect_gt(beta.out[25], 0.93)
  expect_lt(beta.out[25], 1.32)
  expect_gt(beta.out[26], 0.64)
  expect_lt(beta.out[26], 1.16)
  expect_gt(beta.out[27], 0.88)
  expect_lt(beta.out[27], 0.98)
  expect_gt(beta.out[28], -0.001)
  expect_lt(beta.out[28], 0.06)
  expect_gt(beta.out[29], -0.08)
  expect_lt(beta.out[29], 0.001)
  expect_gt(beta.out[30], -0.14)
  expect_lt(beta.out[30], 0.001)
  expect_gt(beta.out[31], -0.08)
  expect_lt(beta.out[31], 0.001)
  expect_gt(beta.out[32], -0.07)
  expect_lt(beta.out[32], 0.12)
  expect_gt(beta.out[33], -0.001)
  expect_lt(beta.out[33], 0.11)

  expect_gt(params.out[1], 4.5)
  expect_lt(params.out[1], 143.8)
  expect_gt(params.out[2], 3.6)
  expect_lt(params.out[2], 29.6)
  expect_gt(params.out[3], 2.4)
  expect_lt(params.out[3], 29.4)
  expect_gt(params.out[4], 0.07)
  expect_lt(params.out[4], 4.2)
  expect_gt(params.out[5], 7.2)
  expect_lt(params.out[5], 87.7)
  expect_gt(params.out[6], 0.03)
  expect_lt(params.out[6], 2.26)
  expect_gt(params.out[7], 6.7)
  expect_lt(params.out[7], 57.1)
  expect_gt(params.out[8], 0.006)
  expect_lt(params.out[8], 1.09)
  expect_gt(params.out[9], 0.005)
  expect_lt(params.out[9], 80.6)
  expect_gt(params.out[10], 0.42)
  expect_lt(params.out[10], 1.13)
  expect_gt(params.out[11], 0.3)
  expect_lt(params.out[11], 1.93)
  expect_gt(params.out[12], 0.33)
  expect_lt(params.out[12], 2.09)
  expect_gt(params.out[13], 0.02)
  expect_lt(params.out[13], 2.31)
  expect_gt(params.out[14], 0.03)
  expect_lt(params.out[14], 44.9)
  expect_gt(params.out[15], 0.01)
  expect_lt(params.out[15], 4.79)
  expect_gt(params.out[16], -0.99)
  expect_gt(params.out[17], -0.5)
  expect_gt(params.out[18], -0.98)
})

test_that("Invalid locs input for prediction", {
  source("sim_multivariate_small_pred.R")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 5)
  pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
    locs.list.otr,
    scaling = c(1, 1),
    common_scale = TRUE, quiet = TRUE,
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
    common_scale = TRUE, quiet = TRUE,
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
    common_scale = TRUE, quiet = TRUE,
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
    common_scale = TRUE, quiet = TRUE,
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
    common_scale = TRUE, quiet = TRUE,
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
    common_scale = TRUE, quiet = TRUE,
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
    common_scale = TRUE, quiet = TRUE,
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
    common_scale = TRUE, quiet = TRUE,
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
    common_scale = TRUE, quiet = TRUE,
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
    common_scale = TRUE, quiet = TRUE,
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
    common_scale = TRUE, quiet = TRUE,
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
  load("sim_multivariate_big_pred.RData")
  pgp.mmodel1 <- new("MultivariateVecchiaModel", n_neighbors = 25)
  pgp.mmodel1 <- prestogp_fit(pgp.mmodel1, y.list.otr, X.st.otr,
    locs.list.otr,
    scaling = c(1, 1),
    common_scale = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )

  pgp.mmodel1.pred <- prestogp_predict(pgp.mmodel1, X.st.otst, locs.list.otst)

  expect_named(pgp.mmodel1.pred$means, c("Y1", "Y2", "Y3"))

  mse1 <- mean((pgp.mmodel1.pred$means[[1]] - y.list.otst[[1]])^2)
  me1 <- mean(pgp.mmodel1.pred$means[[1]] - y.list.otst[[1]])
  mse2 <- mean((pgp.mmodel1.pred$means[[2]] - y.list.otst[[2]])^2)
  me2 <- mean(pgp.mmodel1.pred$means[[2]] - y.list.otst[[2]])
  mse3 <- mean((pgp.mmodel1.pred$means[[3]] - y.list.otst[[3]])^2)
  me3 <- mean(pgp.mmodel1.pred$means[[3]] - y.list.otst[[3]])
  mse <- mean((unlist(pgp.mmodel1.pred$means) - unlist(y.list.otst))^2)
  me <- mean(unlist(pgp.mmodel1.pred$means) - unlist(y.list.otst))

  expect_gt(mse1, 1.69)
  expect_lt(mse1, 1.85)
  expect_gt(me1, -0.07)
  expect_lt(me1, 0.03)
  expect_gt(mse2, 1.05)
  expect_lt(mse2, 1.21)
  expect_gt(me2, -0.33)
  expect_lt(me2, -0.23)
  expect_gt(mse3, 2.75)
  expect_lt(mse3, 3.63)
  expect_gt(me3, 0.11)
  expect_lt(me3, 0.33)
  expect_gt(mse, 1.85)
  expect_lt(mse, 2.4)
  expect_gt(me, -0.1)
  expect_lt(me, 0.11)

  # SCAD fit
  pgp.mmodel2 <- new("MultivariateVecchiaModel", n_neighbors = 25)
  pgp.mmodel2 <- prestogp_fit(pgp.mmodel2, y.list.otr, X.st.otr,
    locs.list.otr,
    scaling = c(1, 1),
    common_scale = TRUE, quiet = TRUE, penalty = "SCAD",
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )

  pgp.mmodel2.pred <- prestogp_predict(pgp.mmodel2, X.st.otst, locs.list.otst)

  expect_named(pgp.mmodel2.pred$means, c("Y1", "Y2", "Y3"))

  mse1 <- mean((pgp.mmodel2.pred$means[[1]] - y.list.otst[[1]])^2)
  me1 <- mean(pgp.mmodel2.pred$means[[1]] - y.list.otst[[1]])
  mse2 <- mean((pgp.mmodel2.pred$means[[2]] - y.list.otst[[2]])^2)
  me2 <- mean(pgp.mmodel2.pred$means[[2]] - y.list.otst[[2]])
  mse3 <- mean((pgp.mmodel2.pred$means[[3]] - y.list.otst[[3]])^2)
  me3 <- mean(pgp.mmodel2.pred$means[[3]] - y.list.otst[[3]])
  mse <- mean((unlist(pgp.mmodel2.pred$means) - unlist(y.list.otst))^2)
  me <- mean(unlist(pgp.mmodel2.pred$means) - unlist(y.list.otst))

  expect_gt(mse1, 1.69)
  expect_lt(mse1, 1.82)
  expect_gt(me1, -0.06)
  expect_lt(me1, 0.06)
  expect_gt(mse2, 1.03)
  expect_lt(mse2, 1.09)
  expect_gt(me2, -0.29)
  expect_lt(me2, -0.23)
  expect_gt(mse3, 2.75)
  expect_lt(mse3, 3.07)
  expect_gt(me3, 0.19)
  expect_lt(me3, 0.31)
  expect_gt(mse, 1.84)
  expect_lt(mse, 1.99)
  expect_gt(me, -0.03)
  expect_lt(me, 0.02)
})
