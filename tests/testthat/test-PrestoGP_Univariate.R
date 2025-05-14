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

test_that("NA's in X", {
  model <- new("VecchiaModel")
  expect_error(
    prestogp_fit(
      model, as.matrix(1:4), as.matrix(c(1:3, NA)),
      as.matrix(1:4)
    ),
    "X must not contain NA's"
  )
})

test_that("NA's in locs", {
  model <- new("VecchiaModel")
  expect_error(
    prestogp_fit(
      model, as.matrix(1:4), as.matrix(1:4),
      as.matrix(c(1:3, NA))
    ),
    "locs must not contain NA's"
  )
})

test_that("NA's in Y", {
  model <- new("VecchiaModel")
  expect_error(
    prestogp_fit(
      model, as.matrix(c(1:3, NA)), as.matrix(1:4),
      as.matrix(1:4)
    ),
    "Y contains NA's and impute.y is FALSE. Set impute.y=TRUE to impute missing Y's."
  )
})

test_that("lod.upper not numeric", {
  model <- new("VecchiaModel")
  expect_error(
    prestogp_fit(
      model, as.matrix(c(1:4)), as.matrix(1:4),
      as.matrix(1:4), impute.y = TRUE, lod.upper = "foo",
    ),
    "lod.upper must be numeric"
  )
})

test_that("lod.lower not numeric", {
  model <- new("VecchiaModel")
  expect_error(
    prestogp_fit(
      model, as.matrix(c(1:4)), as.matrix(1:4),
      as.matrix(1:4), impute.y = TRUE, lod.lower = "foo",
    ),
    "lod.lower must be numeric"
  )
})

test_that("length(lod.upper) != nrow(X)", {
  model <- new("VecchiaModel")
  expect_error(
    prestogp_fit(
      model, as.matrix(c(1:4)), as.matrix(1:4),
      as.matrix(1:4), impute.y = TRUE, lod.upper = 1:2,
    ),
    "Length of lod.upper must equal the number of observations"
  )
})

test_that("length(lod.lower) != nrow(X)", {
  model <- new("VecchiaModel")
  expect_error(
    prestogp_fit(
      model, as.matrix(c(1:4)), as.matrix(1:4),
      as.matrix(1:4), impute.y = TRUE, lod.lower = 1:2,
    ),
    "Length of lod.lower must equal the number of observations"
  )
})

test_that("length(Y.names) != 1", {
  model <- new("VecchiaModel")
  expect_error(
    prestogp_fit(
      model, as.matrix(c(1:4)), as.matrix(1:4),
      as.matrix(1:4), Y.names = c("test1", "test2"),
    ),
    "Length of Y.names must match the number of response variables"
  )
})

test_that("length(X.names) != ncol(X)", {
  model <- new("VecchiaModel")
  expect_error(
    prestogp_fit(
      model, as.matrix(c(1:4)), as.matrix(1:4),
      as.matrix(1:4), X.names = c("test1", "test2"),
    ),
    "Length of X.names must match the number of predictor variables"
  )
})

test_that("Simulated dataset spatial", {
  load("sim_vecchia.RData")
  pgp.model1 <- new("VecchiaModel", n_neighbors = 25)
  pgp.model1 <- prestogp_fit(pgp.model1, y, X, locs, Y.names = "test",
    X.names = c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"),
    scaling = c(1, 1), common_scale = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )

  expect_silent(plot_beta(pgp.model1))
  dev.off()

  expect_true(validObject(pgp.model1))
  show(pgp.model1)

  expect_equal(get_Y(pgp.model1), y)
  expect_equal(get_linear_model(pgp.model1), pgp.model1@linear_model)
  expect_equal(get_neighbors(pgp.model1), pgp.model1@n_neighbors)
  expect_equal(get_scaling(pgp.model1), pgp.model1@scaling)
  expect_equal(get_converged(pgp.model1), pgp.model1@converged)
  expect_equal(get_pen_loglik(pgp.model1), pgp.model1@error)

  beta.out <- get_beta(pgp.model1)
  params.out <- pgp.model1@covparams
  theta.out <- get_theta(pgp.model1)

  expect_length(beta.out, 2)
  expect_length(beta.out[[1]], 10)
  expect_length(beta.out[[2]], 1)
  expect_length(params.out, 5)
  expect_length(theta.out[[1]], 1)
  expect_length(theta.out[[2]], 1)
  expect_length(theta.out[[3]], 1)
  expect_length(theta.out[[4]], 1)

  expect_named(theta.out, c("sigma", "scale", "smoothness", "nuggets"))
  expect_named(theta.out[[1]], "test")
  expect_named(theta.out[[2]], "test")
  expect_named(theta.out[[3]], "test")
  expect_named(theta.out[[4]], "test")
  expect_named(beta.out, c("test", "(Intercept)"))
  expect_named(beta.out[[1]], c("x1", "x2", "x3", "x4", "x5", "x6", "x7",
      "x8", "x9", "x10"))
  expect_named(beta.out[[2]], "(Intercept)")

  expect_identical(as.numeric(theta.out[[1]]), params.out[1])
  expect_identical(as.numeric(theta.out[[2]]), params.out[2])
  expect_identical(as.numeric(theta.out[[3]]), params.out[3])
  expect_identical(as.numeric(theta.out[[4]]), params.out[4])
  expect_equal(as.numeric(beta.out[[2]] - mean(y)),
    as.numeric(pgp.model1@beta[1]))
  expect_identical(as.vector(beta.out[[1]]), as.vector(pgp.model1@beta[2:11]))

  expect_equal(as.numeric(beta.out[[1]]), c(0.96, 1.03, 1.0, 1.0, rep(0, 6)),
    tolerance = 0.02)
  expect_equal(as.numeric(beta.out[[2]] - mean(y)), 0.002, tolerance = 0.003)
  expect_gt(params.out[1], 1.28)
  expect_lt(params.out[1], 1.65)
  expect_gt(params.out[2], 0.35)
  expect_lt(params.out[2], 0.46)
  expect_gt(params.out[3], 0.52)
  expect_lt(params.out[3], 0.62)
  expect_gt(params.out[4], 1.9)
  expect_lt(params.out[4], 2.005)

  pgp.model2 <- new("FullModel")
  pgp.model2 <- prestogp_fit(pgp.model2, y, X, locs,
    scaling = c(1, 1), common_scale = TRUE, quiet = TRUE, verbose = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )

  expect_true(validObject(pgp.model2))

  beta.out2 <- get_beta(pgp.model2)
  params.out2 <- pgp.model2@covparams

  expect_length(beta.out2, 2)
  expect_length(beta.out2[[1]], 10)
  expect_length(beta.out2[[2]], 1)
  expect_length(params.out2, 5)

  expect_named(beta.out2, c("Y", "(Intercept)"))

  expect_equal(beta.out2[[1]], c(0.96, 1.03, 1.005, 0.99, rep(0, 6)),
    tolerance = 0.04)
  expect_equal(as.numeric(beta.out2[[2]] - mean(y)) - 0.0137, 0,
    tolerance = 0.004)
  expect_gt(params.out2[1], 1.50)
  expect_lt(params.out2[1], 1.51)
  expect_gt(params.out2[2], 0.35)
  expect_lt(params.out2[2], 0.40)
  expect_gt(params.out2[3], 0.574)
  expect_lt(params.out2[3], 0.575)
  expect_gt(params.out2[4], 1.92)
  expect_lt(params.out2[4], 1.93)

  # Vecchia and full models should be approximately equal
  expect_equal(as.numeric(beta.out[[2]]), as.numeric(beta.out2[[2]]),
    tolerance = 0.03)
  expect_equal(as.numeric(beta.out[[1]]), as.numeric(beta.out2[[1]]),
    tolerance = 0.02)
  expect_equal(params.out[1], params.out2[1], tolerance = 0.2)
  expect_equal(params.out[2] - params.out2[2], 0, tolerance = 0.15)
  expect_equal(params.out[3] - params.out2[3], 0, tolerance = 0.05)
  expect_equal(params.out[4], params.out2[4], tolerance = 0.1)

  # Missing data
  set.seed(1234)
  y.na <- y
  y.na[sample(seq_along(y), floor(0.1 * length(y)))] <- NA

  pgp.model3 <- new("VecchiaModel", n_neighbors = 25)
  pgp.model3 <- prestogp_fit(pgp.model3, y.na, X, locs,
    scaling = c(1, 1), common_scale = TRUE, quiet = FALSE,
    impute.y = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  beta.out3 <- as.vector(pgp.model3@beta)
  params.out3 <- pgp.model3@covparams

  # Results should be the same after imputation
  expect_equal(as.numeric(beta.out[[2]] - mean(y)), beta.out3[1],
    tolerance = 0.005)
  expect_equal(as.numeric(beta.out[[1]]), beta.out3[-1], tolerance = 0.06)
  expect_equal(params.out[1], params.out3[1], tolerance = 0.9)
  expect_equal(params.out[2] - params.out3[2], 0, tolerance = 0.3)
  expect_equal(params.out[3] - params.out3[3], 0, tolerance = 0.1)
  expect_equal(params.out[4], params.out3[4], tolerance = 0.3)

  # Missing data with lod
  y.lod <- y + 10
  lod.cut <- quantile(y.lod, 0.1)
  y.na.lod <- y.lod
  y.na.lod[y.na.lod <= lod.cut] <- NA
  lodupper <- rep(lod.cut, length(y))
  lodlower <- rep(0, length(y))

  doParallel::registerDoParallel(cores = 2)
  pgp.model4 <- new("VecchiaModel", n_neighbors = 25)
  pgp.model4 <- prestogp_fit(pgp.model4, y.na.lod, X, locs,
    scaling = c(1, 1), common_scale = TRUE, verbose = TRUE, parallel = TRUE,
    impute.y = TRUE, lod.upper = lodupper, lod.lower = lodlower,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  beta.out4 <- as.vector(pgp.model4@beta)
  params.out4 <- pgp.model4@covparams

  # Results should be the same after imputation
  expect_equal(as.numeric(beta.out[[2]] - mean(y)) - beta.out4[1], 0,
    tolerance = 0.004)
  expect_equal(as.numeric(beta.out[[1]]), beta.out4[-1], tolerance = 0.04)
  expect_equal(params.out[1] - params.out4[1], 0, tolerance = 0.9)
  expect_equal(params.out[2] - params.out4[2], 0, tolerance = 0.3)
  expect_equal(params.out[3], params.out4[3], tolerance = 0.4)
  expect_equal(params.out[4], params.out4[4], tolerance = 0.3)

  # Adaptive lasso fit
  pgp.model5 <- new("VecchiaModel", n_neighbors = 25)
  pgp.model5 <- prestogp_fit(pgp.model5, y, X, locs,
    scaling = c(1, 1), common_scale = TRUE, quiet = TRUE, adaptive = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  beta.out5 <- get_beta(pgp.model5)
  params.out5 <- pgp.model5@covparams

  # Adaptive and non-adaptive fits should be approximately equal
  expect_equal(as.numeric(beta.out[[2]]), as.numeric(beta.out5[[2]]),
    tolerance = 0.005)
  expect_equal(as.numeric(beta.out[[1]]), as.numeric(beta.out5[[1]]),
    tolerance = 0.06)
  expect_equal(params.out[1], params.out5[1], tolerance = 0.4)
  expect_equal(params.out[2] - params.out5[2], 0, tolerance = 0.06)
  expect_equal(params.out[3], params.out5[3], tolerance = 0.1)
  expect_equal(params.out[4], params.out5[4], tolerance = 0.1)

  # Relaxed lasso fit
  pgp.model6 <- new("VecchiaModel", n_neighbors = 25)
  pgp.model6 <- prestogp_fit(pgp.model6, y, X, locs,
    scaling = c(1, 1), common_scale = TRUE, quiet = TRUE, penalty = "relaxed",
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  beta.out6 <- get_beta(pgp.model6)
  params.out6 <- pgp.model6@covparams

  expect_equal(as.numeric(beta.out6[[1]]), c(0.979, 1.037, 1.009, 1.015,
      rep(0, 6)), tolerance = 0.003)
  expect_equal(as.numeric(beta.out6[[2]] - mean(y)), 0.002, tolerance = 0.003)
  expect_equal(params.out6[1], 1.44, tolerance = 0.2)
  expect_equal(params.out6[2] - 0.41, 0, tolerance = 0.06)
  expect_equal(params.out6[3] - 0.56, 0, tolerance = 0.06)
  expect_equal(params.out6[4], 1.99, tolerance = 0.1)

  # SCAD fit
  pgp.model7 <- new("VecchiaModel", n_neighbors = 25)
  pgp.model7 <- prestogp_fit(pgp.model7, y, X, locs,
    scaling = c(1, 1), common_scale = TRUE, quiet = TRUE, penalty = "SCAD",
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )

  expect_silent(plot_beta(pgp.model7))
  dev.off()

  beta.out7 <- get_beta(pgp.model7)
  params.out7 <- pgp.model7@covparams

  expect_equal(as.numeric(beta.out7[[1]]), c(0.979, 1.038, 1.009, 1.015,
      rep(0, 6)), tolerance = 0.006)
  expect_equal(as.numeric(beta.out7[[2]] - mean(y)), 0.002, tolerance = 0.003)
  expect_equal(params.out7[1], 1.47, tolerance = 0.2)
  expect_equal(params.out7[2] - 0.405, 0, tolerance = 0.06)
  expect_equal(params.out7[3] - 0.56, 0, tolerance = 0.06)
  expect_equal(params.out7[4], 2.0, tolerance = 0.1)

  # Missing data
  pgp.model8 <- new("VecchiaModel", n_neighbors = 25)
  pgp.model8 <- prestogp_fit(pgp.model8, y.na, X, locs,
    scaling = c(1, 1), common_scale = TRUE, quiet = FALSE,
    impute.y = TRUE, penalty = "SCAD",
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  beta.out8 <- as.vector(pgp.model8@beta)
  params.out8 <- pgp.model8@covparams

  # Results should be the same after imputation
  expect_equal(as.numeric(beta.out7[[2]] - mean(y)), beta.out8[1],
    tolerance = 0.006)
  expect_equal(as.numeric(beta.out7[[1]]), beta.out8[-1], tolerance = 0.06)
  expect_equal(params.out7[1], params.out8[1], tolerance = 0.7)
  expect_equal(params.out7[2] - params.out8[2], 0, tolerance = 0.3)
  expect_equal(params.out7[3] - params.out8[3], 0, tolerance = 0.1)
  expect_equal(params.out7[4], params.out8[4], tolerance = 0.3)

  # Missing data with lod
  cl <- parallel::makeCluster(2)
  pgp.model9 <- new("VecchiaModel", n_neighbors = 25)
  pgp.model9 <- prestogp_fit(pgp.model9, y.na.lod, X, locs,
    scaling = c(1, 1), common_scale = TRUE, verbose = TRUE,
    impute.y = TRUE, lod.upper = lod.cut, lod.lower = 0, penalty = "SCAD",
    cluster = cl, optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  parallel::stopCluster(cl)
  beta.out9 <- as.vector(pgp.model9@beta)
  params.out9 <- pgp.model9@covparams

  # Results should be the same after imputation
  expect_equal(as.numeric(beta.out7[[2]] - mean(y)) - beta.out9[1], 0,
    tolerance = 0.004)
  expect_equal(as.numeric(beta.out7[[1]]), beta.out9[-1], tolerance = 0.02)
  expect_equal(params.out7[1] - params.out9[1], 0, tolerance = 1.0)
  expect_equal(params.out7[2] - params.out9[2], 0, tolerance = 0.2)
  expect_equal(params.out7[3], params.out9[3], tolerance = 0.4)
  expect_equal(params.out7[4], params.out9[4], tolerance = 0.3)

  # MCP fit
  pgp.model10 <- new("VecchiaModel", n_neighbors = 25)
  pgp.model10 <- prestogp_fit(pgp.model10, y, X, locs,
    scaling = c(1, 1), common_scale = TRUE, quiet = TRUE, penalty = "MCP",
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  beta.out10 <- get_beta(pgp.model10)
  params.out10 <- pgp.model10@covparams

  expect_equal(as.numeric(beta.out10[[1]]), c(0.979, 1.038, 1.009, 1.015,
      rep(0, 6)), tolerance = 0.03)
  expect_equal(as.numeric(beta.out10[[2]] - mean(y)), 0.002, tolerance = 0.003)
  expect_equal(params.out10[1], 1.46, tolerance = 0.2)
  expect_equal(params.out10[2] - 0.41, 0, tolerance = 0.05)
  expect_equal(params.out10[3] - 0.565, 0, tolerance = 0.06)
  expect_equal(params.out10[4], 2.0, tolerance = 0.1)
})

test_that("Simulated dataset spatiotemporal", {
  load("sim_vecchia_st.RData")
  Y <- as.matrix(y)
  colnames(Y) <- "test"
  pgp.model1 <- new("VecchiaModel", n_neighbors = 25)
  pgp.model1 <- prestogp_fit(pgp.model1, Y, X, locs,
    scaling = c(1, 1, 2),
    common_scale = FALSE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )

  expect_true(validObject(pgp.model1))

  beta.out <- as.vector(pgp.model1@beta)
  params.out <- pgp.model1@covparams
  theta.out <- get_theta(pgp.model1)

  expect_length(beta.out, 11)
  expect_length(params.out, 6)
  expect_length(theta.out[[1]], 1)
  expect_length(theta.out[[2]], 2)
  expect_length(theta.out[[3]], 1)
  expect_length(theta.out[[4]], 1)

  expect_named(theta.out, c("sigma", "scale", "smoothness", "nuggets"))
  expect_named(theta.out[[1]], "test")
  expect_named(theta.out[[2]], c("test_1", "test_2"))
  expect_named(theta.out[[3]], "test")
  expect_named(theta.out[[4]], "test")

  expect_identical(as.numeric(theta.out[[1]]), params.out[1])
  expect_identical(as.numeric(theta.out[[2]][1]), params.out[2])
  expect_identical(as.numeric(theta.out[[2]][2]), params.out[3])
  expect_identical(as.numeric(theta.out[[3]]), params.out[4])
  expect_identical(as.numeric(theta.out[[4]]), params.out[5])

  expect_equal(beta.out, c(0.002, 0.956, 1.023, 0.961, 1.03, 0, 0.011, 0.056,
      -0.04, -0.026, 0), tolerance = 0.006)
  expect_gt(params.out[1], 0.87)
  expect_lt(params.out[1], 1.2)
  expect_gt(params.out[2], 0.09)
  expect_lt(params.out[2], 0.15)
  expect_gt(params.out[3], 0.9)
  expect_lt(params.out[3], 1.5)
  expect_gt(params.out[4], 1)
  expect_lt(params.out[4], 1.25)
  expect_gt(params.out[5], 0.7)
  expect_lt(params.out[5], 0.77)

  pgp.model2 <- new("FullModel", n_neighbors = 25)
  pgp.model2 <- prestogp_fit(pgp.model2, y, X, locs,
    scaling = c(1, 1, 2),
    common_scale = FALSE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )

  expect_true(validObject(pgp.model2))

  beta.out2 <- as.vector(pgp.model2@beta)
  params.out2 <- pgp.model2@covparams

  expect_length(beta.out2, 11)
  expect_length(params.out2, 6)
  expect_equal(beta.out2, c(0.011, 0.956, 1.023, 0.958, 1.029, 0, 0.013, 0.059,
      -0.04, -0.027, 0), tolerance = 0.02)
  expect_gt(params.out2[1], 1)
  expect_lt(params.out2[1], 1.03)
  expect_gt(params.out2[2], 0.0875)
  expect_lt(params.out2[2], 0.105)
  expect_gt(params.out2[3], 1.1)
  expect_lt(params.out2[3], 1.4)
  expect_gt(params.out2[4], 1.1)
  expect_lt(params.out2[4], 1.4)
  expect_gt(params.out2[5], 0.74)
  expect_lt(params.out2[5], 0.77)

  # Vecchia and full models should be approximately equal
  expect_equal(beta.out, beta.out2, tolerance = 0.02)
  expect_equal(params.out[1], params.out2[1], tolerance = 0.15)
  expect_equal(params.out[2] - params.out2[2], 0, tolerance = 0.06)
  expect_equal(params.out[3] - params.out2[3], 0, tolerance = 0.4)
  expect_equal(params.out[4], params.out2[4], tolerance = 0.4)
  expect_equal(params.out[5], params.out2[5], tolerance = 0.06)
})

test_that("Invalid locs input for prediction", {
  source("sim_vecchia_small_pred.R")
  pgp.model1 <- new("VecchiaModel", n_neighbors = 5)
  pgp.model1 <- prestogp_fit(pgp.model1, y.otr, X.otr,
    locs.otr,
    scaling = c(1, 1),
    common_scale = TRUE, quiet = TRUE,
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
    common_scale = TRUE, quiet = TRUE,
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
    common_scale = TRUE, quiet = TRUE,
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
    common_scale = TRUE, quiet = TRUE,
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
    common_scale = TRUE, quiet = TRUE,
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
    common_scale = TRUE, quiet = TRUE,
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

test_that("full prediction not implemented", {
  source("sim_vecchia_small_pred.R")
  pgp.model1 <- new("FullModel", n_neighbors = 5)
  pgp.model1 <- prestogp_fit(pgp.model1, y.otr, X.otr,
    locs.otr,
    scaling = c(1, 1),
    common_scale = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )
  expect_error(
    prestogp_predict(pgp.model1, X.otst, locs.otst),
    "Prediction is not currently supported for full models"
  )
})

test_that("m too large for prediction", {
  source("sim_vecchia_small_pred.R")
  pgp.model1 <- new("VecchiaModel", n_neighbors = 5)
  pgp.model1 <- prestogp_fit(pgp.model1, y.otr, X.otr,
    locs.otr,
    scaling = c(1, 1),
    common_scale = TRUE, quiet = TRUE,
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
  load("sim_vecchia_pred.RData")
  pgp.model1 <- new("VecchiaModel", n_neighbors = 25)
  pgp.model1 <- prestogp_fit(pgp.model1, y.otr, X.otr, locs.otr,
    scaling = c(1, 1),
    common_scale = TRUE, quiet = TRUE,
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )

  pgp.model1.pred <- prestogp_predict(pgp.model1, X.otst, locs.otst)

  mse <- mean((pgp.model1.pred$means - y.otst)^2)
  me <- mean(pgp.model1.pred$means - y.otst)

  expect_equal(mse, 2.1733, tolerance = 0.0015)
  expect_equal(me - 0.0442, 0, tolerance = 0.002)

  # SCAD fit
  pgp.model2 <- new("VecchiaModel", n_neighbors = 25)
  pgp.model2 <- prestogp_fit(pgp.model2, y.otr, X.otr, locs.otr,
    scaling = c(1, 1),
    common_scale = TRUE, quiet = TRUE, penalty = "SCAD",
    optim.control = list(
      trace = 0, maxit = 5000,
      reltol = 1e-3
    )
  )

  pgp.model2.pred <- prestogp_predict(pgp.model2, X.otst, locs.otst)

  mse <- mean((pgp.model2.pred$means - y.otst)^2)
  me <- mean(pgp.model2.pred$means - y.otst)

  expect_equal(mse, 2.1655, tolerance = 0.001)
  expect_equal(me - 0.04496, 0, tolerance = 0.002)
})
