context("Multivariate Spatial Model")

test_that("Invalid predict Locs", {
  load("small_sim.Rdata")
  model <- new("MultivariateSpatialModel")
  expect_error(prestogp_predict(model, X_test, "locs_test"),
               "The locs parameter must be a matrix.")
})

test_that("Invalid predict X", {
  load("small_sim.Rdata")
  model <- new("MultivariateSpatialModel")
  expect_error(prestogp_predict(model, "X_test", locs_test),
               "X parameter must be a matrix.")
})

test_that("Invalid predict locs (not 2 columns)", {
  load("small_sim.Rdata")
  model <- new("MultivariateSpatialModel")
  expect_error(prestogp_predict(model, matrix(rnorm(100), ncol=10), matrix(rnorm(30), ncol=3)),
               "The locs parameter must have 2 columns.")
})

test_that("locs length mismatch", {
  load("small_sim.Rdata")
  model <- new("MultivariateSpatialModel")
  expect_error(prestogp_predict(model, matrix(rnorm(100), ncol=10), matrix(rnorm(50), ncol=2)),
               "The number of locations must match the number of X observations.")
})

test_that("Simulated univariate dataset spatial", {
  set.seed(7919)
  load("sim_spatial.Rdata")
  model <- new("MultivariateSpatialModel")
  locs_test = locs_test[1:100,]
  locs_train = locs_train[1:100,]
  X_test = X_test[1:100,]
  X_train = X_train[1:100,]
  Y_test = Y_test[1:100]
  Y_train = Y_train[1:100]
  model <- prestogp_fit(model, Y_train, X_train, locs_train)
  prediction <- prestogp_predict(model, X_test, locs_test)
  means <- prediction[[1]]
  mean_sds <- mean(prediction[[2]])
  mse <- mean((Y_test-means)^2)
  expect_equal(4.45, mse, tolerance=10e-2)
  expect_equal(0.11455, mean_sds, tolerance=10e-4)
  expect_equal(11.4, model@covparams[1], tolerance=10e-2)
  expect_equal(59.8, model@covparams[2], tolerance=10e-2)
  expect_equal(0.296, model@covparams[3], tolerance=10e-2)
  expect_equal(0.0016, model@covparams[4], tolerance=10e-4)
})



test_that("Simulated dataset spatial", {
  set.seed(7919)
  load("multivariate_sim_spatial.Rdata")
  model <- new("MultivariateSpatialModel")
  locs_test = locs_test[1:100,]
  locs_train = locs_train[1:100,]
  X_test = X_test[1:100,]
  X_train = X_train[1:100,]
  Y_test = Y_test[1:100]
  Y_train = Y_train[1:100]
  model <- prestogp_fit(model, Y_train, X_train, locs_train)
  prediction <- prestogp_predict(model, X_test, locs_test)
  means <- prediction[[1]]
  mean_sds <- mean(prediction[[2]])
  mse <- mean((Y_test-means)^2)
  expect_equal(4.45, mse, tolerance=10e-2)
  expect_equal(0.11455, mean_sds, tolerance=10e-4)
  expect_equal(11.4, model@covparams[1], tolerance=10e-2)
  expect_equal(59.8, model@covparams[2], tolerance=10e-2)
  expect_equal(0.296, model@covparams[3], tolerance=10e-2)
  expect_equal(0.0016, model@covparams[4], tolerance=10e-4)
})