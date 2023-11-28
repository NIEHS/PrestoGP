context("CreateU Multivariate")

test_that("create.param.sequence", {
  seq = create.param.sequence(1)
  colnames(seq) <- NULL
  expect_equal(2, ncol(seq))
  expect_equal(5, nrow(seq))
  expect_equal(c(1,1), seq[1,])
  expect_equal(c(2,2), seq[2,])
  expect_equal(c(3,3), seq[3,])
  expect_equal(c(4,4), seq[4,])
  expect_equal(c(5,5), seq[5,])

  seq = create.param.sequence(3)
  colnames(seq) <- NULL
  expect_equal(2, ncol(seq))
  expect_equal(5, nrow(seq))
  expect_equal(c(1,3), seq[1,])
  expect_equal(c(4,6), seq[2,])
  expect_equal(c(7,9), seq[3,])
  expect_equal(c(10,12), seq[4,])
  expect_equal(c(13,15), seq[5,])

  seq = create.param.sequence(3, 2)
  colnames(seq) <- NULL
  expect_equal(2, ncol(seq))
  expect_equal(5, nrow(seq))
  expect_equal(c(1,3), seq[1,])
  expect_equal(c(4,9), seq[2,])
  expect_equal(c(10,12), seq[3,])
  expect_equal(c(13,15), seq[4,])
  expect_equal(c(16,18), seq[5,])
})

test_that("max_min_ordering", {
  set.seed(7919)
  load("multivariate_sim_spatial3.Rdata")
  order <- max_min_ordering(locs_train, fields::rdist)
  exp_order <- c(63, 21, 50, 30, 76, 78, 40, 36, 23, 69, 9, 67, 32, 2, 20, 62,
                 22, 31, 74, 39, 35, 58, 68, 54, 41, 3, 80, 46, 6, 88, 12, 47,
                 72, 42, 13, 83, 25, 52, 11, 60, 24, 1, 28, 84, 29, 64, 66, 81,
                 82, 55, 61, 87, 17, 33, 43, 45, 10, 79, 53, 75, 89, 51, 73, 27,
                 26, 77, 44, 38, 65, 16, 19, 37, 57, 70, 15, 4, 5, 86, 14, 49,
                 85, 34, 48, 59, 18, 8, 71, 7, 56, 90)
  order.gpv <- order_maxmin_exact(locs_train)
  expect_equal(exp_order, order)
  expect_equal(order, order.gpv)
})

test_that("knn_indices", {
  set.seed(7919)
  load("multivariate_sim_spatial3.Rdata")
  order <- max_min_ordering(locs_train, fields::rdist)
  ordered_locs <- locs_train[order,,drop=FALSE]
  indices <- knn_indices(ordered_locs, ordered_locs[2,,drop=FALSE], 5, fields::rdist, "rdist")
  expect_equal(c(0.00, 0.00148, 0.0127, 0.243, 0.252, 0.254), indices$distances, tolerance=10e-2)
  expect_equal(c(2, 87, 28, 32, 17, 33), indices$indices)
})

test_that("knn_indices_small_kd_node", {
  set.seed(7919)
  load("multivariate_sim_spatial3.Rdata")
  ordered_locs <- c(c(0.1,0.1), c(0.9,0.9), c(0.9,0.1), c(0.1,0.9))
  ordered_locs <- matrix(data=ordered_locs, ncol=2, nrow=4)
  indices <- knn_indices(ordered_locs, ordered_locs[2,,drop=FALSE], 3, fields::rdist, "rdist")
  expect_equal(c(2, 1, 3, 4), indices$indices)
})
