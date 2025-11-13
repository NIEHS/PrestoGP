test_that("revMat", {
  foo <- matrix(1:12, nrow = 4, byrow = TRUE)
  bar <- matrix(12:1, nrow = 4, byrow = TRUE)
  expect_equal(revMat(foo), bar)
})

test_that("MMatern_cov univariate", {
  source("sim_vecchia.R")

  Sigma.All <- Sigma.All + params.all[4] * diag(ncol(Sigma.All))

  locs.nn <- nn2(locs, k = 25)$nn.idx

  ndx <- sample(seq_len(nrow(locs.nn)), size = 10)

  Sigma.hat <- array(dim = c(ncol(locs.nn), ncol(locs.nn), 10))
  for (i in 1:10) {
    Sigma.hat[, , i] <- MMatern_cov(locs[locs.nn[ndx[i], ], ],
      rep(1, ncol(locs.nn)), c(params.all, 1), 1)
    expect_lt(sum(abs(Sigma.hat[, , i] - Sigma.All[locs.nn[ndx[i], ],
            locs.nn[ndx[i], ]])), 1e-4)
  }
})

test_that("MMatern_cov multivariate", {
  source("sim_multivariate_big.R")

  locs <- NULL
  for (i in 1:3) {
    locs <- rbind(locs, locs.list[[i]])
  }

  npy <- n.spatial.xy^2
  nuggetv <- c(rep(nuggets[1], npy), rep(nuggets[2], npy), rep(nuggets[3], npy))
  y.ndx <- c(rep(1, npy), rep(2, npy), rep(3, npy))
  Sigma.All <- Sigma.All + nuggetv * diag(ncol(Sigma.All))

  locs.nn <- nn2(locs, k = 25)$nn.idx

  ndx <- sample(seq_len(nrow(locs.nn)), size = 10)
  Sigma.hat <- array(dim = c(ncol(locs.nn), ncol(locs.nn), 10))
  for (i in 1:10) {
    Sigma.hat[, , i] <- PrestoGP:::MMatern_cov(locs[locs.nn[ndx[i], ], ],
      y.ndx[locs.nn[ndx[i], ]], params.all, 3)
    expect_lt(sum(abs(Sigma.hat[, , i] - Sigma.All[locs.nn[ndx[i], ],
            locs.nn[ndx[i], ]])), 1e-4)
  }
})
