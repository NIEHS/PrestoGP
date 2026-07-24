#test_that("revMat", {
#  foo <- matrix(1:12, nrow = 4, byrow = TRUE)
#  bar <- matrix(12:1, nrow = 4, byrow = TRUE)
#  expect_equal(revMat(foo), bar)
#})

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

test_that("U2V_cpp", {
  set.seed(1212)

  dim1 <- (0:9)^2 + rnorm(10, 0, 1e-2)
  dim2 <- (1:10)^2 + rnorm(10, 0, 1e-2)

  locs <- as.matrix(expand.grid(dim1, dim2))

  params <- c(3, 1.5, 0.6, 2)

  vec.approx <- vecchia_specify(locs, m = 5)
  U.obj <- createU(vec.approx, covparms = params[1:3], nuggets = params[4])
  vec.mapprox <- vecchia_Mspecify(list(locs), 5)
  U.mobj <- createUMultivariate(vec.mapprox, c(params, 1))

  expect_equal(sum(abs(GPvecchia:::U2V(U.obj) - PrestoGP:::U2V_cpp(U.obj))), 0,
    tolerance = 1e-4)
  expect_equal(sum(abs(GPvecchia:::U2V(U.mobj) - PrestoGP:::U2V_cpp(U.mobj))),
    0, tolerance = 1e-4)

  source("sim_multivariate.R")

  vec.mapprox <- vecchia_Mspecify(locs.list, 25)
  U.mobj <- createUMultivariate(vec.mapprox, c(
    marg.var, ranges,
    marg.smoothness,
    nuggets, rho.vec
  ))
  expect_equal(sum(abs(GPvecchia:::U2V(U.mobj) - PrestoGP:::U2V_cpp(U.mobj))),
    0, tolerance = 1e-4)
})
