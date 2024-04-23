test_that("negloglik_vecchia", {
  set.seed(1234)

  y.orig <- rnorm(100)
  locs.dim1 <- seq(1, 10, length = 10) + rnorm(10, 0, 0.1)
  locs.dim2 <- seq(1, 10, length = 10) + rnorm(10, 0, 0.1)
  locs <- as.matrix(expand.grid(locs.dim1, locs.dim2))
  covmat.true <- Matern(rdist(locs), range = 1, smoothness = 0.5) + diag(100)

  y <- y.orig %*% chol(covmat.true)
  y <- as.vector(y)

  logparams <- create.initial.values.flex(0.9, 0.5, 0.5, 0.1, 1, 1)
  pseq <- create.param.sequence(1)
  vec.approx <- vecchia_specify(locs, 5)
  LL.vecchia <- negloglik_vecchia(logparams, y, vec.approx, pseq)
  expect_equal(194.75, LL.vecchia, tolerance = 1e-2)

  vec.mapprox <- vecchia_Mspecify(list(locs), 5)
  LL.mvecchia <- mvnegloglik(logparams, vec.mapprox, y, pseq, 1)
  # Univariate likelihood should equal multivariate likelihood
  expect_equal(LL.vecchia, LL.mvecchia, tolerance = 1e-1)

  logparams2 <- create.initial.values.flex(0.9, 1, 0.5, 0.1, 1, 1)
  vec.approx2 <- vec.approx
  vec.approx2$locsord <- vec.approx$locsord / 0.5
  LL.vecchia2 <- negloglik_vecchia(logparams2, y, vec.approx2, pseq)
  LL.vecchia.st <- negloglik_vecchia_ST(logparams, y, vec.approx, pseq, 1, 1)
  vec.mapprox2 <- vec.mapprox
  vec.mapprox2$locsord <- vec.mapprox$locsord / 0.5
  LL.mvecchia2 <- mvnegloglik_ST(logparams2, vec.mapprox2, y, pseq, 1, 1, 1)

  # Likelihood should equal spatiotemporal likelihood after scaling locs
  expect_equal(LL.vecchia, LL.vecchia2, tolerance = 1e-3)
  expect_equal(LL.vecchia, LL.vecchia.st, tolerance = 1e-3)
  expect_equal(LL.vecchia, LL.mvecchia2, tolerance = 1e-1)
})

test_that("negloglik_vecchia_ST", {
  set.seed(1234)

  y.orig <- rnorm(125)

  locs.dim1 <- seq(1, 5, length = 5) + rnorm(5, 0, 0.1)
  locs.dim2 <- seq(1, 5, length = 5) + rnorm(5, 0, 0.1)
  locs.dim3 <- seq(1, 5, length = 5) + rnorm(5, 0, 0.1)
  locs <- as.matrix(expand.grid(locs.dim1, locs.dim2, locs.dim3))

  covmat.true <- Matern(rdist(locs), range = 1, smoothness = 0.5) + diag(125)

  y <- y.orig %*% chol(covmat.true)
  y <- as.vector(y)

  logparams <- create.initial.values.flex(0.9, c(2, 3), 0.5, 0.1, 1, 1)
  pseq <- create.param.sequence(1, 2)
  vec.approx <- vecchia_specify(locs, 5)
  LL.vecchia <- negloglik_vecchia_ST(
    logparams, y, vec.approx, pseq,
    c(1, 1, 2), 2
  )
  expect_equal(284.73, LL.vecchia, tolerance = 1e-2)

  logparams2 <- create.initial.values.flex(0.9, 1, 0.5, 0.1, 1, 1)
  pseq2 <- create.param.sequence(1)
  vec.approx2 <- vec.approx
  vec.approx2$locsord[, 1:2] <- vec.approx$locsord[, 1:2] / 2
  vec.approx2$locsord[, 3] <- vec.approx$locsord[, 3] / 3
  LL.vecchia2 <- negloglik_vecchia(logparams2, y, vec.approx2, pseq2)
  expect_equal(LL.vecchia, LL.vecchia2, tolerance = 1e-3)
})

test_that("negloglik.full", {
  set.seed(1234)

  y.orig <- rnorm(100)

  locs.dim <- seq(1, 10, length = sqrt(length(y.orig)))
  locs <- as.matrix(expand.grid(locs.dim, locs.dim))

  covmat.true <- Matern(rdist(locs), range = 1, smoothness = 0.5) + diag(100)

  y <- y.orig %*% chol(covmat.true)
  y <- as.vector(y)

  params.init <- rep(NA, 4)
  params.init[1] <- 0.9 * var(y)
  params.init[2] <- 1
  params.init[3] <- 0.5
  params.init[4] <- 0.1 * var(y)

  params.init <- create.initial.values.flex(
    params.init[1],
    params.init[2],
    params.init[3],
    params.init[4],
    1, 1
  )

  d <- rdist(locs)
  pseq <- create.param.sequence(1)

  res.optim.NM <- optim(
    par = params.init, fn = negloglik.full, d = d, y = y,
    param.seq = pseq, control = list(maxit = 5000)
  )

  LL.full <- negloglik.full(res.optim.NM$par, d, y, pseq)

  params.final <- c(
    exp(res.optim.NM$par[1:2]),
    gtools::inv.logit(res.optim.NM$par[3], 0, 2.5),
    exp(res.optim.NM$par[4])
  )

  pgp.params <- create.initial.values.flex(
    params.final[1],
    params.final[2],
    params.final[3],
    params.final[4],
    1, 1
  )

  LL.full.pgp <- mvnegloglik.full(pgp.params, list(locs), y, pseq)

  vec.approx <- vecchia_specify(locs, 99)

  LL.vecchia <- -1 * vecchia_likelihood(
    y, vec.approx, params.final[1:3],
    params.final[4]
  )

  vec.approx.pgp <- vecchia_Mspecify(list(locs), 99)
  vec.U.pgp <- createUMultivariate(vec.approx.pgp, c(params.final, 1))

  LL.pgp <- -1 * GPvecchia:::vecchia_likelihood_U(y, vec.U.pgp)

  expect_equal(173.315, LL.full, tolerance = 1e-3)
  # Univariate likelihood should equal the multivariate likelihood
  expect_equal(LL.full, LL.full.pgp, tolerance = 1e-3)
  # Full likelihood should equal both the univariate and multivariate
  # Vecchia approximations
  expect_equal(LL.full, LL.vecchia, tolerance = 1e-3)
  expect_equal(LL.full, LL.pgp, tolerance = 1e-3)
})

test_that("negloglik_full_ST", {
  set.seed(1234)

  y.orig <- rnorm(125)

  locs.dim1 <- seq(1, 5, length = 5) + rnorm(5, 0, 0.1)
  locs.dim2 <- seq(1, 5, length = 5) + rnorm(5, 0, 0.1)
  locs.dim3 <- seq(1, 5, length = 5) + rnorm(5, 0, 0.1)
  locs <- as.matrix(expand.grid(locs.dim1, locs.dim2, locs.dim3))

  covmat.true <- Matern(rdist(locs), range = 1, smoothness = 0.5) + diag(125)

  y <- y.orig %*% chol(covmat.true)
  y <- as.vector(y)

  logparams <- create.initial.values.flex(0.9, c(2, 3), 0.5, 0.1, 1, 1)
  pseq <- create.param.sequence(1, 2)
  LL.full <- negloglik_full_ST(
    logparams, locs, y, pseq,
    c(1, 1, 2), 2
  )
  expect_equal(286.84, LL.full, tolerance = 1e-2)

  logparams2 <- create.initial.values.flex(0.9, 1, 0.5, 0.1, 1, 1)
  pseq2 <- create.param.sequence(1)
  locs2 <- locs
  locs2[, 1:2] <- locs[, 1:2] / 2
  locs2[, 3] <- locs[, 3] / 3
  LL.full2 <- negloglik.full(logparams2, rdist(locs2), y, pseq2)
  expect_equal(LL.full, LL.full2, tolerance = 1e-3)
})

test_that("mvnegloglik", {
  load("sim_multivariate_big.RData")
  set.seed(1234)
  P <- 3
  Y <- cbind(runif(10), runif(10), runif(10))
  cor.matrix <- cor(Y)
  cov_mat <- c(cor.matrix[upper.tri(cor.matrix)])
  logparams <- create.initial.values.flex(
    rep(0.9, P), # marginal variance
    rep(0.5, P), # range
    rep(0.5, P), # smoothness
    rep(0.1, P), # nuggets
    cov_mat,
    P
  )
  pseq <- create.param.sequence(P)
  vec.approx <- vecchia_Mspecify(locs.list, 25)
  neg_likelihood <- mvnegloglik(
    logparams, vec.approx,
    unlist(y.list), pseq, P
  )
  expect_equal(78206.41, neg_likelihood, tolerance = 1e-2)
})

test_that("mvnegloglik_ST", {
  source("sim_multivariate_big_st.R")
  P <- 3
  Y <- cbind(runif(10), runif(10), runif(10))
  cor.matrix <- cor(Y)
  cov_mat <- c(cor.matrix[upper.tri(cor.matrix)])
  logparams <- create.initial.values.flex(
    rep(0.9, P), # marginal variance
    rep(c(2, 3), P), # range
    rep(0.5, P), # smoothness
    rep(0.1, P), # nuggets
    cov_mat,
    P
  )
  pseq <- create.param.sequence(P, 2)
  vec.approx <- vecchia_Mspecify(locs.list, 25)
  neg_likelihood <- mvnegloglik_ST(
    logparams, vec.approx,
    unlist(y.list), pseq, P, c(1, 1, 2), 2
  )
  expect_equal(35106.73, neg_likelihood, tolerance = 1e-2)

  vec.approx2 <- vec.approx
  for (i in 1:P) {
    vec.approx2$locsord[, 1:2] <- vec.approx$locsord[, 1:2] / 2
    vec.approx2$locsord[, 3] <- vec.approx$locsord[, 3] / 3
  }
  logparams2 <- logparams
  logparams2[pseq[2, 1]:pseq[2, 2]] <- 0
  neg_likelihood2 <- mvnegloglik_ST(
    logparams2, vec.approx2,
    unlist(y.list), pseq, P, c(1, 1, 2), 2
  )
  expect_equal(neg_likelihood, neg_likelihood2, tolerance = 1e-3)

  logparams <- create.initial.values.flex(
    rep(0.9, P), # marginal variance
    2:7, # range
    rep(0.5, P), # smoothness
    rep(0.1, P), # nuggets
    cov_mat,
    P
  )
  neg_likelihood <- mvnegloglik_ST(
    logparams, vec.approx,
    unlist(y.list), pseq, P, c(1, 1, 2), 2
  )
  expect_equal(36354.9, neg_likelihood, tolerance = 1e-2)

  vec.approx2 <- vec.approx
  vec.approx2$locsord[vec.approx$ondx == 1, 1:2] <-
    vec.approx$locsord[vec.approx$ondx == 1, 1:2] / 2
  vec.approx2$locsord[vec.approx$ondx == 1, 3] <-
    vec.approx$locsord[vec.approx$ondx == 1, 3] / 3
  vec.approx2$locsord[vec.approx$ondx == 2, 1:2] <-
    vec.approx$locsord[vec.approx$ondx == 2, 1:2] / 4
  vec.approx2$locsord[vec.approx$ondx == 2, 3] <-
    vec.approx$locsord[vec.approx$ondx == 2, 3] / 5
  vec.approx2$locsord[vec.approx$ondx == 3, 1:2] <-
    vec.approx$locsord[vec.approx$ondx == 3, 1:2] / 6
  vec.approx2$locsord[vec.approx$ondx == 3, 3] <-
    vec.approx$locsord[vec.approx$ondx == 3, 3] / 7

  neg_likelihood2 <- mvnegloglik_ST(
    logparams2, vec.approx2,
    unlist(y.list), pseq, P, c(1, 1, 2), 2
  )
  expect_equal(neg_likelihood, neg_likelihood2, tolerance = 1e-3)
})

test_that("mvnegloglik.full", {
  source("sim_multivariate_lik.R")

  pseq <- create.param.sequence(3)

  param.marg.var <- 0.9 * unlist(lapply(y.list, var))
  param.marg.scale <- rep(1, 3)
  param.marg.smooth <- rep(0.5, 3)
  param.marg.nugget <- 0.1 * unlist(lapply(y.list, var))
  param.rho <- rep(0.5, 3)
  params.init <- create.initial.values.flex(
    param.marg.var,
    param.marg.scale,
    param.marg.smooth,
    param.marg.nugget,
    param.rho, 3
  )

  params.init.flex.test <- c(
    log(c(param.marg.var, param.marg.scale)),
    gtools::logit(param.marg.smooth, 0, 2.5),
    log(param.marg.nugget), atanh(param.rho)
  )

  # Check create.initial.values.flex:
  expect_equal(params.init, params.init.flex.test, tolerance = 1e-3)

  cov.list.init <- create.cov.upper.flex(
    3,
    param.marg.var,
    param.marg.scale,
    param.marg.smooth,
    param.marg.nugget,
    param.rho)

  cov.mat.init <- cat.covariances(
    locs.list, cov.list.init$variance,
    cov.list.init$range,
    cov.list.init$smoothness,
    cov.list.init$nugget)

  # Check create.cov.upper.flex:
  cov.list.true <- readRDS("covlist.rds")
  expect_equal(cov.list.init$variance,
    cov.list.true$variance, tolerance = 1e-3)
  expect_equal(cov.list.init$range,
    cov.list.true$range, tolerance = 1e-3)
  expect_equal(cov.list.init$smoothness,
    cov.list.true$smoothness, tolerance = 1e-3)
  expect_equal(cov.list.init$nugget,
    cov.list.true$nugget, tolerance = 1e-3)
  # Check cat.covariances:
  cov.mat.true <- readRDS("covmat.rds")
  expect_equal(cov.mat.init, cov.mat.true, tolerance = 1e-3)

  res.optim.NM <- optim(
    par = params.init, fn = mvnegloglik.full,
    locs = locs.list, y = unlist(y.list),
    param.seq = pseq,
    method = "Nelder-Mead",
    control = list(trace = 0, maxit = 10000, reltol = 1e-4)
  )

  LL.full.mv <- mvnegloglik.full(
    res.optim.NM$par, locs.list,
    unlist(y.list), pseq
  )

  param.seq.begin <- pseq[, 1]
  param.seq.end <- pseq[, 2]
  params.init.final.t <- unlog.params(res.optim.NM$par, pseq, 3)

  cov.list <- create.cov.upper.flex(
    3,
    params.init.final.t[param.seq.begin[1]:param.seq.end[1]],
    params.init.final.t[param.seq.begin[2]:param.seq.end[2]],
    params.init.final.t[param.seq.begin[3]:param.seq.end[3]],
    params.init.final.t[param.seq.begin[4]:param.seq.end[4]],
    params.init.final.t[param.seq.begin[5]:param.seq.end[5]]
  )

  cov.mat <- cat.covariances(
    locs.list, cov.list$variance,
    cov.list$range,
    cov.list$smoothness, cov.list$nugget
  )

  LL.full.calc <- -1 * mvtnorm::dmvnorm(unlist(y.list),
    rep(0, length(unlist(y.list))),
    cov.mat,
    log = TRUE
  )

  vec.mapprox <- vecchia_Mspecify(locs.list, length(unlist(y.list)) - 1)
  U.mobj <- createUMultivariate(vec.mapprox, params.init.final.t)

  LL.vecchia.mv <- -1 * GPvecchia:::vecchia_likelihood_U(unlist(y.list), U.mobj)

  expect_equal(541.31, LL.full.mv, tolerance = 1e-3)
  expect_equal(LL.full.calc, LL.full.mv, tolerance = 1e-3)
  # Full likelihood should equal the Vecchia likelihood
  expect_equal(LL.full.mv, LL.vecchia.mv, tolerance = 1e-3)
})
