context("Log Likelihood Functions")

test_that("negloglik_vecchia_ST", {
  set.seed(7919)
  load("small_sim.Rdata")
  return(1) #this test isn't useful at the moment
  params = c(0.5, 0.5, 0.5, 0.5)
  locs <- locs_train
  X <- X_train
  Y <- Y_train
  locs.scaled <- locs/c(params[2], params[2], params[3])
  vecchia_approx = vecchia_specify(locs.scaled,25)
  beta.hat = rep(0,ncol(X))
  Y.hat <- as.matrix(X)%*%beta.hat
  res = as.double(Y-Y.hat)
  ord_locs <- locs[vecchia_approx$ord,]
  locs_mat <- cbind(ord_locs[,1], ord_locs[,2], ord_locs[,3])
  result <- negloglik_vecchia_ST(log(params),locs_mat,res,vecchia_approx)
  expect_equal(1597.9808688, result, tolerance=10e-5)
  vecchia.result<- optim(par=log(params),fn=negloglik_vecchia_ST,
                         locs=locs_mat, res=res, vecchia.approx=vecchia_approx, method = "Nelder-Mead",
                         control=list(trace=0))
  params_optim <- exp(vecchia.result$par)
  result_optim <- negloglik_vecchia_ST(log(params_optim),locs_mat,res,vecchia_approx)
  expect_equal(-6861.65048095, result_optim, tolerance=1e-5)
})

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

test_that("create.initial.values.flex", {
  set.seed(7919)
  P <- 1
  logparams <- create.initial.values.flex(rep(0.9,P), #marginal variance
                                          rep(0.5,P), #range
                                          rep(0.5,P), #smoothness
                                          rep(0.1,P), #nuggets
                                          NULL,
                                          P)
  expect_equal(c(-0.105, -0.693, -1.386, -2.303, 0), logparams, tolerance=1e-2)

  P <- 2
  Y <- cbind(runif(10),runif(10))
  cor.matrix <- cor(Y)
  cov_mat <- c(cor.matrix[upper.tri(cor.matrix)])
  logparams <- create.initial.values.flex(rep(0.9,P), #marginal variance
                                          rep(0.5,P), #range
                                          rep(0.5,P), #smoothness
                                          rep(0.1,P), #nuggets
                                          cov_mat,
                                          P)
  expect_equal(c(-0.105, -0.105, -0.693, -0.693, -1.386, -1.386, -2.303, -2.303, 0.584), logparams, tolerance=1e-2)
})

test_that("mvnegloglik", {
    source("sim_multivariate_big.R")
    P <- 3
    Y <- cbind(runif(10),runif(10), runif(10))
    cor.matrix <- cor(Y)
    cov_mat <- c(cor.matrix[upper.tri(cor.matrix)])
    logparams <- create.initial.values.flex(rep(0.9,P), #marginal variance
                                            rep(0.5,P), #range
                                            rep(0.5,P), #smoothness
                                            rep(0.1,P), #nuggets
                                            cov_mat,
                                            P)
    pseq = create.param.sequence(P)
    vec.approx <- vecchia_Mspecify(locs.list, 25)
    neg_likelihood <- mvnegloglik(logparams, vec.approx,
                                  unlist(y.list), pseq, P)
    expect_equal(34474.4, neg_likelihood, tolerance=1e-2)
})

#TODO implement this test
test_that("cat.covariances", {
  set.seed(7919)
  load("sim_spatial.Rdata")
  P <- 2
  Y <- cbind(runif(10),runif(10))
  cor.matrix <- cor(Y)
  cov_mat <- c(cor.matrix[upper.tri(cor.matrix)])
  param.sequence <- create.param.sequence(P)
  param.sequence.begin <- param.sequence[,1]
  param.sequence.end   <- param.sequence[,2]
  params <- create.initial.values.flex(rep(9.5,P), #marginal variance
                                          rep(15.0,P), #range
                                          rep(0.5,P), #smoothness
                                          rep(0.9,P), #nuggets
                                          cov_mat,
                                          P)
  params <- c(exp(params[1:param.sequence[4,2]]),
              tanh(params[param.sequence[5,1]:param.sequence[5,2]]))
  sig2 <- params[param.sequence.begin[1]:param.sequence.end[1]]
  range <- params[param.sequence.begin[2]:param.sequence.end[2]]
  smoothness <- params[param.sequence.begin[3]:param.sequence.end[3]]
  nugget <- params[param.sequence.begin[4]:param.sequence.end[4]]
  rho <- params[param.sequence.begin[5]:param.sequence.end[5]]

  cov.list  <- create.cov.upper.flex(P,sig2,range,smoothness,nugget,rho)
  cov.mat <- cat.covariances(list(locs_train, locs_train),cov.list$variance,cov.list$range,
                             cov.list$smoothness,cov.list$nugget)
  Y <- unlist(list(Y_train, Y_train))
  N <- length(Y)
  neg_likelihood <- -mvtnorm::dmvnorm(Y,rep(0,N),cov.mat,log=TRUE)
  expect_equal(5905, neg_likelihood, tolerance=10e-2)
})

