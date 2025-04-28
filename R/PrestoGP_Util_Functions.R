################################################################################
# SCAD/MCP penalty values
################################################################################
SCAD_penalty <- function(x, lambda, gamma) {
  idx1 <- abs(x) <= lambda
  idx2 <- abs(x) > lambda & abs(x) <= gamma * lambda
  idx3 <- abs(x) > gamma * lambda
  pen1 <- lambda * abs(x[idx1])
  pen2 <- (2 * gamma * lambda * abs(x[idx2]) - x[idx2]^2 - lambda^2) /
    (2 * (gamma - 1))
  sum(pen1) + sum(pen2) + sum(idx3) * lambda^2 * (gamma + 1) / 2
}

MCP_penalty <- function(x, lambda, gamma) {
  idx <- abs(x) <= (gamma * lambda)
  sum(lambda * abs(x[idx]) - x[idx]^2 / (2 * gamma)) +
    sum(!idx) * gamma * lambda^2
}

#' glmnet_penalty
#'
#' Compute a combination of L1 and L2 penalties
#'
#' @param beta the beta coefficients
#' @param lambda the penalty factor (how big the total penalty is)
#' @param alpha the L1 vs L2 coefficient (0 means all L2, 1 means all L1)
#'
#' @return a penalty score
#' @noRd
#glmnet_penalty <- function(beta, lambda, alpha) {
#  return(lambda * ((1 - alpha) * sqrt(sum(beta^2)) + alpha * sum(abs(beta))))
#}

# Ordinary Spatiotemporal Kriging Prediction with a Local S-T neighborhood.
#
# @param new_coords Locations to predict
# @param obs_coords Training locations
# @param Y_obs Training dependent variable vector
# @param cov.pars Estimated covariance parameters
# @param NN number of nearest neighbors to use
#
# @return A dataframe with predicted means and variances
#Kr_pred <- function(new_coords, obs_coords, Y_obs, cov.pars, NN) {
#  Kr.prediction <- matrix(NA, nrow = nrow(new_coords), ncol = 1)
#  Kr.Var <- matrix(NA, nrow = nrow(new_coords), ncol = 1)
#
#  new_coords_scaled <- cbind(new_coords[, 1] / cov.pars[2], new_coords[, 2] / cov.pars[2], new_coords[, 3] / cov.pars[3])
#  obs_coords_scaled <- cbind(obs_coords[, 1] / cov.pars[2], obs_coords[, 2] / cov.pars[2], obs_coords[, 3] / cov.pars[3])
#  df_new <- new_coords_scaled
#  df_obs <- obs_coords_scaled

#  for (i in seq_len(nrow(new_coords))) {
# print(i)
#    locs.test.i <- rep.row(df_new[i, ], nrow(df_obs))
#    dist.op <- fields::rdist.vec(locs.test.i, df_obs)
#    Sigma.op <- cov.pars[1] * fields::Exponential(dist.op, range = 1)
#    s.op <- sort(dist.op, index.return = TRUE)
#    s.idx <- s.op$ix[1:NN]
# Get the closest "NN" observations
#    dist.oo <- fields::rdist(df_obs[s.idx, ], df_obs[s.idx, ])
#    Sigma.oo <- cov.pars[4] * diag(NN) +
#      cov.pars[1] * fields::Exponential(dist.oo, range = 1)
#    oo <- Y_obs[s.idx]
## kriging predictor (posterior mean)
#    mean_trend <- mean(oo)
#    Kr.prediction[i] <- mean_trend + t(Sigma.op[s.idx]) %*% solve(Sigma.oo, oo - mean_trend)
#
#    Kr.Var[i] <- cov.pars[1] - t(Sigma.op[s.idx]) %*% solve(Sigma.oo, Sigma.op[s.idx])
#  }
#
#  Kr.Data <- data.frame(Kr.prediction, Kr.Var)
#
#  return(Kr.Data)
#}

# ST_Krig_Param_Avg
#
# Spatiotemporal Kriging Maximum Likelihood Estimtion based on an average of multiple random subsets.
#
# @param Y A vector containing training values for the dependent variable.
# @param locs A matrix containing the training spatial coordinates and times.
# @param p The number of data points to sample with each iteration.
# @param k The number of optimization iterations.
#
# @return a vector containing the covariance parameters
#ST_Krig_Param_Avg <- function(Y, locs, p, k = 10) {
#  n <- length(Y)
#  mdl.geo.fit.avg <- list()
#  mdl.geo.fit.avg <- matrix(NA, nrow = k, ncol = 4)
#
#  for (i in 1:k) {
#    set.seed(i)
#    OK.fold <- sample(1:n, p, replace = FALSE)
#    Y.train <- Y[OK.fold]
#    locs.train <- locs[OK.fold, ]
#    D.sample <- rdist(locs.train[, 1:2])
#    t.sample <- rdist(locs.train[, 3])
#    theta.hat <- c(.9 * var(Y.train), mean(D.sample) / 4, mean(t.sample) / 4, 0.1 * var(Y.train)) # var,s-range,t-range,nugget
#    full.result <- optim(
#      par = log(theta.hat), fn = negloglik_full_ST,
#      locs = locs.train, y = Y.train, N = p,
#      control = list(trace = FALSE, maxit = 400)
#    )
#    mdl.geo.fit.avg[i, ] <- full.result$par
#  }
#
#  covparam <- exp(apply(mdl.geo.fit.avg, 2, mean))
#
#  return(covparam)
#}


################################################################################
## replicate rows, helps with vector distance
################################################################################
#rep.row <- function(x, n) {
#  matrix(rep(x, each = n), nrow = n)
#}


################################################################################
## Transform Spatiotemporal data to i.i.d.
################################################################################
transform_iid <- function(data, vecchia.approx, covparms, nuggets) {
  # compute required matrices
  U.obj <- createU(vecchia.approx, covparms, nuggets)
  V.ord <- U2V(U.obj)
  U.z <- U.obj$U[!U.obj$latent, ]
  U.y <- U.obj$U[U.obj$latent, ]

  # compute transformed data in parts
  part1.ord <- crossprod(U.z, data[U.obj$ord.z, ]) # C.hat^-1
  temp1 <- U.y %*% part1.ord
  revord <- rev(seq_len(nrow(temp1)))
  temp2 <- solve(V.ord, temp1[revord, ])
  part2.rev <- solve(Matrix::t(V.ord), temp2)
  part2.ord <- crossprod(U.y, part2.rev[revord, ])
  transform.ord <- part1.ord - part2.ord

  # return to original ordering
  # orig.order <- order(U.obj$ord)
  # transformed.data <- transform.ord[orig.order, ]
  transform.ord
}

################################################################################
## Transform Multivariate Spatiotemporal data to i.i.d.
################################################################################
transform_miid <- function(data, vecchia.approx, params) {
  # compute required matrices
  U.obj <- createUMultivariate(vecchia.approx, params)
  V.ord <- U2V(U.obj)
  U.z <- U.obj$U[!U.obj$latent, ]
  U.y <- U.obj$U[U.obj$latent, ]

  # compute transformed data in parts
  part1.ord <- crossprod(U.z, data[U.obj$ord.z, ]) # C.hat^-1
  temp1 <- U.y %*% part1.ord #
  revord <- rev(seq_len(nrow(temp1)))
  temp2 <- solve(V.ord, temp1[revord, ])
  part2.rev <- solve(Matrix::t(V.ord), temp2)
  part2.ord <- crossprod(U.y, part2.rev[revord, ])
  transform.ord <- part1.ord - part2.ord

  # return to original ordering
  # orig.order <- order(U.obj$ord)
  # transformed.data <- transform.ord[orig.order, ]
  transform.ord
}

######  GPvecchia local function
###### compute V for posterior inference - needed for transform.iid   #######
U2V <- function(U.obj) {
  U.y <- U.obj$U[U.obj$latent, ]

  if (U.obj$cond.yz == "zy") {
    V.ord <- revMat(U.y[, U.obj$latent, drop = FALSE])
  } else if (U.obj$ord.pred != "obspred") {
    W <- Matrix::tcrossprod(U.y)
    W.rev <- revMat(W)
    res <- try(V.ord <- Matrix::t(Matrix::chol(W.rev)), silent = TRUE)
    if (inherits(res, "try-error")) {
      W.rev <- Matrix::nearPD(W.rev)
      V.ord <- Matrix::t(Matrix::chol(W.rev))
    }
  } else { # for obspred ordering

    last.obs <- max(which(!U.obj$latent))
    latents.before <- sum(U.obj$latent[1:last.obs])
    latents.after <- sum(U.obj$latent[-(1:last.obs)])

    # pred columns are unchanged
    V.pr <- revMat(U.y[, (last.obs + 1):ncol(U.y), drop = FALSE])

    # have to compute cholesky for obs block
    U.oo <- U.y[1:latents.before, 1:last.obs]
    A <- Matrix::tcrossprod(U.oo)
    A.rev <- revMat(A)
    res <- try(V.oor <- Matrix::t(Matrix::chol(A.rev)), silent = TRUE)
    if (inherits(res, "try-error")) {
      A.rev <- Matrix::nearPD(A.rev)
      V.oor <- Matrix::t(Matrix::chol(A.rev))
    }

    # combine the blocks into one matrix
    zeromat.sparse <- Matrix::sparseMatrix(c(), c(), dims = c(latents.after, latents.before))
    V.or <- rbind(zeromat.sparse, V.oor)

    V.ord <- methods::as(cbind(V.pr, V.or), "triangularMatrix")
  }

  V.ord
}

###########################################################
## Reverse order of matrix rows,cols
revMat <- function(mat) {
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    mat.out <- mat
  } else {
    row_seq <- rev(seq_len(nrow(mat)))
    col_seq <- rev(seq_len(ncol(mat)))
    mat.out <- mat[row_seq, col_seq, drop = FALSE]
  }
  mat.out
}

#' Multivariate Vecchia prediction
#'
#' This function is used to make predictions based on multivariate Vecchia
#' models. It is a multivariate version of
#' \code{\link[GPvecchia]{vecchia_prediction}}.
#'
#' @param z The observed data.
#' @param vecchia.approx A Vecchia object returned by
#' \code{\link{vecchia_Mspecify}}.
#' @param covparms Vector of covariance parameters. See
#' \code{\link{create_param_sequence}} or the examples below for details
#' about the format of this vector.
#' @param var.exact Should prediction variances by computed exactly, or is a
#' (faster) approximation acceptable? See
#' \code{\link[GPvecchia]{vecchia_prediction}}.
#' @param return.values Values that should be returned. Possible values
#' include "mean", "meanvar", "meanmat", and "all". See
#' \code{\link[GPvecchia]{vecchia_prediction}}. Defaults to "mean".
#'
#' @return The posterior means/variances/V matrices at the observed and
#' unobserved locations. See \code{\link[GPvecchia]{vecchia_prediction}}.
#'
#' @seealso \code{\link[GPvecchia]{vecchia_prediction}},
#' \code{\link{vecchia_Mspecify}}, \code{\link{create_param_sequence}}
#'
#' @references
#' \itemize{
#' \item Katzfuss, M., and Guinness, J. "A general framework for Vecchia
#' approximations of Gaussian processes", Statistical Science (2021)
#' 36(1):124-141.
#' \item Katzfuss, M., Guinness, J., Gong, W. and Zilber, D. "Vecchia
#' approximations of Gaussian-process predictions", Journal of Agricultural,
#' Biological and Environmental Statistics (2020) 25:383-414.
#' }
#'
#' @export
#' @examples
#' data(soil)
#' soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
#' locs <- as.matrix(soil[,1:2])
#' locsm <- list()
#' locsm[[1]] <- locsm[[2]] <- locs
#' locsp <- locsm
#' locsp[[1]] <- locsp[[1]] + 0.5
#' locsp[[2]] <- locsp[[2]] - 0.5
#' soil.vap <- vecchia_Mspecify(locsm, m=10, locs.list.pred=locsp)
#'
#' pseq <- create_param_sequence(2)
#' # Initialize the vector of covariance parameters
#' params <- rep(NA, pseq[5,2])
#' # Sigma parameters:
#' params[pseq[1,1]:pseq[1,2]] <- c(100, 80)
#' # Scale parameters:
#' params[pseq[2,1]:pseq[2,2]] <- c(60, 50)
#' # Smoothness parameters:
#' params[pseq[3,1]:pseq[3,2]] <- c(0.5, 0.5)
#' # Nuggets:
#' params[pseq[4,1]:pseq[4,2]] <- c(30, 30)
#' # Correlation:
#' params[pseq[5,1]:pseq[5,2]] <- -0.9
#'
#' soil.yhat <- vecchia_Mprediction(rnorm(nrow(locs)), soil.vap, params)
vecchia_Mprediction <- function(z, vecchia.approx, covparms, var.exact = NULL, return.values = "mean") {
  removeNAs <- getFromNamespace("removeNAs", "GPvecchia")
  removeNAs()
  U.obj <- createUMultivariate(vecchia.approx, covparms)
  V.ord <- U2V(U.obj)
  #    if (length(U.obj$zero.nugg) > 0)
  #        warning("Rows/cols of V have been removed for data with zero noise")
  vecchia_mean <- getFromNamespace("vecchia_mean", "GPvecchia")
  vecchia.mean <- vecchia_mean(z, U.obj, V.ord)
  return.list <- list(
    mu.pred = vecchia.mean$mu.pred, mu.obs = vecchia.mean$mu.obs,
    var.pred = NULL, var.obs = NULL, V.ord = NULL, U.obj = NULL
  )
  if (return.values == "meanmat" || return.values == "all") {
    return.list$V.ord <- V.ord
    return.list$U.obj <- U.obj
  }
  if (return.values == "meanvar" || return.values == "all") {
    if (is.null(var.exact)) {
      var.exact <- (sum(!vecchia.approx$obs) < 4 * 10000)
    }
    vecchia_var <- getFromNamespace("vecchia_var", "GPvecchia")
    vars.vecchia <- vecchia_var(U.obj, V.ord, exact = var.exact)
    return.list$var.pred <- vars.vecchia$vars.pred
    return.list$var.obs <- vars.vecchia$vars.obs
  }
  return.list
}

# noise_locs
#
# Adds a small amount of noise to the locs to avoid numeric issues related
# to duplicated locs
noise_locs <- function(locs, eps = 1e-4) {
  ee <- min(apply(locs, 2, stats::sd))
  n <- nrow(locs)
  p <- ncol(locs)
  locs <- locs + matrix(
    ee * eps *
      stats::rnorm(n * p),
    n, p
  )
  locs
}

# eliminate_dupes
#
# Eliminates duplicate locs by adding noise if needed
eliminate_dupes <- function(locs, locs.pred = NULL) {
  locs.all <- NULL
  for (i in seq_along(locs)) {
    locs.all <- rbind(locs.all, locs[[i]])
  }
  if (!is.null(locs.pred)) {
    for (i in seq_along(locs.pred)) {
      locs.all <- rbind(locs.all, locs.pred[[i]])
    }
  }
  if (sum(duplicated(locs.all)) > 0) {
    for (i in seq_along(locs)) {
      locs[[i]] <- noise_locs(locs[[i]])
    }
    if (!is.null(locs.pred)) {
      for (i in seq_along(locs.pred)) {
        locs.pred[[i]] <- noise_locs(locs.pred[[i]])
      }
    }
  }
  list(locs = locs, locs.pred = locs.pred)
}

lod_reg_mi <- function(y, X, lodu, lodl, miss, n.mi = 10, eps = 0.01,
  maxit = 10, penalty, alpha, parallel, cluster, foldid, verbose) {
  lodu <- lodu[miss]
  lodl <- lodl[miss]
  last.coef <- rep(Inf, ncol(X) + 1)
  if (penalty == "lasso" || penalty == "relaxed") {
    relax <- penalty == "relaxed"
    cur.glmnet <- cv.glmnet(X, y, parallel = parallel, foldid = foldid,
      alpha = alpha, relax = relax)
    cur.coef <- as.matrix(predict(cur.glmnet, type = "coefficients",
        s = "lambda.min", gamma = "gamma.min"))
  } else {
    cur.ncvreg <- cv.ncvreg.wrap(X, y, cluster = cluster, foldid = foldid,
      penalty = penalty, alpha = alpha)
    cur.coef <- as.matrix(coef(cur.ncvreg))
  }
  itn <- 0
  while (max(abs(cur.coef - last.coef)) > eps && itn <= maxit) {
    itn <- itn + 1
    last.coef <- cur.coef
    miss.means <- X[miss, ] %*% last.coef[-1] + last.coef[1]
    obs.means <- X[!miss, ] %*% last.coef[-1] + last.coef[1]
    cur.resid <- obs.means - y[!miss]
    cur.sd <- sqrt(sum(cur.resid^2) / (sum(!miss) - length(last.coef)))
    coef.mat <- matrix(nrow = n.mi, ncol = length(cur.coef))
    for (i in 1:n.mi) {
      y[miss] <- rtruncnorm(sum(miss), a = lodl, b = lodu, mean = miss.means,
        sd = cur.sd)
      if (penalty == "lasso" || penalty == "relaxed") {
        cur.glmnet <- cv.glmnet(X, y, parallel = parallel, foldid = foldid,
          alpha = alpha, relax = relax)
        coef.mat[i, ] <- as.matrix(predict(cur.glmnet, type = "coefficients",
            s = "lambda.min", gamma = "gamma.min"))
      } else {
        cur.ncvreg <- cv.ncvreg.wrap(X, y, cluster = cluster, foldid = foldid,
          penalty = penalty, alpha = alpha)
        coef.mat[i, ] <- as.matrix(coef(cur.ncvreg))
      }
    }
    cur.coef <- colMeans(coef.mat)
    if (verbose) {
      cat("LOD imputation iteration", itn, "complete", "\n")
      cat("Current coefficients:", "\n")
      print(cur.coef)
    }
  }
  miss.means <- X[miss, ] %*% cur.coef[-1] + cur.coef[1]
  obs.means <- X[!miss, ] %*% cur.coef[-1] + cur.coef[1]
  cur.resid <- obs.means - y[!miss]
  cur.sd <- sqrt(sum(cur.resid^2) / (sum(!miss) - length(last.coef)))
  y.impute <- etruncnorm(a = lodl, b = lodu, mean = miss.means, sd = cur.sd)
  list(coef = cur.coef, y.impute = y.impute)
}

MMatern_cov <- function(locs, y_ndx, covparams, P) {
  param.seq <- create_param_sequence(P)
  sigma <- covparams[param.seq[1, 1]:param.seq[1, 2]]
  rangep <- covparams[param.seq[2, 1]:param.seq[2, 2]]
  smoothness <- covparams[param.seq[3, 1]:param.seq[3, 2]]
  nuggets <- covparams[param.seq[4, 1]:param.seq[4, 2]]

  rho <- covparams[param.seq[5, 1]:param.seq[5, 2]]
  rho.mat <- matrix(0, nrow = P, ncol = P)
  rho.mat[upper.tri(rho.mat, diag = FALSE)] <- rho
  rho.mat <- rho.mat + t(rho.mat)
  diag(rho.mat) <- 1

  Sigma.hat <- matrix(nrow = nrow(locs), ncol = nrow(locs))
  for (i in seq_len(P)) {
    for (j in i:P) {
      which.i <- which(y_ndx == i)
      which.j <- which(y_ndx == j)
      if (sum(which.i) > 0 && sum(which.j) > 0) {
        smooth.ii <- smoothness[i]
        smooth.jj <- smoothness[j]
        smooth.ij <- (smooth.ii + smooth.jj) / 2
        alpha.ii <- 1 / rangep[i]
        alpha.jj <- 1 / rangep[j]
        alpha.ij <- sqrt((alpha.ii^2 + alpha.jj^2) / 2)
        Sigma.hat[which.i, which.j] <- rho.mat[i, j] * sqrt(sigma[i]) *
          sqrt(sigma[j]) * alpha.ii^smooth.ii * alpha.jj^smooth.jj *
          gamma(smooth.ij) / (alpha.ij^(2 * smooth.ij) *
            sqrt(gamma(smooth.ii) * gamma(smooth.jj))) *
          Matern(rdist(locs[which.i, drop = FALSE],
              locs[which.j, drop = FALSE]),
            alpha = alpha.ij, smoothness = smooth.ij)
        if (i != j) {
          Sigma.hat[which.j, which.i] <- t(Sigma.hat[which.i, which.j])
        } else {
          Sigma.hat[which.i, which.i] <- Sigma.hat[which.i, which.i] +
            nuggets[i] * diag(nrow = length(which.i))
        }
      }
    }
  }
  Sigma.hat
}

rtmvn_snn2 <- function(y, cens_lb, cens_ub, mask_cens, NN, cov_array) {
  ind_cens <- which(mask_cens)
  k <- 1
  for (i in ind_cens) {
    NN_row <- NN[i, ]
    covmat_sub <- cov_array[, , k]
    k <- k + 1
    mask_cens_sub <- mask_cens[NN_row]
    y_sub <- y[NN_row]
    cens_lb_sub <- cens_lb[NN_row]
    cens_ub_sub <- cens_ub[NN_row]
    n_cens_sub <- sum(mask_cens_sub)
    if (n_cens_sub == length(NN_row)) {
      cond_covmat_sub_cens <- covmat_sub
      cond_mean_sub_cens <- rep(0, n_cens_sub)
    } else {
      tmp_mat <- covmat_sub[mask_cens_sub, !mask_cens_sub] %*%
        solve(covmat_sub[!mask_cens_sub, !mask_cens_sub])
      cond_covmat_sub_cens <- covmat_sub[mask_cens_sub, mask_cens_sub] -
        tmp_mat %*% covmat_sub[!mask_cens_sub, mask_cens_sub]
      cond_covmat_sub_cens[lower.tri(cond_covmat_sub_cens)] <-
        t(cond_covmat_sub_cens)[lower.tri(cond_covmat_sub_cens)]
      cond_mean_sub_cens <- as.vector(tmp_mat %*% y_sub[!mask_cens_sub])
    }
    samp_cens_sub <- t(TruncatedNormal::rtmvnorm(
      1, cond_mean_sub_cens,
      cond_covmat_sub_cens,
      cens_lb_sub[mask_cens_sub], cens_ub_sub[mask_cens_sub]
    ))
    y[i] <- samp_cens_sub[1]
    mask_cens[i] <- FALSE
  }
  y
}

cv.ncvreg.wrap <- function(X, y, cluster, foldid, penalty, ...) {
  if (!is.null(cluster)) {
    if (!is.null(foldid)) {
      junk <- cv.ncvreg(X, y, cluster = cluster, fold = foldid,
        penalty = penalty, ...)
    } else {
      junk <- cv.ncvreg(X, y, cluster = cluster, penalty = penalty, ...)
    }
  } else {
    if (!is.null(foldid)) {
      junk <- cv.ncvreg(X, y, fold = foldid, penalty = penalty, ...)
    } else {
      junk <- cv.ncvreg(X, y, penalty = penalty, ...)
    }
  }
  junk
}
