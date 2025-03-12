# Spatiotemporal Vecchia negative log likelihood.

negloglik_vecchia_ST <- function(logparms, res, vecchia.approx, param.seq, scaling, nscale) {
  parms <- unlog.params(logparms, param.seq, 1)
  locs.scaled <- vecchia.approx$locsord
  for (j in 1:nscale) {
    locs.scaled[, scaling == j] <- locs.scaled[, scaling == j] /
      parms[param.seq[2, 1] + j - 1]
  }
  vecchia.approx$locsord <- locs.scaled
  -vecchia_likelihood(
    res, vecchia.approx, c(
      parms[1], 1,
      parms[param.seq[3, 1]]
    ),
    parms[param.seq[4, 1]]
  )
}

# Spatial Vecchia negative log likelihood.

negloglik_vecchia <- function(logparms, res, vecchia.approx, param.seq) {
  parms <- unlog.params(logparms, param.seq, 1)
  -vecchia_likelihood(
    res, vecchia.approx, c(parms[1], parms[2], parms[3]),
    parms[4]
  )
}

# Spatiotemporal Full Kriging negative log likelihood.

negloglik_full_ST <- function(logparms, locs, y, param.seq, scaling, nscale) {
  parms <- unlog.params(logparms, param.seq, 1)
  locs.scaled <- locs
  for (j in 1:nscale) {
    locs.scaled[, scaling == j] <- locs.scaled[, scaling == j] /
      parms[param.seq[2, 1] + j - 1]
  }
  d <- fields::rdist(locs.scaled)
  N <- nrow(d)
  cov.mat <- parms[1] * fields::Matern(d,
    range = 1,
    smoothness = parms[param.seq[3, 1]]
  ) +
    parms[param.seq[4, 1]] * diag(N)
  -1 * mvtnorm::dmvnorm(y, rep(0, N), cov.mat, log = TRUE)
}

# Spatial Full Kriging negative log likelihood

negloglik.full <- function(logparams, d, y, param.seq) {
  params <- unlog.params(logparams, param.seq, 1)
  N <- nrow(d)
  cov.mat <- params[1] * fields::Matern(d,
    range = params[2],
    smoothness = params[3]
  ) +
    params[4] * diag(N)
  -1 * mvtnorm::dmvnorm(y, rep(0, N), cov.mat, log = TRUE)
}

#' Evaluation of the multivariate Vecchia likelihood
#'
#' This function is used to evaluate the multivariate Vecchia likelihood.
#'
#' @param z The observed data.
#' @param vecchia.approx A Vecchia object returned by
#' \code{\link{vecchia_Mspecify}}.
#' @param covparams Vector of covariance parameters. See
#' \code{\link{create_param_sequence}} or the examples below for details
#' about the format of this vector.
#'
#' @return The log likelihood implied by the multivariate Vecchia
#' approximation.
#'
#' @seealso \code{\link[GPvecchia]{vecchia_likelihood}},
#' \code{\link{vecchia_Mspecify}}, \code{\link{create_param_sequence}}
#'
#' @references
#' \itemize{
#' \item Katzfuss, M., and Guinness, J. "A general framework for Vecchia
#' approximations of Gaussian processes", Statistical Science (2021)
#' 36(1):124-141.
#' }
#'
#' @export
#' @examples
#' data(soil)
#' soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
#' locs <- as.matrix(soil[,1:2])
#' locsm <- list()
#' locsm[[1]] <- locsm[[2]] <- locs
#' soil.va <- vecchia_Mspecify(locsm, m=10)
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
#' vecchia_Mlikelihood(rnorm(nrow(locs)), soil.va, params)
vecchia_Mlikelihood <- function(z, vecchia.approx, covparams) {
  U.obj <- createUMultivariate(vecchia.approx, covparams)
  vecchia_likelihood_U <- getFromNamespace("vecchia_likelihood_U", "GPvecchia")
  res <- try(junk <- vecchia_likelihood_U(z, U.obj), silent = TRUE)
  if (inherits(res, "try-error")) {
    junk <- -1 * Inf
  }
  junk
}


##############################################################################
### Flexible Multivariate Matern Negative Loglikelihood Function ###########

mvnegloglik <- function(logparams, vecchia.approx, y, param.seq, P) {
  #  Input-
  #  logparams: A numeric vector of length (4*P)+(4*choose(P,2)).
  #             To construct these parameters we unlist a list of the 7 covariance
  #             categories- in order: (1) marginal variances, (2) Marginal ranges,
  #             (3) Marginal smoothness, (4) Nuggets, and
  #             (5) cross-covariance correlation. These seven parameters are to be
  #             created in a list. The variance, range, smoothness, and nugget
  #             have P terms, and the correlations have choose(P,2)
  #             terms. Use unlist() to create the vector of parameters.
  #
  # locs: list of the location coordinates, each outcome in y is a separate cell
  #       in the list
  # y    :  multivariate outcome, each out outcome in a separate entry in a list
  #  param.seq: The vector of parameter index sequences created by the function
  #            create_param_sequence - Used to identify the beginning and end
  #           index locations of each parameter.

  # P <- length(y)
  # transform the postively constrained parameters from log-space to normal-space
  params <- unlog.params(logparams, param.seq, P)
  -1 * vecchia_Mlikelihood(y, vecchia.approx, params)
}

##############################################################################
### Flexible Spatiotemporal Multivariate Matern Negative Loglikelihood Function ###########

mvnegloglik_ST <- function(logparams, vecchia.approx, y, param.seq, P, scaling, nscale) {
  #  Input-
  #  logparams: A numeric vector of length (4*P)+(4*choose(P,2)).
  #             To construct these parameters we unlist a list of the 7 covariance
  #             categories- in order: (1) marginal variances, (2) Marginal ranges,
  #             (3) Marginal smoothness, (4) Nuggets, and
  #             (5) cross-covariance correlation. These seven parameters are to be
  #             created in a list. The variance, range, smoothness, and nugget
  #             have P terms, and the correlations have choose(P,2)
  #             terms. Use unlist() to create the vector of parameters.
  #
  # locs: list of the location coordinates, each outcome in y is a separate cell
  #       in the list
  # y    :  multivariate outcome, each out outcome in a separate entry in a list
  #  param.seq: The vector of parameter index sequences created by the function
  #            create_param_sequence - Used to identify the beginning and end
  #           index locations of each parameter.

  # P <- length(y)
  # transform the postively constrained parameters from log-space to normal-space
  params <- unlog.params(logparams, param.seq, P)
  locs.scaled <- vecchia.approx$locsord
  for (i in 1:P) {
    for (j in 1:nscale) {
      locs.scaled[vecchia.approx$ondx == i, scaling == j] <-
        locs.scaled[vecchia.approx$ondx == i, scaling == j] /
          params[param.seq[2, 1] + nscale * (i - 1) + j - 1]
    }
  }
  vecchia.approx$locsord <- locs.scaled

  -1 * vecchia_Mlikelihood(y, vecchia.approx, c(
    params[1:param.seq[1, 2]],
    rep(1, param.seq[2, 2] - param.seq[2, 1] + 1),
    params[param.seq[3, 1]:param.seq[5, 2]]
  ))
}

##############################################################################
### Full Multivariate Matern Negative Loglikelihood Function ###########

mvnegloglik.full <- function(logparams, locs, y, param.seq) {
  #  Input-
  #  logparams: A numeric vector of length (4*P)+(4*choose(P,2)).
  #             To construct these parameters we unlist a list of the 7 covariance
  #             categories- in order: (1) marginal variances, (2) Marginal ranges,
  #             (3) Marginal smoothness, (4) Nuggets, and
  #             (5) cross-covariance correlation. These seven parameters are to be
  #             created in a list. The variance, range, smoothness, and nugget
  #             have P terms, and the correlations have choose(P,2)
  #             terms. Use unlist() to create the vector of parameters.
  #
  # locs: list of the location coordinates, each outcome in y is a separate cell
  #       in the list
  # y    :  multivariate outcome, each out outcome in a separate entry in a list
  #  param.seq: The vector of parameter index sequences created by the function
  #            create_param_sequence - Used to identify the beginning and end
  #           index locations of each parameter.

  # P <- length(y)
  # transform the positively constrained parameters from log-space to normal-space
  P <- length(locs)
  params <- unlog.params(logparams, param.seq, P)
  sig2 <- params[param.seq[1, 1]:param.seq[1, 2]]
  range <- params[param.seq[2, 1]:param.seq[2, 2]]
  smoothness <- params[param.seq[3, 1]:param.seq[3, 2]]
  nugget <- params[param.seq[4, 1]:param.seq[4, 2]]
  rho <- NA
  if (P > 1) {
    rho <- params[param.seq[5, 1]:param.seq[5, 2]] # Rho is estimated and constrained between -1 and 1
  }

  # set up the covariance parameters for multivariate input

  # This function takes the marginal and cross-covariance parameters and
  # aligns them in lists
  cov.list <- create.cov.upper.flex(P, sig2, range, smoothness, nugget, rho)

  cov.mat <- cat.covariances(
    locs, cov.list$variance, cov.list$range,
    cov.list$smoothness, cov.list$nugget
  )

  # print(sum(cov.mat))
  # Y <- unlist(y)
  # Y <- matrix(unlist(y), ncol=P)
  N <- length(y)
  # N <-nrow(Y)
  -mvtnorm::dmvnorm(y, rep(0, N), cov.mat, log = TRUE) # + lambda*norm(unlist(cov.list),type = "F")
}

##############################################################################
create.cov.upper.flex <- function(P, marg.var, marg.range, marg.smooth, nugget, R.corr) {
  # Create the symmetrical marginal+cross-covariance flexible Matern from the
  # given parameters. Output is a list of the 4 Matern parameters as matrices
  sig2.mat <- diag(marg.var, P, P)
  range.mat <- diag(marg.range, P, P)
  smoothness.mat <- diag(marg.smooth, P, P)
  nugget.mat <- diag(nugget, P, P)
  if (P > 1) {
    combs <- gtools::combinations(P, 2)
    for (iter in seq_len(nrow(combs))) {
      i <- combs[iter, 1]
      j <- combs[iter, 2]

      smoothness.mat[i, j] <- (marg.smooth[i] + marg.smooth[j]) / 2
      range.mat[i, j] <- 1 / sqrt(((1 / marg.range[i])^2 + (1 / marg.range[j])^2) / 2)

      s1 <- sqrt(marg.var[i] * marg.var[j])
      s2 <- ((1 / marg.range[i])^marg.smooth[i] * (1 / marg.range[j])^marg.smooth[j]) /
        ((1 / range.mat[i, j])^(2 * smoothness.mat[i, j]))
      s3 <- gamma(smoothness.mat[i, j]) / (sqrt(gamma(marg.smooth[i])) * sqrt(gamma(marg.smooth[j])))
      s4 <- R.corr[iter]
      sig2.mat[i, j] <- s1 * s2 * s3 * s4
    }
  }

  list(
    "variance" = sig2.mat,
    "range" = range.mat,
    "smoothness" = smoothness.mat,
    "nugget" = nugget.mat
  )
}


##############################################################################
### Calculate the Matern marginal and cross-covariance super-matrix #########

cat.covariances <- function(locs.list, sig2, range, smoothness, nugget) {
  # cat.covariance: This functions takes the locations and flexible Matern parameters
  #                 and pieces together the marginal+cross-covariance matrix
  #
  #   Inputs: locs- P dimension list of the locations,
  #           sig2 : marginal variances as matrix
  #           range : marginal ranges as matrix
  #           smoothness: marginal smoothness as matrix
  #           nugget as matrix

  dims <- sapply(locs.list, nrow)
  total.dims <- sum(dims)
  cov.mat.out <- matrix(NA, nrow = total.dims, ncol = total.dims)

  l <- length(locs.list)
  combs <- gtools::combinations(l, 2, repeats.allowed = TRUE)
  for (iter in seq_len(nrow(combs))) {
    i <- combs[iter, 1]
    j <- combs[iter, 2]
    # d <- fields::rdist.earth(locs.list[[i]],locs.list[[j]],miles = FALSE)
    d <- fields::rdist(locs.list[[i]], locs.list[[j]])
    # Calculate the covariance matrix - if/then based on its location in the super-matrix
    N <- nrow(d)
    if (i == j) { # To accommodate varying size outcomes- the nugget is not included on cross-covariances
      cov.mat.ij <- sig2[i, j] * fields::Matern(d,
        range = range[i, j], smoothness =
          smoothness[i, j]
      ) +
        nugget[i, j] * diag(N)
    } else {
      cov.mat.ij <- sig2[i, j] * fields::Matern(d,
        range = range[i, j], smoothness =
          smoothness[i, j]
      )
    }


    if (combs[iter, 1] == 1) {
      row.idx <- 1:dims[1]
    } else {
      start.row <- cumsum(dims)[combs[iter, 1] - 1] + 1
      row.idx <- start.row:(start.row + (dims[combs[iter, 1]] - 1))
    }

    if (combs[iter, 2] == 1) {
      col.idx <- 1:dims[1]
    } else {
      start.col <- cumsum(dims)[combs[iter, 2] - 1] + 1
      col.idx <- start.col:(start.col + (dims[combs[iter, 2]] - 1))
    }

    cov.mat.out[row.idx, col.idx] <- cov.mat.ij
    if (!identical(row.idx, col.idx)) {
      cov.mat.out[col.idx, row.idx] <- t(cov.mat.ij)
    }
  }


  cov.mat.out
}

##############################################################################
### Create the likelihood initial values                    #########

create.initial.values.flex <- function(marg.var, marg.range, marg.smooth, nugget, R.corr, P) {
  # Log-transform the covariance parameters and arrange in the proper order
  # for the likelihood function
  logparams.init <- c(
    log(marg.var), log(marg.range),
    gtools::logit(marg.smooth, 0, 2.5), log(nugget)
  )
  if (P > 1) {
    logparams.init <- c(logparams.init, atanh(R.corr))
  } else {
    logparams.init <- c(logparams.init, 0)
  }
  logparams.init
}

##############################################################################
### Transform the log Matern parameters back to the original #########

unlog.params <- function(logparams, param.seq, P) {
  params <- c(
    exp(logparams[1:param.seq[2, 2]]),
    gtools::inv.logit(logparams[param.seq[3, 1]:param.seq[3, 2]], 0, 2.5),
    exp(logparams[param.seq[4, 1]:param.seq[4, 2]])
  )
  if (P > 1) {
    params <- c(params, tanh(logparams[param.seq[5, 1]:param.seq[5, 2]]))
  } else {
    params <- c(params, 1)
  }
  params
}
