#' Multivariate Vecchia PrestoGP model class
#'
#' This class is used to create multivariate models with a likelihood function
#' conditioned on a subset of the observations (i.e., Vecchia models). See
#' \code{\link{PrestoGPModel}} for a description of the slots.
#'
#' @seealso \code{\link{PrestoGPModel}}
#'
#' @export MultivariateVecchiaModel
#'
#' @examples
#' pgp.mmodel <- new("MultivariateVecchiaModel", n_neighbors = 25)
#' @include PrestoGP_Model.R
MultivariateVecchiaModel <- setClass("MultivariateVecchiaModel",
  contains = "PrestoGPModel"
)

validityMultivariateVecchiaModel <- function(object) {
  TRUE
}
setValidity(
  "MultivariateVecchiaModel",
  validityMultivariateVecchiaModel
)

setMethod("initialize", "MultivariateVecchiaModel", function(.Object, n_neighbors = 25, ...) {
  .Object@n_neighbors <- n_neighbors
  .Object@min_m <- 3
  .Object <- callNextMethod()
  validObject(.Object)
  .Object
})

#' @rdname get_Y
setMethod("get_Y", "MultivariateVecchiaModel",
  function(model) {
    Y.all <- model@Y_train
    Y.out <- list()
    for (i in seq_along(model@locs_train)) {
      cur.y <- seq_len(nrow(model@locs_train[[i]]))
      Y.out[[i]] <- Y.all[cur.y] + model@Y_bar[i]
      Y.all <- Y.all[-cur.y]
    }
    Y.out
  }
)

#' @rdname prestogp_predict
setMethod("prestogp_predict", "MultivariateVecchiaModel",
  function(model, X, locs, m = NULL, ordering.pred = c("obspred", "general"), pred.cond = c("independent", "general"), return.values = c("mean", "meanvar")) {
    # validate parameters
    ordering.pred <- match.arg(ordering.pred)
    pred.cond <- match.arg(pred.cond)
    return.values <- match.arg(return.values)
    pred.list <- check_input_pred(model, X, locs)
    X <- pred.list$X
    Y_bar <- pred.list$Y_bar
    if (is.null(m)) { # m defaults to the value used for training
      m <- model@n_neighbors
    }
    if (m < model@min_m) {
      stop(paste("m must be at least ", model@min_m, sep = ""))
    }
    if (m >= nrow(model@X_train)) {
      warning("Conditioning set size m chosen to be >=n. Changing to m=n-1")
      m <- nrow(model@X_train) - 1
    }

    if (model@penalty == "lasso" || model@penalty == "relaxed") {
      # Vecchia prediction at new locations
      Vecchia.Pred <- predict(model@linear_model, newx = X, s = "lambda.min",
        gamma = "gamma.min")
      # Vecchia trend prediction at observed data
      Vecchia.hat <- predict(model@linear_model, newx = model@X_train,
        s = "lambda.min", gamma = "gamma.min")
    } else {
      # Vecchia prediction at new locations
      Vecchia.Pred <- predict(model@linear_model, X)
      # Vecchia trend prediction at observed data
      Vecchia.hat <- predict(model@linear_model, model@X_train)
    }

    # Test set prediction
    res <- model@Y_train - Vecchia.hat

    locs.train.scaled <- scale_locs(model, model@locs_train)
    locs.scaled <- scale_locs(model, locs)
    if (model@common_scale & (length(model@locs_train) > 1)) {
      locs.nd <- eliminate_dupes(locs.train.scaled, locs.scaled)
      locs.train.scaled <- locs.nd$locs
      locs.scaled <- locs.nd$locs.pred
    }
    vec.approx.test <- vecchia_Mspecify(locs.train.scaled, m,
      locs.list.pred = locs.scaled,
      ordering.pred = ordering.pred,
      pred.cond = pred.cond
    )

    ## carry out prediction
    if (!model@common_scale) {
      params <- model@covparams
      param.seq <- model@param_sequence
      pred <- vecchia_Mprediction(res, vec.approx.test,
        c(
          params[1:param.seq[1, 2]],
          rep(1, param.seq[2, 2] - param.seq[2, 1] + 1),
          params[param.seq[3, 1]:param.seq[5, 2]]
        ),
        return.values = return.values
      )
    } else {
      pred <- vecchia_Mprediction(res, vec.approx.test, model@covparams,
        return.values = return.values
      )
    }

    # prediction function can return both mean and sds
    # returns a list with elements mu.pred,mu.obs,var.pred,var.obs,V.ord
    Vec.mean.all <- Y_bar + pred$mu.pred + Vecchia.Pred # residual + mean trend
    ndx.out <- NULL
    for (i in seq_along(locs)) {
      ndx.out <- c(ndx.out, rep(i, nrow(locs[[i]])))
    }
    Vec.mean <- list()
    for (i in seq_along(locs)) {
      Vec.mean[[i]] <- Vec.mean.all[ndx.out == i]
    }
    names(Vec.mean) <- names(model@locs_train)
    if (return.values == "mean") {
      return.list <- list(means = Vec.mean)
    } else {
      warning("Variance estimates do not include model fitting variance and are anticonservative. Use with caution.")
      Vec.sds <- list()
      for (i in seq_along(locs)) {
        Vec.sds[[i]] <- sqrt(pred$var.pred[ndx.out == i] +
            model@covparams[model@param_sequence[4, 1] + i - 1])
      }
      return.list <- list(means = Vec.mean, sds = Vec.sds)
    }

    return.list
  }
)

setMethod("check_input", "MultivariateVecchiaModel", function(model, Y, X, locs, Y.names, X.names, center.y, impute.y, lod.upper, lod.lower) {
  if (!is.list(locs)) {
    stop("locs must be a list for multivariate models")
  }
  if (!is.list(Y)) {
    stop("Y must be a list for multivariate models")
  }
  if (!is.list(X)) {
    stop("X must be a list for multivariate models")
  }
  if (length(locs) != length(Y)) {
    stop("locs and Y must have the same length")
  }
  if (length(locs) != length(X)) {
    stop("locs and X must have the same length")
  }
  if (!is.null(names(Y))) {
    names(locs) <- names(Y)
  }
  if (!is.null(Y.names)) {
    if (length(Y.names) != length(locs)) {
      stop("Length of Y.names must match the number of response variables")
    } else {
      names(locs) <- Y.names
    }
  }
  if (is.null(names(locs))) {
    names(locs) <- paste0("Y", seq_along(locs))
  }
  if (!is.null(lod.upper)) {
    if (!is.list(lod.upper)) {
      stop("lod.upper must be a list for multivariate models")
    }
    if (length(locs) != length(lod.upper)) {
      stop("locs and lod.upper must have the same length")
    }
  }
  if (!is.null(lod.lower)) {
    if (!is.list(lod.lower)) {
      stop("lod.lower must be a list for multivariate models")
    }
    if (length(locs) != length(lod.lower)) {
      stop("locs and lod.lower must have the same length")
    }
  }
  for (i in seq_along(locs)) {
    if (!is.matrix(locs[[i]])) {
      stop("Each locs must be a matrix")
    }
    if (i > 1) {
      if (ncol(locs[[i]]) != ncol(locs[[1]])) {
        stop("All locs must have the same number of columns")
      }
    }
    if (!is.matrix(Y[[i]]) & !is.numeric(Y[[i]])) {
      stop("Each Y must be a numeric vector or matrix")
    }
    if (!is.matrix(X[[i]])) {
      stop("Each X must be a matrix")
    }
    if (is.vector(Y[[i]])) {
      Y[[i]] <- as.matrix(Y[[i]])
    }
    if (ncol(Y[[i]]) != 1) {
      stop("Each Y must have only 1 column")
    }
    if (nrow(Y[[i]]) != nrow(locs[[i]])) {
      stop("Each Y must have the same number of rows as locs")
    }
    if (nrow(Y[[i]]) != nrow(X[[i]])) {
      stop("Each Y must have the same number of rows as X")
    }
    if (sum(is.na(X[[i]])) > 0) {
      stop("X must not contain NA's")
    }
    if (sum(is.na(locs[[i]])) > 0) {
      stop("locs must not contain NA's")
    }
    if (sum(is.na(Y[[i]])) > 0 & !impute.y) {
      stop("Y contains NA's and impute.y is FALSE. Set impute.y=TRUE to impute missing Y's.")
    }
    if (!is.null(lod.upper[[i]])) {
      if (!is.numeric(lod.upper[[i]])) {
        stop("Each lod.upper must be numeric")
      }
      if (length(lod.upper[[i]]) != nrow(X[[i]]) &
          length(lod.upper[[i]]) != 1) {
        stop("Length of each lod.upper must equal the number of observations")
      }
    }
    if (!is.null(lod.lower[[i]])) {
      if (!is.numeric(lod.lower[[i]])) {
        stop("Each lod.lower must be numeric")
      }
      if (length(lod.lower[[i]]) != nrow(X[[i]]) &
          length(lod.lower[[i]]) != 1) {
        stop("Length of each lod.lower must equal the number of observations")
      }
    }
  }
  if (!is.null(X.names)) {
    if (!is.list(X.names)) {
      stop("X.names must be a list for multivariate models")
    } else if (length(X.names) != length(locs)) {
      stop("Length of X.names must match the number of response variables")
    } else {
      for (i in seq_along(X.names)) {
        if (length(X.names[[i]]) != ncol(X[[i]])) {
          stop("Length of each X.names must match the number of predictors")
        } else {
          colnames(X[[i]]) <- X.names[[i]]
        }
      }
    }
  }
  Y_bar <- rep(NA, length(Y))
  Y_obs <- NULL
  for (i in seq_along(Y)) {
    Y_obs <- c(Y_obs, !is.na(Y[[i]]))
    if (!is.null(lod.upper)) {
      if (length(lod.upper[[i]]) > 1) {
        Y[[i]][is.na(Y[[i]])] <- lod.upper[[i]][is.na(Y[[i]])]
      } else {
        Y[[i]][is.na(Y[[i]])] <- lod.upper[[i]]
      }
    } else if (!is.null(lod.lower)) {
      if (length(lod.lower[[i]]) > 1) {
        Y[[i]][is.na(Y[[i]])] <- lod.lower[[i]][is.na(Y[[i]])]
      } else {
        Y[[i]][is.na(Y[[i]])] <- lod.lower[[i]]
      }
    }
    if (center.y) {
      Y_bar[i] <- mean(Y[[i]], na.rm = TRUE)
      Y[[i]] <- Y[[i]] - Y_bar[i]
      Y[[i]][is.na(Y[[i]])] <- 0
    } else {
      Y[[i]][is.na(Y[[i]])] <- mean(Y[[i]], na.rm = TRUE)
      Y_bar[i] <- 0
    }
  }
  model@Y_bar <- Y_bar
  model@locs_train <- locs
  model@Y_train <- as.matrix(unlist(Y))
  model@Y_obs <- Y_obs
  model@X_ndx <- ncol(X[[1]]) + 1
  if (length(X) == 1) {
    model@X_train <- X[[1]]
  } else {
    if (is.null(colnames(X[[1]]))) {
      colnames(X[[1]]) <- paste0(names(locs)[1], "_", seq_len(ncol(X[[1]])))
    }
    for (i in 2:length(X)) {
      X[[i]] <- cbind(rep(1, nrow(X[[i]])), X[[i]])
      if (is.null(colnames(X[[i]]))) {
        colnames(X[[i]]) <- paste0(names(locs)[i], "_",
          (seq_len(ncol(X[[i]])) - 1))
      }
      model@X_ndx <- c(model@X_ndx, model@X_ndx[i - 1] + ncol(X[[i]]))
    }
    model@X_train <- psych::superMatrix(X)
  }
  invisible(model)
})

setMethod("check_input_pred", "MultivariateVecchiaModel", function(model, X, locs) {
  if (!is.list(locs)) {
    stop("locs must be a list for multivariate models")
  }
  if (!is.list(X)) {
    stop("X must be a list for multivariate models")
  }
  if (length(locs) != length(X)) {
    stop("locs and X must have the same length")
  }
  if (length(locs) != length(model@locs_train)) {
    stop("Training and test set locs must have the same length")
  }
  for (i in seq_along(locs)) {
    if (!is.matrix(locs[[i]])) {
      stop("Each locs must be a matrix")
    }
    if (ncol(locs[[i]]) != ncol(model@locs_train[[1]])) {
      stop("All locs must have the same number of columns as locs_train")
    }
    if (!is.matrix(X[[i]])) {
      stop("Each X must be a matrix")
    }
    if (nrow(X[[i]]) != nrow(locs[[i]])) {
      stop("Each X must have the same number of rows as locs")
    }
  }
  if (length(X) == 1) {
    X <- X[[1]]
    Y_bar <- rep(model@Y_bar[1], nrow(X))
  } else {
    Y_bar <- NULL
    for (i in seq_along(X)) {
      Y_bar <- c(Y_bar, rep(model@Y_bar[i], nrow(X[[i]])))
    }
    for (i in 2:length(X)) {
      X[[i]] <- cbind(rep(1, nrow(X[[i]])), X[[i]])
    }
    X <- psych::superMatrix(X)
  }
  if (ncol(X) != ncol(model@X_train)) {
    stop("X and X_train must have the same number of predictors")
  }
  list(X = X, Y_bar = Y_bar)
})

setMethod("impute_y", "MultivariateVecchiaModel", function(model) {
  if (model@penalty == "lasso" || model@penalty == "relaxed") {
    # Vecchia prediction at missing values
    Vecchia.Pred <- predict(model@linear_model,
      newx = model@X_train[!model@Y_obs, ], s = "lambda.min",
      gamma = "gamma.min")
    # Vecchia trend prediction at observed data
    Vecchia.hat <- predict(model@linear_model,
      newx = model@X_train[model@Y_obs, ], s = "lambda.min",
      gamma = "gamma.min")
  } else {
    # Vecchia prediction at missing values
    Vecchia.Pred <- predict(model@linear_model, model@X_train[!model@Y_obs, ])
    # Vecchia trend prediction at observed data
    Vecchia.hat <- predict(model@linear_model, model@X_train[model@Y_obs, ])
  }

  # Test set prediction
  res <- model@Y_train[model@Y_obs] - Vecchia.hat

  locs.scaled <- scale_locs(model, model@locs_train)
  all.obs <- model@Y_obs
  locs.otr <- list()
  locs.otst <- list()
  for (i in seq_along(model@locs_train)) {
    nl <- nrow(locs.scaled[[i]])
    cur.obs <- all.obs[1:nl]
    locs.otr[[i]] <- locs.scaled[[i]][cur.obs, ]
    locs.otst[[i]] <- locs.scaled[[i]][!cur.obs, ]
    all.obs <- all.obs[-(1:nl)]
  }

  if (model@common_scale & (length(model@locs_train) > 1)) {
    locs.nd <- eliminate_dupes(locs.otr, locs.otst)
    locs.otr <- locs.nd$locs
    locs.otst <- locs.nd$locs.pred
  }
  vec.approx.test <- vecchia_Mspecify(locs.otr, model@n_neighbors,
    locs.list.pred = locs.otst,
    ordering.pred = "obspred",
    pred.cond = "independent"
  )

  ## carry out prediction
  if (!model@common_scale) {
    params <- model@covparams
    param.seq <- model@param_sequence
    pred <- vecchia_Mprediction(res, vec.approx.test,
      c(
        params[1:param.seq[1, 2]],
        rep(1, length(model@locs_train)),
        params[param.seq[3, 1]:param.seq[5, 2]]
      ),
      return.values = "mean"
    )
  } else {
    pred <- vecchia_Mprediction(res, vec.approx.test, model@covparams,
      return.values = "mean"
    )
  }

  model@Y_train[!model@Y_obs] <- pred$mu.pred + Vecchia.Pred
  invisible(model)
})

setMethod("impute_y_lod", "MultivariateVecchiaModel", function(model, lodu,
  lodl, n.mi = 10, eps = 0.01, maxit = 0, family, nfolds, foldid, parallel,
  cluster, verbose) {
  y <- model@Y_train
  X <- model@X_train
  miss <- !model@Y_obs
  vecchia.approx <- model@vecchia_approx
  params <- model@covparams
  param.seq <- model@param_sequence
  P <- length(model@locs_train)
  if (!model@common_scale) {
    olocs.scaled <- vecchia.approx$locsord
    for (i in 1:vecchia.approx$P) {
      for (j in 1:model@nscale) {
        olocs.scaled[model@vecchia_approx$ondx == i, model@scaling == j] <-
          olocs.scaled[
            model@vecchia_approx$ondx == i,
            model@scaling == j
          ] /
            model@covparams[param.seq[2, 1] + model@nscale * (i - 1) + j - 1]
      }
    }
    vecchia.approx$locsord <- olocs.scaled
    params <- c(params[1:param.seq[1, 2]], rep(1, vecchia.approx$P),
      params[param.seq[3, 1]:param.seq[5, 2]])
  }

  locs.scaled <- scale_locs(model, model@locs_train)
  locs.all <- NULL
  y_ndx <- NULL
  for (i in seq_len(P)) {
    locs.all <- rbind(locs.all, locs.scaled[[i]])
    y_ndx <- c(y_ndx, rep(i, nrow(locs.scaled[[i]])))
  }
  locs.nn <- nn2(locs.all, k = model@n_neighbors + 1)$nn.idx

  Sigma.hat <- array(dim = c(ncol(locs.nn), ncol(locs.nn), sum(miss)))
  k <- 1
  for (i in which(miss)) {
    Sigma.hat[, , k] <- MMatern_cov(locs.all[locs.nn[i, ], , drop = FALSE],
      y_ndx[locs.nn[i, ]], params, P)
    k <- k + 1
  }

  cur.coef <- as.vector(model@beta)
  last.coef <- rep(Inf, ncol(X) + 1)
  itn <- 0
  while (max(abs(cur.coef - last.coef)) > eps & itn < maxit) {
    itn <- itn + 1

    yhat.ni <- X %*% cur.coef[-1]
    yhat.ni <- yhat.ni + mean(y[!miss]) - mean(yhat.ni[!miss])

    coef.mat <- matrix(nrow = n.mi, ncol = (ncol(X) + 1))
    if (parallel) {
      yi <- foreach(i = seq_len(n.mi), .combine = cbind) %dopar% {
        out <- rtmvn_snn2(y - yhat.ni, lodl - yhat.ni,
          lodu - yhat.ni, miss, locs.nn, Sigma.hat)
        out + yhat.ni
      }
    } else {
      yi <- matrix(nrow = length(yhat.ni), ncol = n.mi)
      for (i in seq_len(n.mi)) {
        yi[, i] <- rtmvn_snn2(y - yhat.ni, lodl - yhat.ni,
          lodu - yhat.ni, miss, locs.nn, Sigma.hat)
        yi[, i] <- yi[, i] + yhat.ni
      }
    }
    tiid <- transform_miid(cbind(yi, X), vecchia.approx, params)
    yt <- tiid[, seq_len(ncol(yi))]
    Xt <- as.matrix(tiid[, -(seq_len(ncol(yi)))])

    for (i in seq_len(n.mi)) {
      if (model@penalty == "lasso" || model@penalty == "relaxed") {
        cur.glmnet <- cv.glmnet(as.matrix(Xt), as.matrix(yt[, i]),
          alpha = model@alpha, family = family, nfolds = nfolds,
          foldid = foldid, parallel = parallel,
          relax = model@penalty == "relaxed")
        coef.mat[i, ] <- as.matrix(coef(cur.glmnet, s = "lambda.min",
            gamma = "gamma.min"))
      } else {
        cur.ncvreg <- cv.ncvreg.wrap(as.matrix(Xt), as.matrix(yt[, i]),
          cluster = cluster, foldid = foldid, penalty = model@penalty,
          alpha = model@alpha, family = family, nfolds = nfolds)
        coef.mat[i, ] <- as.matrix(coef(cur.ncvreg))
      }
    }
    last.coef <- cur.coef
    cur.coef <- colMeans(coef.mat)
    if (verbose) {
      cat("LOD imputation iteration", itn, "complete", "\n")
      cat("Current coefficients:", "\n")
      print(cur.coef)
    }
  }
  yhat.ni <- X %*% cur.coef[-1]
  yhat.ni <- yhat.ni + mean(y[!miss]) - mean(yhat.ni[!miss])

  if (parallel) {
    y.na.mat <- foreach(i = seq_len(100), .combine = rbind) %dopar% {
      yi <- rtmvn_snn2(y - yhat.ni, lodl - yhat.ni,
        lodu - yhat.ni, miss, locs.nn, Sigma.hat)
      yi <- yi + yhat.ni
      yi[miss]
    }
  } else {
    y.na.mat <- matrix(nrow = 100, ncol = sum(miss))

    for (i in seq_len(100)) {
      yi <- rtmvn_snn2(y - yhat.ni, lodl - yhat.ni,
        lodu - yhat.ni, miss, locs.nn, Sigma.hat)
      yi <- yi + yhat.ni
      y.na.mat[i, ] <- yi[miss]
    }
  }

  model@Y_train[miss] <- colMeans(y.na.mat)
  invisible(model)
})

setMethod("specify", "MultivariateVecchiaModel", function(model) {
  locs <- model@locs_train
  locs.scaled <- scale_locs(model, locs)
  model@vecchia_approx <- vecchia_Mspecify(locs.scaled, model@n_neighbors)
  if (!model@common_scale) {
    olocs.scaled <- model@vecchia_approx$locsord
    for (i in seq_along(locs)) {
      for (j in 1:model@nscale) {
        olocs.scaled[model@vecchia_approx$ondx == i, model@scaling == j] <-
          olocs.scaled[model@vecchia_approx$ondx == i, model@scaling == j] *
            model@covparams[model@param_sequence[2, 1] + model@nscale * (i - 1) + j - 1]
      }
    }
    model@vecchia_approx$locsord <- olocs.scaled
  }
  invisible(model)
})

#' estimate_theta
#'
#' Estimate covariance parameters during a single round of model optimization
#
#' @param model The model to estimate theta for
#' @param locs the locations matrix
#'
#' @return a model with an updated covariance parameters estimate
#' @noRd
setMethod("estimate_theta", "MultivariateVecchiaModel", function(model, locs, optim.control, method) {
  P <- length(locs)
  if (model@common_scale) {
    vecchia.result <- optim(
      par = model@logparams,
      fn = mvnegloglik,
      vecchia.approx = model@vecchia_approx,
      y = model@res,
      P = P,
      param.seq = model@param_sequence,
      method = method,
      control = optim.control
    )
  } else {
    vecchia.result <- optim(
      par = model@logparams,
      fn = mvnegloglik_ST,
      vecchia.approx = model@vecchia_approx,
      y = model@res,
      P = P,
      param.seq = model@param_sequence,
      scaling = model@scaling,
      nscale = model@nscale,
      method = method,
      control = optim.control
    )
  }

  model@LL_Vecchia_krig <- vecchia.result$value
  model@logparams <- vecchia.result$par
  model <- transform_covariance_parameters(model)
  invisible(model)
})

setMethod("transform_data", "MultivariateVecchiaModel", function(model, Y, X) {
  vecchia.approx <- model@vecchia_approx
  if (!model@common_scale) {
    params <- model@covparams
    param.seq <- model@param_sequence
    olocs.scaled <- vecchia.approx$locsord
    for (i in 1:vecchia.approx$P) {
      for (j in 1:model@nscale) {
        olocs.scaled[model@vecchia_approx$ondx == i, model@scaling == j] <-
          olocs.scaled[
            model@vecchia_approx$ondx == i,
            model@scaling == j
          ] /
            model@covparams[param.seq[2, 1] + model@nscale * (i - 1) + j - 1]
      }
    }
    vecchia.approx$locsord <- olocs.scaled
    transformed.data <- transform_miid(cbind(Y, as.matrix(X)),
      vecchia.approx = vecchia.approx,
      c(
        params[1:param.seq[1, 2]],
        rep(1, vecchia.approx$P),
        params[param.seq[3, 1]:param.seq[5, 2]]
      )
    )
  } else {
    transformed.data <- transform_miid(cbind(Y, as.matrix(X)),
      vecchia.approx = vecchia.approx,
      model@covparams
    )
  }

  model@y_tilde <- Matrix(transformed.data[, 1])
  model@X_tilde <- Matrix(transformed.data[, -1], sparse = FALSE)
  invisible(model)
})
