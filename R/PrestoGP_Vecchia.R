#' Univariate Vecchia PrestoGP model class
#'
#' This class is used to create univariate models with a likelihood function
#' conditioned on a subset of the observations (i.e., Vecchia models). See
#' \code{\link{PrestoGPModel}} for a description of the slots.
#'
#' @seealso \code{\link{PrestoGPModel}}
#'
#' @export VecchiaModel
#'
#' @examples
#' pgp.model <- new("VecchiaModel", n_neighbors = 25)
#' @include PrestoGP_Model.R
VecchiaModel <- setClass("VecchiaModel",
  contains = "PrestoGPModel"
)

validityVecchiaModel <- function(object) {
  TRUE
}
setValidity("VecchiaModel", validityVecchiaModel)

setMethod("initialize", "VecchiaModel", function(.Object, n_neighbors = 25, ...) {
  .Object@n_neighbors <- n_neighbors
  .Object@min_m <- 3
  .Object <- callNextMethod()
  validObject(.Object)
  .Object
})

#' @rdname get_Y
setMethod("get_Y", "VecchiaModel",
  function(model) {
    as.vector(model@Y_train + model@Y_bar)
  }
)

#' @rdname prestogp_predict
setMethod("prestogp_predict", "VecchiaModel",
  function(model, X, locs, m = NULL, ordering.pred = c("obspred", "general"), pred.cond = c("independent", "general"), return.values = c("mean", "meanvar")) {
    # validate parameters
    ordering.pred <- match.arg(ordering.pred)
    pred.cond <- match.arg(pred.cond)
    return.values <- match.arg(return.values)
    model <- check_input_pred(model, X, locs)
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

    locs.train.scaled <- scale_locs(model, model@locs_train)[[1]]
    locs.scaled <- scale_locs(model, list(locs))[[1]]
    vec.approx.test <- vecchia_specify(locs.train.scaled, m, locs.pred = locs.scaled, ordering.pred = ordering.pred, pred.cond = pred.cond)

    ## carry out prediction
    if (!model@common_scale) {
      pred <- vecchia_prediction(
        res,
        vec.approx.test,
        c(model@covparams[1], 1, model@covparams[3]),
        model@covparams[4],
        return.values = return.values
      )
    } else {
      pred <- vecchia_prediction(
        res,
        vec.approx.test,
        c(model@covparams[1], model@covparams[2], model@covparams[3]),
        model@covparams[4], return.values = return.values
      )
    }

    # prediction function can return both mean and sds
    # returns a list with elements mu.pred,mu.obs,var.pred,var.obs,V.ord
    # residual + mean trend
    Vec.mean <- model@Y_bar + pred$mu.pred + Vecchia.Pred
    if (return.values == "mean") {
      return.list <- list(means = Vec.mean)
    } else {
      warning("Variance estimates do not include model fitting variance and are anticonservative. Use with caution.")
      Vec.sds <- sqrt(pred$var.pred + model@covparams[4])
      return.list <- list(means = Vec.mean, sds = Vec.sds)
    }

    return.list
  }
)

setMethod("check_input", "VecchiaModel", function(model, Y, X, locs, Y.names, X.names, center.y, impute.y, lod.upper, lod.lower) {
  if (!is.matrix(locs)) {
    stop("locs must be a matrix")
  }
  if (!is.matrix(Y) & !is.numeric(Y)) {
    stop("Y must be a numeric vector or matrix")
  }
  if (!is.matrix(X)) {
    stop("X must be a matrix")
  }
  if (is.vector(Y)) {
    Y <- as.matrix(Y)
  }
  if (ncol(Y) != 1) {
    stop("Y must have only 1 column")
  }
  if (nrow(Y) != nrow(locs)) {
    stop("Y must have the same number of rows as locs")
  }
  if (nrow(Y) != nrow(X)) {
    stop("Y must have the same number of rows as X")
  }
  if (sum(is.na(X)) > 0) {
    stop("X must not contain NA's")
  }
  if (sum(is.na(locs)) > 0) {
    stop("locs must not contain NA's")
  }
  if (sum(is.na(Y)) > 0 & !impute.y) {
    stop("Y contains NA's and impute.y is FALSE. Set impute.y=TRUE to impute missing Y's.")
  }
  if (!is.null(lod.upper)) {
    if (!is.numeric(lod.upper)) {
      stop("lod.upper must be numeric")
    }
    if (length(lod.upper) != nrow(X) & length(lod.upper) != 1) {
      stop("Length of lod.upper must equal the number of observations")
    }
  }
  if (!is.null(lod.lower)) {
    if (!is.numeric(lod.lower)) {
      stop("lod.lower must be numeric")
    }
    if (length(lod.lower) != nrow(X) & length(lod.lower) != 1) {
      stop("Length of lod.lower must equal the number of observations")
    }
  }
  model@Y_obs <- !is.na(as.vector(Y))
  if (!is.null(lod.upper)) {
    if (length(lod.upper) > 1) {
      Y[!model@Y_obs] <- lod.upper[!model@Y_obs]
    } else {
      Y[!model@Y_obs] <- lod.upper
    }
  } else if (!is.null(lod.lower)) {
    if (length(lod.lower) > 1) {
      Y[!model@Y_obs] <- lod.lower[!model@Y_obs]
    } else {
      Y[!model@Y_obs] <- lod.lower
    }
  }
  if (center.y) {
    model@Y_bar <- mean(Y, na.rm = TRUE)
    Y <- Y - model@Y_bar
    Y[is.na(Y)] <- 0
  } else {
    model@Y_bar <- 0
    Y[is.na(Y)] <- mean(Y, na.rm = TRUE)
  }
  model@X_train <- X
  model@X_ndx <- ncol(X) + 1
  model@Y_train <- Y
  model@locs_train <- list(locs)
  if (!is.null(colnames(Y))) {
    names(model@locs_train) <- colnames(Y)
  }
  if (!is.null(Y.names)) {
    if (length(Y.names) > 1) {
      stop("Length of Y.names must match the number of response variables")
    } else {
      names(model@locs_train) <- Y.names
    }
  }
  if (!is.null(X.names)) {
    if (length(X.names) != ncol(model@X_train)) {
      stop("Length of X.names must match the number of predictor variables")
    } else {
      colnames(model@X_train) <- X.names
    }
  }
  invisible(model)
})

setMethod("check_input_pred", "VecchiaModel", function(model, X, locs) {
  if (!is.matrix(locs)) {
    stop("locs must be a matrix")
  }
  if (!is.matrix(X)) {
    stop("X must be a matrix")
  }
  if (ncol(locs) != ncol(model@locs_train[[1]])) {
    stop("locs must have the same number of columns as locs_train")
  }
  if (nrow(X) != nrow(locs)) {
    stop("X must have the same number of rows as locs")
  }
  if (ncol(X) != ncol(model@X_train)) {
    stop("X and X_train must have the same number of predictors")
  }
  invisible(model)
})

setMethod("impute_y", "VecchiaModel", function(model) {
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

  locs.scaled <- scale_locs(model, model@locs_train)[[1]]
  vec.approx.test <- vecchia_specify(locs.scaled[model@Y_obs, ],
    model@n_neighbors,
    locs.pred = locs.scaled[!model@Y_obs, ],
    ordering.pred = "obspred", pred.cond = "independent")

  ## carry out prediction
  if (!model@common_scale) {
    pred <- vecchia_prediction(
      res,
      vec.approx.test,
      c(model@covparams[1], 1, model@covparams[3]),
      model@covparams[4],
      return.values = "mean"
    )
  } else {
    pred <- vecchia_prediction(
      res,
      vec.approx.test,
      c(model@covparams[1], model@covparams[2], model@covparams[3]),
      model@covparams[4], return.values = "mean"
    )
  }

  model@Y_train[!model@Y_obs] <- pred$mu.pred + Vecchia.Pred
  invisible(model)
})

setMethod("impute_y_lod", "VecchiaModel", function(model, lodu, lodl, n.mi = 10,
  eps = 0.01, maxit = 1, family, nfolds, foldid, parallel, cluster, verbose) {
  y <- model@Y_train
  X <- model@X_train
  miss <- !model@Y_obs
  vecchia.approx <- model@vecchia_approx
  params <- model@covparams
  if (!model@common_scale) {
    vecchia.approx$locsord <- scale_locs(
      model,
      list(vecchia.approx$locsord)
    )[[1]]
    params <- c(params[1], 1, params[3:4])
  }

  locs.scaled <- scale_locs(model, model@locs_train)[[1]]
  locs.nn <- nn2(locs.scaled, k = model@n_neighbors + 1)$nn.idx

  Sigma.hat <- array(dim = c(ncol(locs.nn), ncol(locs.nn), sum(miss)))
  k <- 1
  for (i in which(miss)) {
    Sigma.hat[, , k] <- MMatern_cov(locs.scaled[locs.nn[i, ], , drop = FALSE],
      rep(1, ncol(locs.nn)), params, 1)
    k <- k + 1
  }

  cur.coef <- as.vector(model@beta)
  last.coef <- rep(Inf, ncol(X) + 1)
  itn <- 0
  while (max(abs(cur.coef - last.coef)) > eps & itn <= maxit) {
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
    tiid <- transform_iid(cbind(yi, X), vecchia.approx, params[1:3],
      params[4])
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

#' specify
#'
#' Specify the conditioning set using m nearest neighbors.
#'
#' @param model The model to specify
#'
#' @return a model with a specified conditioning set
#' @noRd
setMethod("specify", "VecchiaModel", function(model) {
  locs.scaled <- scale_locs(model, model@locs_train)
  model@vecchia_approx <- vecchia_specify(locs.scaled[[1]], model@n_neighbors)
  if (!model@common_scale) {
    olocs.scaled <- model@vecchia_approx$locsord
    for (j in 1:model@nscale) {
      olocs.scaled[, model@scaling == j] <- olocs.scaled[, model@scaling == j] *
        model@covparams[model@param_sequence[2, 1] + j - 1]
    }
    model@vecchia_approx$locsord <- olocs.scaled
  }
  invisible(model)
})

#' estimate_theta
#'
#' Estimate covariance parameters during a single round of model optimization
#'
#' @param model The model to estimate theta for
#' @param locs the locations matrix
#'
#' @return a model with an updated covariance parameters estimate
#' @noRd
setMethod("estimate_theta", "VecchiaModel", function(model, locs, optim.control, method) {
  if (model@common_scale) {
    vecchia.result <- optim(
      par = model@logparams,
      fn = negloglik_vecchia,
      res = model@res,
      vecchia.approx = model@vecchia_approx,
      param.seq = model@param_sequence,
      method = method,
      control = optim.control
    )
  } else {
    vecchia.result <- optim(
      par = model@logparams,
      fn = negloglik_vecchia_ST,
      res = model@res,
      vecchia.approx = model@vecchia_approx,
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

#' transform_data
#'
#' Transform data to be independent, identically distributed
#'
#' @param model The model to estimate theta for
#' @param Y the dependent variable matrix
#' @param X the independent variable matrix
#'
#' @return a model with i.i.d. data
#' @noRd
setMethod("transform_data", "VecchiaModel", function(model, Y, X) {
  vecchia.approx <- model@vecchia_approx
  params <- model@covparams
  if (!model@common_scale) {
    vecchia.approx$locsord <- scale_locs(
      model,
      list(vecchia.approx$locsord)
    )[[1]]
    transformed.data <- transform_iid(cbind(Y, as.matrix(X)),
      vecchia.approx = vecchia.approx, c(params[1], 1, params[3]),
      params[4])
  } else {
    transformed.data <- transform_iid(cbind(Y, as.matrix(X)),
      vecchia.approx = vecchia.approx,
      params[1:3], params[4]
    )
  }

  model@y_tilde <- Matrix(transformed.data[, 1], sparse = FALSE)
  model@X_tilde <- Matrix(transformed.data[, -1], sparse = FALSE)
  invisible(model)
})
