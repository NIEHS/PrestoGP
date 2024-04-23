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

    # Vecchia prediction at new locations
    Vecchia.Pred <- predict(model@linear_model, newx = X, s = model@linear_model$lambda[model@lambda_1se_idx])
    # Vecchia trend prediction at observed data
    Vecchia.hat <- predict(model@linear_model, newx = model@X_train, s = model@linear_model$lambda[model@lambda_1se_idx])

    # Test set prediction
    res <- model@Y_train - Vecchia.hat

    locs.train.scaled <- scale_locs(model, model@locs_train)
    locs.scaled <- scale_locs(model, locs)
    if (model@apanasovich & (length(model@locs_train) > 1)) {
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
    if (!model@apanasovich) {
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
    Vec.mean <- Y_bar + pred$mu.pred + Vecchia.Pred # residual + mean trend
    if (return.values == "mean") {
      return.list <- list(means = Vec.mean)
    } else {
      warning("Variance estimates do not include model fitting variance and are anticonservative. Use with caution.")
      Vec.sds <- pred$var.pred
      ndx.out <- NULL
      for (i in seq_along(locs)) {
        ndx.out <- c(ndx.out, rep(i, nrow(locs[[i]])))
      }
      for (i in seq_along(locs)) {
        Vec.sds[ndx.out == i] <- sqrt(Vec.sds[ndx.out == i] +
            model@covparams[model@param_sequence[4, i]])
      }
      return.list <- list(means = Vec.mean, sds = Vec.sds)
    }

    return(return.list)
  }
)

setMethod("check_input", "MultivariateVecchiaModel", function(model, Y, X, locs, center.y, impute.y, lod) {
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
  if (!is.null(lod)) {
    if (!is.list(lod)) {
      stop("lod must be a list for multivariate models")
    }
    if (length(locs) != length(lod)) {
      stop("locs and lod must have the same length")
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
    if (!is.null(lod[[i]])) {
      if (!is.numeric(lod[[i]])) {
        stop("Each lod must be numeric")
      }
      if (length(lod[[i]]) != nrow(X[[i]]) & length(lod[[i]]) != 1) {
        stop("Length of each lod must equal the number of observations")
      }
    }
  }
  Y_bar <- rep(NA, length(Y))
  Y_obs <- NULL
  for (i in seq_along(Y)) {
    Y_obs <- c(Y_obs, !is.na(Y[[i]]))
    if (!is.null(lod)) {
      Y[[i]][is.na(Y[[i]])] <- mean(lod[[i]]) / 2
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
  if (length(X) == 1) {
    model@X_train <- X[[1]]
  } else {
    for (i in 2:length(X)) {
      X[[i]] <- cbind(rep(1, nrow(X[[i]])), X[[i]])
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
  return(list(X = X, Y_bar = Y_bar))
})

setMethod("impute_y", "MultivariateVecchiaModel", function(model) {
  if (sum(!model@Y_obs) > 0) {
    Vecchia.Pred <- predict(model@linear_model,
      newx = model@X_train[!model@Y_obs, ],
      s = model@linear_model$lambda[model@lambda_1se_idx])
    Vecchia.hat <- predict(model@linear_model,
      newx = model@X_train[model@Y_obs, ],
      s = model@linear_model$lambda[model@lambda_1se_idx])

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

    if (model@apanasovich & (length(model@locs_train) > 1)) {
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
    if (!model@apanasovich) {
      params <- model@covparams
      param.seq <- model@param_sequence
      pred <- vecchia_Mprediction(res, vec.approx.test,
        c(
          params[1:param.seq[1, 2]],
          rep(1, param.seq[2, 2] - param.seq[2, 1] + 1),
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
  }
  invisible(model)
})

setMethod("specify", "MultivariateVecchiaModel", function(model) {
  locs <- model@locs_train
  locs.scaled <- scale_locs(model, locs)
  model@vecchia_approx <- vecchia_Mspecify(locs.scaled, model@n_neighbors)
  if (!model@apanasovich) {
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
  if (model@apanasovich) {
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
  if (!model@apanasovich) {
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
        rep(1, param.seq[2, 2] - param.seq[2, 1] + 1),
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

setMethod("theta_names", "MultivariateVecchiaModel", function(model) {
  c("Marginal Variance", "Range", "Smoothness", "Nugget")
})
