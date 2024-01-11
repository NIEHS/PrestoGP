#' Model created with a likelihood function conditioned on a subset of the observations.
#'
#' @slot model PrestoGPModel.
#'
#' @export VecchiaModel
#'
#' @examples
#' @include PrestoGP_Model.R
#' @noRd
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

#' Make predictions using a previously fit spatiotemporal model.
#'
#' @param model A model object returned by prestogp_fit.
#' @param X Independent variable matrix for prediction.
#' @param locs Locations matrix for which the value is predicted.
#' @param m The number of neighbors to condition on (if not provided, training value of m will be used).
#'
#' @return predicted values for the dependent variable
#' @export
#'
#' @examples
#'
#' ...
#' model <- VecchiaModel()
#' model <- prestogp_fit(model, logNO2, X, locs)
#' prediction <- prestogp_predict(model, X.test, locs.test)
#' Vec.mean <- prediction[[1]]
#' Vec.sds <- prediction[[2]]
setMethod("prestogp_predict", "VecchiaModel", function(model, X, locs, m = NULL, ordering.pred = c("obspred", "general"), pred.cond = c("independent", "general"), return.values = c("mean", "meanvar")) {
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

  # Vecchia prediction at new locations
  # Vecchia.Pred <- predict(model@Vecchia_SCAD_fit[[1]], X = X, which = model@lambda_1se_idx[[1]])
  Vecchia.Pred <- predict(model@linear_model, newx = X, s = model@linear_model$lambda[model@lambda_1se_idx])
  # Vecchia trend prediction at observed data
  # Vecchia.hat <- predict(model@Vecchia_SCAD_fit[[1]], X = model@X_train, which = model@lambda_1se_idx[[1]])
  Vecchia.hat <- predict(model@linear_model, newx = model@X_train, s = model@linear_model$lambda[model@lambda_1se_idx])

  # Test set prediction
  res <- model@Y_train - Vecchia.hat

  locs.train.scaled <- scale_locs(model, model@locs_train)[[1]]
  locs.scaled <- scale_locs(model, list(locs))[[1]]
  vec.approx.test <- vecchia_specify(locs.train.scaled, m, locs.pred = locs.scaled, ordering.pred=ordering.pred, pred.cond=pred.cond)

    ## carry out prediction
    if (!model@apanasovich) {
        pred <- vecchia_prediction(res, vec.approx.test, c(model@covparams[1], 1, model@covparams[3]), model@covparams[4], return.values=return.values)
    } else {
        pred <- vecchia_prediction(res, vec.approx.test, c(model@covparams[1], model@covparams[2], model@covparams[3]), model@covparams[4], return.values=return.values)
    }

  # prediction function can return both mean and sds
  # returns a list with elements mu.pred,mu.obs,var.pred,var.obs,V.ord
  Vec.mean <- pred$mu.pred + Vecchia.Pred # residual + mean trend
  if (return.values == "mean") {
      return.list <- list(means = Vec.mean)
  } else {
      warning("Variance estimates do not include model fitting variance and are anticonservative. Use with caution.")
      Vec.sds <- sqrt(pred$var.pred + model@covparams[4])
      return.list <- list(means = Vec.mean, sds = vec.sds)
  }

  return(return.list)
})

setMethod("check_input", "VecchiaModel", function(model, Y, X, locs) {
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
  model@X_train <- X
  model@Y_train <- Y
  model@locs_train <- list(locs)
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

#' specify
#'
#' Specify the conditioning set using m nearest neighbors.
#'
#' @param model The model to specify
#' @param locs the locations matrix
#' @param m the number of neighbors to condition on
#'
#' @return a model with a specified conditioning set
setMethod("specify", "VecchiaModel", function(model) {
  locs.scaled <- scale_locs(model, model@locs_train)
  model@vecchia_approx <- vecchia_specify(locs.scaled[[1]], model@n_neighbors)
  if (!model@apanasovich) {
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
setMethod("estimate_theta", "VecchiaModel", function(model, locs, optim.control, method) {
  if (model@apanasovich) {
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
setMethod("transform_data", "VecchiaModel", function(model, Y, X) {
  vecchia.approx <- model@vecchia_approx
  params <- model@covparams
  if (!model@apanasovich) {
    param.seq <- model@param_sequence
    vecchia.approx$locsord <- scale_locs(
      model,
      list(vecchia.approx$locsord)
    )[[1]]
    transformed.data <- transform_iid(cbind(Y, as.matrix(X)),
      vecchia.approx = vecchia.approx,
      c(
        params[1:param.seq[1, 2]],
        rep(
          1,
          param.seq[2, 2] - param.seq[2, 1] + 1
        ),
        params[param.seq[3, 1]]
      ),
      params[param.seq[4, 1]]
    )
  } else {
    transformed.data <- transform_iid(cbind(Y, as.matrix(X)),
      vecchia.approx = vecchia.approx,
      params[1:3], params[4]
    )
  }

  model@y_tilde <- Matrix(transformed.data[, 1])
  model@X_tilde <- Matrix(transformed.data[, -1], sparse = FALSE)
  invisible(model)
})

#' theta_names
#'
#' Return a vector specifying the names of the different covariance parameters (used for show method)
#'
#' @param model The model to estimate theta for
#'
#' @return a vector with the namems of the covariance parameterss
setMethod("theta_names", "VecchiaModel", function(model) {
  c("Marginal Variance", "Spatial Range", "Temporal Range", "Nugget")
})
