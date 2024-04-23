#' Univariate full PrestoGP model class
#'
#' This class is used to create univariate models created with a likelihood
#' function conditioned on all observations. See \code{\link{PrestoGPModel}}
#' for a description of the slots.
#'
#' @seealso \code{\link{PrestoGPModel}}
#'
#' @export FullModel
#'
#' @examples
#' pgp.model <- new("FullModel")
#' @include PrestoGP_Vecchia.R
FullModel <- setClass("FullModel",
  contains = "VecchiaModel"
)

validityFullModel <- function(object) {
  TRUE
}
setValidity("FullModel", validityFullModel)

setMethod("initialize", "FullModel", function(.Object, n_neighbors = 0, ...) {
  .Object <- callNextMethod()
  .Object@n_neighbors <- 0
  .Object@min_m <- 0
  validObject(.Object)
  .Object
})

#' @rdname prestogp_predict
setMethod("prestogp_predict", "FullModel", function(model,
  X, locs, m = NULL, ordering.pred = c("obspred", "general"), pred.cond = c("independent", "general"), return.values = c("mean", "meanvar")) {
  stop("Prediction is not currently supported for full models")
})

setMethod("specify", "FullModel", function(model) {
  invisible(model)
})

setMethod("compute_residuals", "FullModel", function(model, Y, Y.hat) {
  model@res <- as.double(Y - Y.hat)
  invisible(model)
})

setMethod("estimate_theta", "FullModel", function(model, locs, optim.control, method) {
  if (model@apanasovich) {
    full.result <- optim(
      par = model@logparams,
      fn = negloglik.full,
      d = fields::rdist(locs),
      y = model@res,
      param.seq = model@param_sequence,
      method = method,
      control = optim.control
    )
  } else {
    full.result <- optim(
      par = model@logparams,
      fn = negloglik_full_ST,
      locs = locs,
      y = model@res,
      param.seq = model@param_sequence,
      scaling = model@scaling,
      nscale = model@nscale,
      method = method,
      control = optim.control
    )
  }

  model@LL_Vecchia_krig <- full.result$value
  model@logparams <- full.result$par
  model <- transform_covariance_parameters(model)
  invisible(model)
})

setMethod("transform_data", "FullModel", function(model, Y, X) {
  n <- nrow(model@Y_train)
  params <- model@covparams
  locs.scaled <- scale_locs(model, model@locs_train)[[1]]
  if (!model@apanasovich) {
    param.seq <- model@param_sequence
    Omega.full <- params[1] * Matern(rdist(locs.scaled),
      range = 1,
      smoothness = params[param.seq[3, 1]]
    ) +
      params[param.seq[4, 1]] * diag(n)
  } else {
    Omega.full <- params[1] * Matern(rdist(locs.scaled),
      range = params[2],
      smoothness = params[3]
    ) +
      params[4] * diag(n)
  }

  Omega.lc <- t(chol(Omega.full))
  model@y_tilde <- Matrix(solve(Omega.lc, Y))
  model@X_tilde <- Matrix(solve(Omega.lc, X), sparse = FALSE)
  invisible(model)
})
