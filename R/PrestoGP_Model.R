#' PrestoGPModel superclass
#'
#' This is the superclass for PrestoGP models. All other types of PrestoGP
#' model classes (e.g. \code{\link{VecchiaModel-class}} and
#' \code{\link{FullModel-class}} are inherited from this class. Normally
#' users should not create objects of this class. Instead, they should use
#' the appropriate inherited class for the type of model they are fitting.
#'
#' @slot covparams A numeric vector containing the parameters for the Matern
#' model.
#' @slot logparams A numeric vector containing the transformed versions of the
#' Matern parameters (used internally for likelihood calculations).
#' @slot beta A column matrix containing the regression coefficients.
#' @slot vecchia_approx The output of the Vecchia specify function. See
#' \code{\link[GPvecchia]{vecchia_specify}} and \code{\link{vecchia_Mspecify}}.
#' @slot X_train A matrix containing the original predictors. This will be a
#' "super matrix" for multivariate models. See
#' \code{\link[psych]{superMatrix}}.
#' @slot Y_train A column matrix containing the original response values.
#' After fitting the model, missing values will be replaced with the imputed
#' values if impute.y is TRUE. See \code{link{prestogp_fit}}.
#' @slot X_ndx A vector used to find the elements of beta corresponding to
#' specific outcomes. The \emph{i}th element of X_ndx is the index of the
#' last element of beta corresponding to predictors for outcome \emph{i}.
#' @slot Y_bar A vector containing the means of each outcome.
#' @slot Y_obs A logical vector used to track which values of Y are non-missing.
#' @slot X_tilde The matrix of transformed predictors.
#' @slot y_tilde The column matrix containing the transformed response values.
#' @slot res A numeric vector of the residuals.
#' @slot locs_train A list containing the location coordinates. Each element
#' of the list corresponds to a different outcome. (The list will have length
#' 1 for univariate models.)
#' @slot penalty The type of penalized regression used. Should be one
#' of "lasso", "relaxed", "MCP", or "SCAD".
#' @slot linear_model The glmnet or ncvreg model. See
#' \code{\link[glmnet]{glmnet}}, \code{\link[ncvreg]{ncvreg}},
#' \code{\link[glmnet]{cv.glmnet}}, and \code{\link[ncvreg]{cv.ncvreg}}.
#' @slot converged Did the model fitting process converge (boolean)?
#' @slot LL_Vecchia_krig The value of the negative log likelihood function
#' after optimization.
#' @slot error Penalized model error. See References for details.
#' @slot n_neighbors Number of neighbors to condition on for the Vecchia
#' approximation. Ignored for full models.
#' @slot min_m Minimum permissible number of neighbors.
#' @slot alpha Parameter alpha for glmnet or ncvreg. See
#' \code{\link[glmnet]{glmnet}} or \code{\link[ncvreg]{ncvreg}}.
#' @slot scaling The indices of the scale parameters. See
#' \code{link{prestogp_fit}}.
#' @slot nscale The number of scale parameters in the model.
#' @slot common_scale Do all columns of locs have the same scale parameter?
#' @slot param_sequence Records the indices of the various Matern parameters.
#' See \code{\link{create_param_sequence}}.
#'
#' @seealso \code{\link{VecchiaModel-class}}, \code{\link{FullModel-class}},
#' \code{\link{MultivariateVecchiaModel-class}}, \code{\link{prestogp_fit}}
#'
#' @references
#' \itemize{
#' \item Apanasovich, T.V., Genton, M.G. and Sun, Y. "A valid Matérn class of
#' cross-covariance functions for multivariate random fields with any number
#' of components", Journal of the American Statistical Association (2012)
#' 107(497):180-193.
#' \item Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
#' land-use regression for ground-level nitrogen dioxide", The Annals of
#' Applied Statistics (2021) 15(2):688-710.
#' }
#'
#' @examples
#' pgp.vmodel <- new("VecchiaModel", n_neighbors = 25)
#' pgp.fmodel <- new("FullModel")
#' pgp.mmodel <- new("MultivariateVecchiaModel", n_neighbors = 25)
#'
#' @export PrestoGPModel
PrestoGPModel <- setClass("PrestoGPModel",
  slots = list(
    covparams = "numeric",
    beta = "matrix", # "numeric", #the beta matrix
    vecchia_approx = "list", # the output of vecchia_specify
    y_tilde = "dgeMatrix", # iid transformed dependent variables matrix
    X_tilde = "dgeMatrix", # iid transformed independent variables matrix
    res = "numeric", # residuals
    penalty = "character", # type of model (lasso, SCAD, etc.)
    linear_model = "ANY", # the linear model
    X_train = "matrix", # the original independent variable matrix
    Y_train = "matrix", # the original dependent variable matrix
    X_ndx = "numeric", # used to map X_train to elements of the original X's
    Y_bar = "numeric", # the means of each outcome variable
    Y_obs = "logical", # a logical vector indicating whether the ith observation was observed
    locs_train = "list", # the location / temporal matrix
    converged = "logical", # a logical variable that is true if the model fitting process has converged
    LL_Vecchia_krig = "numeric", # the value of the negative log likelihood function after optimization
    error = "numeric", # negative log likelihood + SCAD penalty likelihood
    n_neighbors = "numeric", # the number of neighbors to condition on for the Vecchia approximation
    min_m = "numeric", # the minimum m required by the specific model type (full vs Vecchia)
    alpha = "numeric", # the alpha ratio of ridge to lasso penalty
    scaling = "numeric", # the indices of the scale parameters,
    nscale = "numeric", # the number of scale parameters
    common_scale = "logical", # is there a common scale parameter?
    param_sequence = "matrix", # maps the indices of the various Matern parameters
    logparams = "numeric" # transformed version of the Matern parameters
  )
)


validityPrestoGPModel <- function(object) {
  #  if(object@n_neighbors < object@min_m){
  #    stop(paste("N_neighbors must be at least ", object@min_m, ".", sep=""))
  #  }
  TRUE
}
setValidity("PrestoGPModel", validityPrestoGPModel)

setMethod("initialize", "PrestoGPModel", function(.Object, ...) {
  #  .Object@linear_model <- structure(list(), class = "cv.glmnet")
  .Object@alpha <- 1 # 0.5
  .Object <- callNextMethod()
  validObject(.Object)
  .Object
})

#' Extract the Y values for a PrestoGP model.
#'
#' This method extracts the values of Y for a PrestoGP model. If imputation
#' was used when fitting the model, the missing Y's will be replaced with
#' their imputed values.
#'
#' @param model The PrestoGP model object
#'
#' @return A vector or list containing the values of Y. Any missing Y's will
#' be replaced with their imputed values.
#'
#' @seealso \code{\link{PrestoGPModel-class}}, \code{\link{prestogp_fit}}
#'
#' @references
#' \itemize{
#' \item Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
#' land-use regression for ground-level nitrogen dioxide", The Annals of
#' Applied Statistics (2021) 15(2):688-710.
#' }
#'
#' @rdname get_Y
#' @export
#'
#' @examples
#' data(soil)
#' soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
#' y <- soil[,4]                   # predict moisture content
#' X <- as.matrix(soil[,5:9])
#' locs <- as.matrix(soil[,1:2])
#'
#' soil.vm <- new("VecchiaModel", n_neighbors = 10)
#' soil.vm <- prestogp_fit(soil.vm, y, X, locs)
#' get_Y(soil.vm)
setGeneric("get_Y", function(model) standardGeneric("get_Y"))
setGeneric("get_theta", function(model) standardGeneric("get_theta"))
setGeneric("get_beta", function(model) standardGeneric("get_beta"))
setGeneric(
  "get_linear_model", function(model) standardGeneric("get_linear_model"))
setGeneric("get_neighbors", function(model) standardGeneric("get_neighbors"))
setGeneric("get_scaling", function(model) standardGeneric("get_scaling"))
setGeneric("get_converged", function(model) standardGeneric("get_converged"))
setGeneric("get_pen_loglik", function(model) standardGeneric("get_pen_loglik"))
setGeneric("plot_beta", function(model, ...) standardGeneric("plot_beta"))
setGeneric(
  "prestogp_fit",
  function(model, Y, X, locs, Y.names = NULL, X.names = NULL, scaling = NULL,
    common_scale = NULL, covparams = NULL, beta.hat = NULL, tol = 0.999999,
    max_iters = 100, center.y = NULL, impute.y = FALSE, lod = NULL,
    quiet = FALSE, verbose = FALSE, optim.method = "Nelder-Mead",
    optim.control = list(trace = 0, reltol = 1e-3, maxit = 5000),
    penalty = c("lasso", "relaxed", "MCP", "SCAD"), alpha = 1,
    family = c("gaussian", "binomial"), nfolds = 10, foldid = NULL,
    parallel = FALSE, cluster = NULL, adaptive = FALSE) {
    standardGeneric("prestogp_fit")
  }
)

#' Prediction for PrestoGP models
#'
#' After fitting a PrestoGP model, this method can be used to make predictions
#' on an independent test set (consisting of both new locations and new values
#' of the predictor variables).
#'
#' @param model A PrestoGP model object obtained after running
#' \code{link{prestogp_fit}}.
#' @param X The values of the predictor variable(s) for which prediction will
#' be performed. Should be a matrix for univariate models or a list for
#' multivariate models.
#' @param locs The locations where prediction will be performed. Should be a
#' matrix for univariate models or a list for multivariate models.
#' @param m The number of neighbors to condition on. If not specified, it will
#' default to the value of m used to fit the model.
#' @param ordering.pred Should "obspred" or "general" ordering be used for
#' prediction? See \code{\link[GPvecchia]{vecchia_specify}} or
#' \code{\link{vecchia_Mspecify}}. Defaults to "obspred".
#' @param pred.cond Should prediction conditioning be "general" or
#' "indepedent"? See \code{\link[GPvecchia]{vecchia_specify}} or
#' \code{\link{vecchia_Mspecify}}. Defaults to "independent".
#' @param return.values Values that should be returned. Possible values
#' include "mean" and "meanvar". See
#' \code{\link[GPvecchia]{vecchia_prediction}} or
#' \code{\link{vecchia_Mprediction}}. Defaults to "mean".
#'
#' @details It is important to note that the variance estimates produced by
#' this function assume that the Matern covariance parameters are fixed and
#' known. It does not consider any variance resulting from estimating the
#' covariance parameters or the regression coefficients in PrestoGP models.
#' These variance estimates will be anticonservative in such cases (and a
#' warning will be returned when these estimates are calculated).
#'
#' Prediction is currently not implemented for full models. This function
#' will return an error if it is applied to a full model.
#'
#' @return A list containing the estimated mean values and (if requested)
#' standard deviations for the predictions.
#'
#' @seealso \code{\link{PrestoGPModel-class}},\code{\link{prestogp_fit}},
#' \code{\link[GPvecchia]{vecchia_prediction}},
#' \code{\link{vecchia_Mprediction}}
#'
#' @references
#' \itemize{
#' \item Katzfuss, M., and Guinness, J. "A general framework for Vecchia
#' approximations of Gaussian processes", Statistical Science (2021)
#' 36(1):124-141.
#' \item Katzfuss, M., Guinness, J., Gong, W. and Zilber, D. "Vecchia
#' approximations of Gaussian-process predictions", Journal of Agricultural,
#' Biological and Environmental Statistics (2020) 25:383-414.
#' \item Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
#' land-use regression for ground-level nitrogen dioxide", The Annals of
#' Applied Statistics (2021) 15(2):688-710.
#' }
#'
#' @rdname prestogp_predict
#' @export
#'
#' @examples
#' data(soil)
#' soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
#' y <- soil[,4]                   # predict moisture content
#' X <- as.matrix(soil[,5:9])
#' locs <- as.matrix(soil[,1:2])
#'
#' # Create training and test sets
#' n <- length(y)
#' otr <- rep(FALSE, n)
#' otr[sample(1:n, size=floor(n/2))] <- TRUE
#' otst <- !otr
#'
#' # Fit the model on the training set
#' soil.vm <- new("VecchiaModel", n_neighbors = 10)
#' soil.vm <- prestogp_fit(soil.vm, y[otr], X[otr,], locs[otr,])
#'
#' # Perform predictions on the test set
#' soil.yhat <- prestogp_predict(soil.vm, X[otst,], locs[otst,])
setGeneric(
  "prestogp_predict",
  function(model, X = "matrix", locs = "matrix", m = "numeric",
    ordering.pred = c("obspred", "general"),
    pred.cond = c("independent", "general"),
    return.values = c("mean", "meanvar")) {
    standardGeneric("prestogp_predict")
  }
)
setGeneric("calc_covparams",
  function(model, locs, Y, covparams) standardGeneric("calc_covparams"))
setGeneric("specify", function(model, ...) standardGeneric("specify"))
setGeneric("compute_residuals",
  function(model, Y, Y.hat, family) standardGeneric("compute_residuals"))
setGeneric("transform_data",
  function(model, Y, X) standardGeneric("transform_data"))
setGeneric("estimate_theta",
  function(model, locs, optim.control, method) {
    standardGeneric("estimate_theta")})
setGeneric("estimate_betas",
  function(model, family, nfolds, foldid, parallel, cluster, penalty,
    adaptive) {
    standardGeneric("estimate_betas")})
setGeneric("impute_y", function(model) standardGeneric("impute_y"))
setGeneric("impute_y_lod",
  function(model, lod, n.mi = 10, eps = 0.01, maxit = 5, family, nfolds,
    foldid, parallel, cluster, verbose) {
    standardGeneric("impute_y_lod")})
setGeneric("compute_error",
  function(model, y, X) standardGeneric("compute_error"))
setGeneric("scale_locs", function(model, locs) standardGeneric("scale_locs"))
setGeneric("transform_covariance_parameters",
  function(model) standardGeneric("transform_covariance_parameters"))
setGeneric("check_input",
  function(model, Y, X, locs, Y.names, X.names, center.y, impute.y, lod) {
    standardGeneric("check_input")})
setGeneric("check_input_pred",
  function(model, X, locs) standardGeneric("check_input_pred"))

#' show method for PrestoGP models
#'
#' This method prints a summary of a PrestoGP model and its parameters.
#'
#' @param object The PrestoGP model object
#'
#' @seealso \code{\link{PrestoGPModel-class}}, \code{\link{prestogp_fit}}
#'
#' @export
#'
#' @examples
#' data(soil)
#' soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
#' y <- soil[,4]                   # predict moisture content
#' X <- as.matrix(soil[,5:9])
#' locs <- as.matrix(soil[,1:2])
#'
#' soil.vm <- new("VecchiaModel", n_neighbors = 10)
#' soil.vm <- prestogp_fit(soil.vm, y, X, locs)
#' show(soil.vm)

setMethod(
  "show", "PrestoGPModel",
  function(object) {
    cat("Matern covariance parameters (theta):", "\n")
    print(get_theta(object))
    cat("Regression coefficients (beta):", "\n")
    cur.coef <- get_beta(object)
    for (i in seq_along(cur.coef)) {
      cur.coef[[i]] <- cur.coef[[i]][cur.coef[[i]] != 0]
    }
    print(cur.coef)
    cat("Model type:", is(object)[1], "\n")
    if (object@n_neighbors > 0) {
      cat("Nearest neighbors:", object@n_neighbors, "\n")
    }
    cat("Scaling:", object@scaling, "\n")
    cat("Penalized likelihood:", object@error, "\n")
    cat("MSE:", mean(object@res^2), "\n")
    if (!object@converged) {
      warning("Model fitting did not converge")
    }
  }
)

#' Extract the Matern (theta) parameters from a PrestoGP model.
#'
#' This method extracts a list containing the Matern variance/covariance
#' parameters from a PrestoGP model.
#'
#' @param model The PrestoGP model object
#'
#' @return A list containing the following elements:
#' \describe{
#' \item{sigma:}{A vector containing the estimates of the variance parameter
#' sigma for each outcome.}
#' \item{scale:}{A vector containing the estimates of the the scale
#' parameter(s) for each outcome.}
#' \item{smoothness:}{A vector containing the estimates of the smoothness
#' parameter for each outcome.}
#' \item{nuggets:}{A vector containing the estimates of the nugget variance
#' for each outcome.}
#' \item{correlation:}{A matrix containing the estimates of the correlations
#' between outcomes. Omitted for univariate models.}
#' }
#'
#' @seealso \code{\link{PrestoGPModel-class}}, \code{\link{prestogp_fit}}
#'
#' @references
#' \itemize{
#' \item Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
#' land-use regression for ground-level nitrogen dioxide", The Annals of
#' Applied Statistics (2021) 15(2):688-710.
#' \item Porcu, E., Bevilacqua, M., Schaback, R., and Oates RJ. "The Matérn
#' Model: A Journey Through Statistics, Numerical Analysis and Machine
#' Learning", Statistical Science (2024) 39(3):469-492.
#' }
#'
#' @aliases get_theta
#' @export
#'
#' @examples
#' data(soil)
#' soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
#' y <- soil[,4]                   # predict moisture content
#' X <- as.matrix(soil[,5:9])
#' locs <- as.matrix(soil[,1:2])
#'
#' soil.vm <- new("VecchiaModel", n_neighbors = 10)
#' soil.vm <- prestogp_fit(soil.vm, y, X, locs)
#' get_theta(soil.vm)

setMethod(
  "get_theta", "PrestoGPModel",
  function(model) {
    pseq <- model@param_sequence
    junk <- list()
    junk$sigma <- model@covparams[pseq[1, 1]:pseq[1, 2]]
    junk$scale <- model@covparams[pseq[2, 1]:pseq[2, 2]]
    junk$smoothness <- model@covparams[pseq[3, 1]:pseq[3, 2]]
    junk$nuggets <- model@covparams[pseq[4, 1]:pseq[4, 2]]
    names(junk$sigma) <- names(model@locs_train)
    names(junk$smoothness) <- names(model@locs_train)
    names(junk$nuggets) <- names(model@locs_train)
    if (model@nscale > 1) {
      scale.names <- rep(NA, model@nscale * length(model@locs_train))
      for (i in seq_along(names(model@locs_train))) {
        for (j in seq_len(model@nscale)) {
          scale.names[(i - 1) * model@nscale + j] <-
            paste0(names(model@locs_train)[i], "_", j)
        }
      }
      names(junk$scale) <- scale.names
    } else {
      names(junk$scale) <- names(model@locs_train)
    }
    if (length(model@locs_train) > 1) {
      P <- length(model@locs_train)
      rho <- model@covparams[pseq[5, 1]:pseq[5, 2]]
      rho.mat <- matrix(0, nrow = P, ncol = P)
      rho.mat[upper.tri(rho.mat, diag = FALSE)] <- rho
      rho.mat <- rho.mat + t(rho.mat)
      diag(rho.mat) <- 1
      colnames(rho.mat) <- names(model@locs_train)
      rownames(rho.mat) <- names(model@locs_train)
      junk$correlation <- rho.mat
    }
    junk
  }
)

#' Extract the regression coefficients (beta) for a PrestoGP model.
#'
#' This method extracts a list containing the regression coefficients for
#' a PrestoGP model.
#'
#' @param model The PrestoGP model object
#'
#' @details It is important to note that the intercepts are estimated on the
#' transformed data. They are not useful for making predictions on the
#' original (untransformed) data. Use \code{\link{prestogp_predict}} to
#' make predictions based on a PrestoGP model.
#'
#' @return A list containing the regression coefficients. Each element of the
#' list corresponds to an outcome variable. The final element of the list
#' contains the intercept(s).
#'
#' @seealso \code{\link{PrestoGPModel-class}}, \code{\link{prestogp_fit}},
#' \code{\link{prestogp_predict}}
#'
#' @references
#' \itemize{
#' \item Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
#' land-use regression for ground-level nitrogen dioxide", The Annals of
#' Applied Statistics (2021) 15(2):688-710.
#' }
#'
#' @aliases get_beta
#' @export
#'
#' @examples
#' data(soil)
#' soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
#' y <- soil[,4]                   # predict moisture content
#' X <- as.matrix(soil[,5:9])
#' locs <- as.matrix(soil[,1:2])
#'
#' soil.vm <- new("VecchiaModel", n_neighbors = 10)
#' soil.vm <- prestogp_fit(soil.vm, y, X, locs)
#' get_beta(soil.vm)

setMethod(
  "get_beta", "PrestoGPModel",
  function(model) {
    junk <- list()
    P <- length(model@locs_train)
    length(junk) <- P + 1
    for (i in seq_len(P)) {
      if (i == 1) {
        indx <- 1
        junk[[P + 1]] <- model@beta[1] + model@Y_bar[1]
        names(junk[[P + 1]]) <- "(Intercept)"
      } else {
        indx <- model@X_ndx[i - 1] + 1
        new.names <- c(names(junk[[P + 1]]),
          names(model@locs_train)[i])
        junk[[P + 1]] <- c(junk[[P + 1]], model@beta[indx] + model@Y_bar[i] -
            model@Y_bar[1])
        names(junk[[P + 1]]) <- new.names
      }
      pred.ndx <- (indx + 1):model@X_ndx[i]
      junk[[i]] <- model@beta[pred.ndx]
      names(junk[[i]]) <- colnames(model@X_train)[pred.ndx - 1]
    }
    if (!is.null(names(model@locs_train))) {
      names(junk) <- c(names(model@locs_train), "(Intercept)")
    } else {
      names(junk) <- c("Y", "(Intercept)")
    }
    junk
  }
)

#' Extract the fitted linear model for a PrestoGP model
#'
#' This method return the fitted linear model (of class cv.glmnet or cv.ncvreg)
#' for a PrestoGP model.
#'
#' @param model The PrestoGP model object
#'
#' @details It is important to note that the model is fit to the
#' transformed data. The CV error rate and predicted values of Y will not be
#' correct for the original (untransformed) data. This method should be
#' used primarily for examining the coefficient paths and generating plots.
#'
#' @return The fitted linear model (of class cv.glmnet or cv.ncvreg).
#'
#' @seealso \code{\link{PrestoGPModel-class}}, \code{\link{prestogp_fit}},
#' \code{\link[glmnet]{cv.glmnet}}, \code{\link[ncvreg]{cv.ncvreg}}
#'
#' @references
#' \itemize{
#' \item Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
#' land-use regression for ground-level nitrogen dioxide", The Annals of
#' Applied Statistics (2021) 15(2):688-710.
#' }
#'
#' @aliases get_linear_model
#' @export
#'
#' @examples
#' data(soil)
#' soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
#' y <- soil[,4]                   # predict moisture content
#' X <- as.matrix(soil[,5:9])
#' locs <- as.matrix(soil[,1:2])
#'
#' soil.vm <- new("VecchiaModel", n_neighbors = 10)
#' soil.vm <- prestogp_fit(soil.vm, y, X, locs)
#' get_linear_model(soil.vm)

setMethod(
  "get_linear_model", "PrestoGPModel",
  function(model) {
    model@linear_model
  }
)

#' Extract the number of neighbors for a PrestoGP model.
#'
#' This method extracts the number of nearest neighbors that were used to
#' compute the PrestoGP model.
#'
#' @param model The PrestoGP model object
#'
#' @return The (integer) number of neighbors used to fit the model.
#'
#' @seealso \code{\link{PrestoGPModel-class}}, \code{\link{prestogp_fit}}
#'
#' @references
#' \itemize{
#' \item Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
#' land-use regression for ground-level nitrogen dioxide", The Annals of
#' Applied Statistics (2021) 15(2):688-710.
#' }
#'
#' @aliases get_neighbors
#' @export
#'
#' @examples
#' data(soil)
#' soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
#' y <- soil[,4]                   # predict moisture content
#' X <- as.matrix(soil[,5:9])
#' locs <- as.matrix(soil[,1:2])
#'
#' soil.vm <- new("VecchiaModel", n_neighbors = 10)
#' soil.vm <- prestogp_fit(soil.vm, y, X, locs)
#' get_neighbors(soil.vm)

setMethod(
  "get_neighbors", "PrestoGPModel",
  function(model) {
    model@n_neighbors
  }
)

#' Extract the scaling parameter for a PrestoGP model.
#'
#' This method extracts the scaling parameter for a PrestoGP model.
#'
#' @param model The PrestoGP model object
#'
#' @return A vector containing the scaling parameter. See
#' \code{\link{prestogp_fit}} for an explanation of this parameter.
#'
#' @seealso \code{\link{PrestoGPModel-class}}, \code{\link{prestogp_fit}}
#'
#' @references
#' \itemize{
#' \item Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
#' land-use regression for ground-level nitrogen dioxide", The Annals of
#' Applied Statistics (2021) 15(2):688-710.
#' }
#'
#' @aliases get_scaling
#' @export
#'
#' @examples
#' data(soil)
#' soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
#' y <- soil[,4]                   # predict moisture content
#' X <- as.matrix(soil[,5:9])
#' locs <- as.matrix(soil[,1:2])
#'
#' soil.vm <- new("VecchiaModel", n_neighbors = 10)
#' soil.vm <- prestogp_fit(soil.vm, y, X, locs)
#' get_scaling(soil.vm)

setMethod(
  "get_scaling", "PrestoGPModel",
  function(model) {
    model@scaling
  }
)

#' Did the PrestoGP model converge?
#'
#' This method returns a boolean value that indicates whether or not a
#' PrestoGP model converged.
#'
#' @param model The PrestoGP model object
#'
#' @return A boolean value indicating whether or not the model converged.
#'
#' @seealso \code{\link{PrestoGPModel-class}}, \code{\link{prestogp_fit}}
#'
#' @references
#' \itemize{
#' \item Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
#' land-use regression for ground-level nitrogen dioxide", The Annals of
#' Applied Statistics (2021) 15(2):688-710.
#' }
#'
#' @aliases get_converged
#' @export
#'
#' @examples
#' data(soil)
#' soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
#' y <- soil[,4]                   # predict moisture content
#' X <- as.matrix(soil[,5:9])
#' locs <- as.matrix(soil[,1:2])
#'
#' soil.vm <- new("VecchiaModel", n_neighbors = 10)
#' soil.vm <- prestogp_fit(soil.vm, y, X, locs)
#' get_converged(soil.vm)

setMethod(
  "get_converged", "PrestoGPModel",
  function(model) {
    model@converged
  }
)

#' Extract the penalized log likelihood for a PrestoGP model.
#'
#' This method extracts the (negative) penalized log likelihood for a PrestoGP
#' model.
#'
#' @param model The PrestoGP model object
#'
#' @return The value of the (negative) penalized log likelihood. See References
#' for details.
#'
#' @seealso \code{\link{PrestoGPModel-class}}, \code{\link{prestogp_fit}}
#'
#' @references
#' \itemize{
#' \item Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
#' land-use regression for ground-level nitrogen dioxide", The Annals of
#' Applied Statistics (2021) 15(2):688-710.
#' }
#'
#' @aliases get_pen_loglik
#' @export
#'
#' @examples
#' data(soil)
#' soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
#' y <- soil[,4]                   # predict moisture content
#' X <- as.matrix(soil[,5:9])
#' locs <- as.matrix(soil[,1:2])
#'
#' soil.vm <- new("VecchiaModel", n_neighbors = 10)
#' soil.vm <- prestogp_fit(soil.vm, y, X, locs)
#' get_pen_loglik(soil.vm)

setMethod(
  "get_pen_loglik", "PrestoGPModel",
  function(model) {
    model@error
  }
)

#' Plots the glide path for the coefficients for a PrestoGP model.
#'
#' This method generates a plot showing the coefficients of the model for
#' different values of the tuning parameter. It is a wrapper for
#' \code{\link[glmnet]{plot.glmnet}} or \code{\link[ncvreg]{plot.ncvreg}}.
#'
#' @param model The PrestoGP model object
#'
#' @param ... Additional parameters to \code{\link[glmnet]{plot.glmnet}} or
#' \code{\link[ncvreg]{plot.ncvreg}}.
#'
#' @seealso \code{\link{PrestoGPModel-class}}, \code{\link{prestogp_fit}},
#' \code{\link[glmnet]{plot.glmnet}}, \code{\link[ncvreg]{plot.ncvreg}}
#'
#' @references
#' \itemize{
#' \item Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
#' land-use regression for ground-level nitrogen dioxide", The Annals of
#' Applied Statistics (2021) 15(2):688-710.
#' }
#'
#' @aliases plot_beta
#' @export
#'
#' @examples
#' data(soil)
#' soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
#' y <- soil[,4]                   # predict moisture content
#' X <- as.matrix(soil[,5:9])
#' locs <- as.matrix(soil[,1:2])
#'
#' soil.vm <- new("VecchiaModel", n_neighbors = 10)
#' soil.vm <- prestogp_fit(soil.vm, y, X, locs)
#' plot_beta(soil.vm)

setMethod(
  "plot_beta", "PrestoGPModel",
  function(model, ...) {
    if (model@penalty == "lasso" || model@penalty == "relaxed") {
      plot(model@linear_model$glmnet.fit, ...)
    } else {
      plot(model@linear_model$fit, ...)
    }
  }
)

#' Train a PrestoGP model.
#'
#' This method fits a PrestoGP model given a set of locations and predictor
#' and outcome variables.
#'
#' @param model The PrestoGP model object being fit.
#' @param Y The values of the response variable(s). Should be a matrix or
#' vector for univariate models or a list for multivariate models.
#' @param X The values of the predictor variable(s). Should be a matrix for
#' univariate models or a list for multivariate models.
#' @param locs The values of the locations. Should be a matrix for univariate
#' models or a list for multivariate models.
#' @param Y.names The name(s) for the response variable(s). Should be a vector
#' with length equal to the number of response variables. Defaults to NULL.
#' @param X.names The names for the predictor variables. Should be a vector
#' for univariate models or a list for multivariate models.
#' @param scaling A vector of consecutive positive integers that is used to
#' specify which columns of locs should have the same scaling parameter. For
#' example, in a spatiotemporal model with two spatial measures and a time
#' measure, the value of scaling would be c(1, 1, 2). The length of scaling
#' must match the number of columns of locs. If it is not specified, all
#' columns of locs will have a common scale parameter.
#' @param common_scale Do all columsn of locs have a common scales parameter?
#' See Details for the effects of this parameter. Defaults to TRUE if there
#' is only one scale parameter for each outcome and FALSE otherwise.
#' @param covparams The initial covariance parameters estimate (optional).
#' @param beta.hat The initial beta parameters estimates (optional).
#' @param tol The model is considered converged when error is not less than
#' tol*previous_error (optional). Defaults to 0.999999.
#' @param max_iters Maximum number of iterations for the model fitting
#' procedure. Defaults to 100.
#' @param center.y Should the Y's be mean centered before fitting the model?
#' Defaults to TRUE for gaussian models and FALSE for binomial models.
#' @param impute.y Should missing Y's be imputed? Defaults to FALSE.
#' @param lod Limit of detection value(s). Any Y value less than lod is
#' assumed to be missing when performing missing data imputation. Should be
#' numeric for univariate models and a list for multivariate models, where
#' each element of the list corresponds to an outcome. The ith element of
#' the lod is the limit of detection for observation i. Alternatively, one
#' can specify a single limit of detection that is assumed to be the same for
#' all observations. If not specified, it is assumed that no limit of
#' detection exists. Ignored if impute.y is FALSE.
#' @param quiet If FALSE, the penalized log likelihood and the model MSE
#' will be printed for each iteration of the model fitting procedure. No
#' intermediate output will be printed if TRUE. Defaults to FALSE.
#' @param verbose If TRUE, the estimated theta/beta parameters will be printed
#' for each iteration of the model fitting procedure. Defaults to FALSE.
#' Ignored if quiet is TRUE.
#' @param optim.method Optimization method to be used for the maximum
#' likelihood estimation that is passed to optim. Defaults to "Nelder-Mead".
#' See \code{\link[stats]{optim}}.
#' @param optim.control Control parameter that is passed to optim. See
#' \code{\link[stats]{optim}}.
#' @param penalty The type of penalized regression to be used. Should be one
#' of "lasso", "relaxed", "MCP", or "SCAD". Note that "lasso" and "relaxed"
#' will fit the model using glmnet and "MCP" and "SCAD" will fit the model
#' using ncvreg. See \code{\link[glmnet]{glmnet}} or
#' \code{\link[ncvreg]{ncvreg}}. Defaults to "lasso".
#' @param alpha The elastic net mixing parameter. 'alpha=1' corresponds to
#' lasso (or SCAD/MCP) penalty; 'alpha=0' corresponds to ridge regression.
#' See \code{\link[glmnet]{glmnet}} or \code{\link[ncvreg]{ncvreg}}. Defaults
#' to 1.
#' @param family Family parameter for the glmnet or ncvreg model. Currently
#' only "gaussian" and "binomial" are supported. Defaults to "gaussian". See
#' \code{\link[glmnet]{glmnet}} or \code{\link[ncvreg]{ncvreg}}.
#' @param nfolds Number of cross-validation folds for cv.glmnet. Defaults to
#' 10. See \code{\link[glmnet]{cv.glmnet}}.
#' @param foldid Optional vector of values between 1 and "nfolds" specifying
#' what fold each observation should be assigned to in the cv.glmnet
#' cross-validation. See \code{\link[glmnet]{cv.glmnet}} and
#' \code{\link[ncvreg]{cv.ncvreg}}.
#' @param parallel Should parallel "foreach" be used to speed up the model
#' fitting procedure where possible? Defaults to FALSE. Specifically,
#' parallelization will be used for imputation and fitting the cv.glmnet
#' object. See \code{\link[glmnet]{cv.glmnet}}. Note that this
#' only applies to glmnet models (where penalty="lasso" or
#' penalty="relaxed"). Models using ncvreg (where penalty="MCP" or
#' penalty="SCAD") require a cluster argument for parallelization (see below).
#' @param cluster A cluster for running cv.ncvreg in parallel. See
#' \code{\link[ncvreg]{cv.ncvreg}} and \code{\link[parallel]{makeCluster}}.
#' This must be specified to run cv.ncvreg in parallel. It is ignored for
#' glmnet models.
#' @param adaptive Should adaptive lasso be used? Defaults to FALSE. It is
#' ignored for SCAD and MCP models.
#'
#' @details If common_scale is TRUE, multivariate models will use the Matern
#' cross-covariance function described in Apanasovich et al. (2012). This
#' model can only be used if each outcome has only a single scale parameter.
#' If common_scale is FALSE, each column of locs will be divded by the
#' corresponding element of scaling, and the cross-covariance will be computed
#' under the assumption that all scale parameters are equal to 1. For
#' univariate models, the Vecchia approximation will not be recomputed in
#' each iteration of the model fitting procedure if common_scale is TRUE,
#' but the parameter is otherwise ignored.
#'
#' @return A PrestoGPModel object with slots updated based on the results of
#' the model fitting procedure. See \code{\link{PrestoGPModel-class}} for
#' details.
#'
#' @seealso \code{\link{PrestoGPModel-class}}, \code{\link[glmnet]{glmnet}}
#'
#' @references
#' \itemize{
#' \item Apanasovich, T.V., Genton, M.G. and Sun, Y. "A valid Matérn class of
#' cross-covariance functions for multivariate random fields with any number
#' of components", Journal of the American Statistical Association (2012)
#' 107(497):180-193.
#' \item Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
#' land-use regression for ground-level nitrogen dioxide", The Annals of
#' Applied Statistics (2021) 15(2):688-710.
#' \item Zou, H. "The adaptive lasso and its oracle properties", Journal of
#' the American Statistical Association (2006) 101(476):1418-1429.
#' }
#'
#' @aliases prestogp_fit
#' @export
#'
#' @examples
#' data(soil)
#' soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
#' y <- soil[,4]                   # predict moisture content
#' X <- as.matrix(soil[,5:9])
#' locs <- as.matrix(soil[,1:2])
#'
#' # Vecchia model
#' soil.vm <- new("VecchiaModel", n_neighbors = 10)
#' soil.vm <- prestogp_fit(soil.vm, y, X, locs)
#'
#' # Impute missing y's
#' miss <- sample(1:nrow(soil), 20)
#' y.miss <- y
#' y.miss[miss] <- NA
#' soil.vm2 <- new("VecchiaModel", n_neighbors = 10)
#' soil.vm2 <- prestogp_fit(soil.vm, y, X, locs, impute.y = TRUE)
#'
#' # Impute y's missing due to limit of detection
#' soil.lod <- quantile(y, 0.1)
#' y.lod <- y
#' y.lod[y.lod <= soil.lod] <- NA
#' soil.vm3 <- new("VecchiaModel", n_neighbors = 10)
#' soil.vm3 <- prestogp_fit(soil.vm, y, X, locs, impute.y = TRUE, lod = soil.lod)
#'
#' # Full model
#' soil.fm <- new("FullModel")
#' soil.fm <- prestogp_fit(soil.fm, y, X, locs)
#'
#' # Multivariate model
#' ym <- list()
#' ym[[1]] <- soil[,5]             # predict two nitrogen concentration levels
#' ym[[2]] <- soil[,7]
#' Xm <- list()
#' Xm[[1]] <- Xm[[2]] <- as.matrix(soil[,c(4,6,8,9)])
#' locsm <- list()
#' locsm[[1]] <- locsm[[2]] <- locs
#'
#' soil.mvm <-  new("MultivariateVecchiaModel", n_neighbors = 10)
#' soil.mvm <- prestogp_fit(soil.mvm, ym, Xm, locsm)
#'
#' # Space/elevation model
#' data(soil250, package="geoR")
#' y2 <- soil250[,7]               # predict pH level
#' X2 <- as.matrix(soil250[,c(4:6,8:22)])
#' # columns 1+2 are location coordinates; column 3 is elevation
#' locs2 <- as.matrix(soil250[,1:3])
#'
#' soil.vm2 <- new("VecchiaModel", n_neighbors = 10)
#' # fit separate scale parameters for location and elevation
#' soil.vm2 <- prestogp_fit(soil.vm2, y2, X2, locs2, scaling = c(1, 1, 2))

setMethod(
  "prestogp_fit", "PrestoGPModel",
  function(model, Y, X, locs, Y.names = NULL, X.names = NULL, scaling = NULL,
    common_scale = NULL, covparams = NULL, beta.hat = NULL, tol = 0.999999,
    max_iters = 100, center.y = NULL, impute.y = FALSE, lod = NULL,
    quiet = FALSE, verbose = FALSE, optim.method = "Nelder-Mead",
    optim.control = list(trace = 0, reltol = 1e-3, maxit = 5000),
    penalty = c("lasso", "relaxed", "MCP", "SCAD"), alpha = 1,
    family = c("gaussian", "binomial"),
    nfolds = 10, foldid = NULL, parallel = FALSE, cluster = NULL,
    adaptive = FALSE) {
    penalty <- match.arg(penalty)
    model@penalty <- penalty
    family <- match.arg(family)
    if (is.null(center.y)) {
      if (family == "gaussian") {
        center.y <- TRUE
      } else {
        center.y <- FALSE
      }
    }
    model <- check_input(model, Y, X, locs, Y.names, X.names, center.y,
      impute.y, lod)
    if (is.null(lod)) {
      lodv <- rep(Inf, nrow(model@Y_train))
    } else {
      if (is.numeric(lod)) {
        lod <- list(lod)
      }
      lodv <- NULL
      for (i in seq_along(lod)) {
        if (length(lod[[i]] == 1)) {
          lod[[i]] <- rep(lod[[i]], nrow(model@locs_train[[i]]))
        }
        lodv <- c(lodv, lod[[i]] - model@Y_bar[i])
      }
    }
    if (!is.null(beta.hat)) {
      if (!is.vector(beta.hat) | !is.numeric(beta.hat)) {
        stop("beta.hat parameter must be a numeric vector")
      }
      if (length(beta.hat) != (ncol(model@X_train) +
            length(model@locs_train))) {
        stop("Length of beta.hat must match the number of predictors")
      }
      beta.hat <- as.matrix(beta.hat)
    }
    if (!is.numeric(alpha)) {
      stop("tol must be numeric")
    }
    if (length(alpha) != 1) {
      stop("tol must be a scalar")
    }
    if (alpha < 0 || alpha > 1) {
      stop("alpha must satisfy 0<=alpha<=1")
    }
    if (alpha == 0 && (penalty == "SCAD" || penalty == "MCP")) {
      stop("alpha must be positive for SCAD/MCP penalties")
    }
    model@alpha <- alpha
    if (quiet & verbose) {
      verbose <- FALSE
    }
    if (!is.numeric(tol)) {
      stop("tol must be numeric")
    }
    if (length(tol) != 1) {
      stop("tol must be a scalar")
    }
    if (tol <= 0 | tol > 1) {
      stop("tol must satisfy 0<tol<=1")
    }
    if (is.null(scaling)) {
      scaling <- rep(1, ncol(model@locs_train[[1]]))
      nscale <- 1
    } else {
      if (length(scaling) != ncol(model@locs_train[[1]])) {
        stop("Length of scaling must equal ncol of locs")
      }
      nscale <- length(unique(scaling))
      if (sum(sort(unique(scaling)) == 1:nscale) < nscale) {
        stop("scaling must consist of sequential integers starting at 1")
      }
    }
    if (is.null(common_scale)) {
      if (nscale == 1) {
        common_scale <- TRUE
      } else {
        common_scale <- FALSE
      }
    }
    if (common_scale & nscale > 1) {
      stop("common_scale must be FALSE if there are multiple scale parameters")
    }
    model@scaling <- scaling
    model@nscale <- nscale
    model@common_scale <- common_scale
    if (model@common_scale & (length(model@locs_train) > 1)) {
      model@locs_train <- eliminate_dupes(model@locs_train)$locs
    }
    if (!is.null(covparams)) {
      if (!is.vector(covparams) | !is.numeric(covparams)) {
        stop("covparams must be a numeric vector")
      }
    }
    model <- calc_covparams(model, locs, Y, covparams)
    if (model@n_neighbors < model@min_m) {
      stop(paste("m must be at least ", model@min_m, sep = ""))
    }
    if (model@n_neighbors >= nrow(model@Y_train)) {
      warning("Conditioning set size m chosen to be >=n. Changing to m=n-1")
      model@n_neighbors <- nrow(model@Y_train) - 1
    }

    model <- specify(model)

    if (is.null(beta.hat)) {
      if (!quiet) {
        cat("\n")
      }
      if (sum(!model@Y_obs) > 0 & min(lodv) < Inf) {
        if (!quiet) {
          cat("Imputing missing y's and estimating initial beta...", "\n")
        }
        cur.mi <- lod_reg_mi(model@Y_train, model@X_train, lodv, !model@Y_obs,
          penalty = model@penalty, alpha = model@alpha,
          parallel = parallel, cluster = cluster, foldid = foldid,
          verbose = verbose)
        model@Y_train[!model@Y_obs] <- cur.mi$y.impute
        beta.hat <- as.matrix(cur.mi$coef)
        if (!quiet) {
          cat("Estimation of initial beta complete", "\n")
        }
      } else {
        if (!quiet) {
          cat("Estimating initial beta...", "\n")
        }
        if (model@penalty == "lasso" || model@penalty == "relaxed") {
          beta0.glmnet <- cv.glmnet(model@X_train, model@Y_train,
            parallel = parallel, foldid = foldid, alpha = model@alpha,
            relax = model@penalty == "relaxed")
          beta.hat <- as.matrix(predict(beta0.glmnet,
              type = "coefficients", s = "lambda.min", gamma = "gamma.min"))
        } else {
          beta0.ncvreg <- cv.ncvreg.wrap(model@X_train, model@Y_train,
            cluster = cluster, foldid = foldid, alpha = model@alpha,
            penalty = model@penalty)
          beta.hat <- as.matrix(coef(beta0.ncvreg))
        }
        if (!quiet) {
          cat("Estimation of initial beta complete", "\n")
        }
      }
    }
    Y.hat <- beta.hat[1, 1] + model@X_train %*% beta.hat[-1, ]
    model@beta <- beta.hat
    if (family == "binomial") {
      Y.hat <- exp(Y.hat) / (1 + exp(Y.hat))
    }

    # Begining algorithm (Algorithm 1 from Messier and Katzfuss 2020)
    model@converged <- FALSE
    prev.error <- 1e10
    iter <- 1
    if (!quiet) {
      cat("\n")
    }
    while (!model@converged && (iter <= max_iters)) {
      if (!quiet) {
        cat("Beginning iteration", iter, "\n")
      }
      model <- compute_residuals(model, model@Y_train, Y.hat, family)
      if (!quiet) {
        cat("Estimating theta...", "\n")
      }
      model <- estimate_theta(model, locs, optim.control, optim.method)
      if (!quiet) {
        cat("Estimation of theta complete", "\n")
        if (verbose) {
          print(get_theta(model))
        }
      }
      # transform data to iid
      if (!model@common_scale) {
        model <- specify(model)
      }

      if (sum(!model@Y_obs) > 0 & min(lodv) < Inf) {
        if (!quiet) {
          cat("Imputing missing y's...", "\n")
        }
        model <- impute_y_lod(model, lodv, family = family, nfolds = nfolds,
          foldid = foldid, parallel = parallel, cluster = cluster,
          verbose = verbose)
        if (!quiet) {
          cat("Imputation complete", "\n")
        }
      }

      if (!quiet) {
        cat("Estimating beta...", "\n")
      }
      model <- transform_data(model, model@Y_train, model@X_train)
      model <- estimate_betas(model, family, nfolds, foldid, parallel,
        cluster, model@penalty, adaptive)
      if (!quiet) {
        cat("Estimation of beta complete", "\n")
        if (verbose) {
          print(get_beta(model))
        }
      }
      min.error <- compute_error(model)
      ### Check min-error against the previous error and tolerance
      if (min.error < prev.error * tol) {
        prev.error <- min.error
        model@error <- prev.error
        #        beta.hat <- sparseToDenseBeta(model@linear_model)
        #        model@beta <- beta.hat
        if (model@penalty == "lasso" || model@penalty == "relaxed") {
          beta.hat <- as.matrix(predict(model@linear_model, s = "lambda.min",
              gamma = "gamma.min", type = "coefficients"))
          Y.hat <- as.matrix(predict(model@linear_model, newx = model@X_train,
              s = "lambda.min", gamma = "gamma.min", type = "response"))
        } else {
          beta.hat <- as.matrix(coef(model@linear_model))
          Y.hat <- as.matrix(predict(model@linear_model, model@X_train))
        }
        model@beta <- beta.hat
        if (sum(!model@Y_obs) > 0 & sum(lodv < Inf) == 0) {
          if (!quiet) {
            cat("Imputing missing y's...", "\n")
          }
          model <- impute_y(model)
          if (!quiet) {
            cat("Imputation complete", "\n")
          }
        }
        covparams.iter <- model@covparams
        Vecchia.SCAD.iter <- model@linear_model
      } else {
        model@converged <- TRUE
        model@beta <- beta.hat
        model@covparams <- covparams.iter
        model@linear_model <- Vecchia.SCAD.iter
        model@error <- prev.error
      }
      if (!quiet) {
        cat("Iteration", iter, "complete", "\n")
        cat("Current penalized negative log likelihood:", model@error, "\n")
        cat("Current MSE:", mean(model@res^2), "\n")
      }
      iter <- iter + 1
    }
    if (!model@converged) {
      warning("Model fitting did not converge")
    }
    model
  }
)

#' estimate_betas
#'
#' Estimate the beta coefficients for a model (not called by user)
#'
#' @param model the model to estimate coeffients for
#'
#' @return A model with updated coefficients
#' @noRd
setMethod("estimate_betas", "PrestoGPModel", function(model, family, nfolds,
  foldid, parallel, cluster, penalty, adaptive) {
  if (penalty == "lasso" || penalty == "relaxed") {
    if (adaptive) {
      ridge.model <- cv.glmnet(as.matrix(model@X_tilde),
        as.matrix(model@y_tilde), alpha = 0, family = family,
        nfolds = nfolds, foldid = foldid, parallel = parallel)
      pen.factor <- 1 / abs(as.numeric(coef(ridge.model,
            s = ridge.model$lambda.min))[-1])
    } else {
      pen.factor <- rep(1, ncol(model@X_tilde))
    }
    model@linear_model <- cv.glmnet(as.matrix(model@X_tilde),
      as.matrix(model@y_tilde), alpha = model@alpha, family = family,
      nfolds = nfolds, foldid = foldid, parallel = parallel,
      penalty.factor = pen.factor, relax = penalty == "relaxed")
  } else {
    model@linear_model <- cv.ncvreg.wrap(as.matrix(model@X_tilde),
      as.matrix(model@y_tilde), cluster = cluster, foldid = foldid,
      penalty = penalty, alpha = model@alpha, family = family, nfolds = nfolds)
  }
  invisible(model)
})

#' sparseToDenseBeta
#'
#' Convert the sparse beta coefficients matrix from glmnet to a dense matrix
#'
#' @param linear_model the glmnet model
#'
#' @return A dense matrix
#sparseToDenseBeta <- function(linear_model) {
#  coefs <- coef(linear_model)
#  if (!is.list(coefs)) {
#    coefs <- list(coefs)
#  }
#  beta_construct <- matrix(data = 0, nrow = coefs[[1]]@Dim[1], ncol = length(coefs))
# coefs[[1]]@Dim[1]+2s because dgCMatrix is 0 offset, and we want to include intercept
#  for (i in seq_along(coefs)) {
#    for (j in seq_along(coefs[[i]]@i)) {
#      k <- coefs[[i]]@i[j]
# beta_construct[k+1,i] <- coefs[[i]]@x[j]
#      beta_construct[k + 1, i] <- coefs[[i]]@x[j]
#    }
#  }
# show(beta_construct)
#  beta <- matrix(beta_construct, nrow = coefs[[1]]@Dim[1], ncol = length(coefs))
#  beta
#}

#' compute_error
#'
#' Compute the error (log likelihood using the GP log likelihood and penalty from the beta coefficients)
#'
#' @param model the PrestoGP model object
#'
#' @return The total error used in the main optimization loop
#' @noRd
setMethod("compute_error", "PrestoGPModel", function(model) {
  ### Betas
  # beta.iter <- sparseToDenseBeta(model@linear_model)
  # beta.iter <- as.numeric(coef(model@linear_model))
  # beta.iter <- do.call(cbind, coef(model@linear_model))
  if (model@penalty == "lasso" || model@penalty == "relaxed") {
    beta.iter <- as.matrix(predict(model@linear_model, type = "coefficients",
        s = "lambda.min", gamma = "gamma.min"))
  } else {
    beta.iter <- as.matrix(coef(model@linear_model))
  }

  ### Get SCAD penalty values
  # LL.vecchia.beta <- SCAD_Penalty_Loglike(beta.iter,lambda.iter)

  ### Compute log-likelihood
  # error <- model@LL_Vecchia_krig + LL.vecchia.beta[model@lambda_1se_idx[[1]]]
  alpha <- model@alpha
  if (model@penalty == "lasso") {
    error <- model@LL_Vecchia_krig + model@linear_model$lambda.min *
      ((1 - alpha) * sqrt(sum(beta.iter^2)) + alpha * sum(abs(beta.iter)))
  } else if (model@penalty == "relaxed") {
    error <- model@LL_Vecchia_krig + model@linear_model$relaxed$lambda.min *
      model@linear_model$relaxed$gamma.min * ((1 - alpha) *
        sqrt(sum(beta.iter^2)) + alpha * sum(abs(beta.iter)))
  } else if (model@penalty == "SCAD") {
    error <- model@LL_Vecchia_krig + model@linear_model$lambda.min *
      (1 - alpha) * sqrt(sum(beta.iter^2)) + alpha * SCAD_penalty(beta.iter,
      model@linear_model$lambda.min, model@linear_model$fit$gamma)
  } else if (model@penalty == "MCP") {
    error <- model@LL_Vecchia_krig + model@linear_model$lambda.min *
      (1 - alpha) * sqrt(sum(beta.iter^2)) + alpha * MCP_penalty(beta.iter,
      model@linear_model$lambda.min, model@linear_model$fit$gamma)
  }
  # Min error (stopping criterion) is the log-likelihood
  error
})

#' calc_covparams
#'
#' Set initial value of covarariance parameters.
#'
#' @param model The model to set the covariance parameters of
#' @param locs the locations matrix
#' @param Y the dependent variable matrix
#'
#' @return a model with initial covariance parameters
#' @noRd
setMethod("calc_covparams", "PrestoGPModel", function(model, locs, Y, covparams) {
  if (!is.list(locs)) {
    P <- 1
    locs <- list(locs)
    Y <- list(Y)
  } else {
    P <- length(locs)
  }
  pseq <- create_param_sequence(P, model@nscale)
  if (is.null(covparams)) {
    col.vars <- rep(NA, P)
    D.sample.bar <- rep(NA, model@nscale * P)
    for (i in 1:P) {
      col.vars[i] <- var(Y[[i]], na.rm = TRUE)
      # N <- length(Y[[i]])
      # TODO find a better way to compute initial spatial range
      for (j in 1:model@nscale) {
        # d.sample <- sample(1:N, max(2, ceiling(N / 50)), replace = FALSE)
        D.sample <- rdist(locs[[i]][, model@scaling == j])
        D.sample.bar[(i - 1) * model@nscale + j] <- mean(D.sample) / 4
      }
    }
    model@logparams <- create.initial.values.flex(
      c(0.9 * col.vars), # marginal variance
      D.sample.bar, # range
      rep(0.5, P), # smoothness
      c(.1 * col.vars), # nuggets
      rep(0, choose(P, 2)),
      P
    )
  } else {
    if (P == 1) {
      if (length(covparams) != pseq[4, 2]) {
        stop("Incorrect number of parameters in covparams")
      }
    } else {
      if (length(covparams) != pseq[5, 2]) {
        stop("Incorrect number of parameters in covparams")
      }
    }
    init.var <- covparams[pseq[1, 1]:pseq[1, 2]]
    init.range <- covparams[pseq[2, 1]:pseq[2, 2]]
    init.smooth <- covparams[pseq[3, 1]:pseq[3, 2]]
    init.nugget <- covparams[pseq[4, 1]:pseq[4, 2]]
    if (P > 1) {
      init.corr <- covparams[pseq[5, 1]:pseq[5, 2]]
    } else {
      init.corr <- 0
    }
    if (sum(init.var <= 0) > 0) {
      stop("Initial variance estimates must be positive")
    }
    if (sum(init.range <= 0) > 0) {
      stop("Initial range estimates must be positive")
    }
    if (sum(init.nugget <= 0) > 0) {
      stop("Initial nugget estimates must be positive")
    }
    if (sum(init.smooth <= 0) > 0 | sum(init.smooth >= 2.5) > 0) {
      stop("Initial smoothness estimates must be between 0 and 2.5")
    }
    if (sum(init.corr < -1) > 0 | sum(init.corr > 1) > 0) {
      stop("Initial correlation estimates must be between -1 and 1")
    }
    model@logparams <- create.initial.values.flex(
      init.var, init.range,
      init.smooth, init.nugget,
      init.corr, P
    )
  }
  model@param_sequence <- pseq
  model <- transform_covariance_parameters(model)
  invisible(model)
})

#' scale_locs
#'
#' Scale the locations matrix by the covariance parameters
#'
#' @param model The model with locations to scale
#' @param locs the locations matrix
#'
#' @return a matrix with scaled locations
#' @noRd
setMethod("scale_locs", "PrestoGPModel", function(model, locs) {
  locs.out <- locs
  if (!model@common_scale) {
    for (i in seq_along(locs)) {
      for (j in 1:model@nscale) {
        locs.out[[i]][, model@scaling == j] <-
          locs[[i]][, model@scaling == j] /
            model@covparams[model@param_sequence[2, 1] + model@nscale * (i - 1) + j - 1]
      }
    }
  }
  locs.out
})

setMethod("transform_covariance_parameters", "PrestoGPModel", function(model) {
  P <- length(model@locs_train)
  if (P > 1) {
    model@covparams <- c(
      exp(model@logparams[1:model@param_sequence[2, 2]]),
      gtools::inv.logit(
        model@logparams[model@param_sequence[3, 1]:model@param_sequence[3, 2]],
        0, 2.5
      ),
      exp(model@logparams[model@param_sequence[4, 1]:model@param_sequence[4, 2]]),
      tanh(model@logparams[model@param_sequence[5, 1]:model@param_sequence[5, 2]])
    )
  } else {
    model@covparams <- c(
      exp(model@logparams[1:model@param_sequence[2, 2]]),
      gtools::inv.logit(
        model@logparams[model@param_sequence[3, 1]:model@param_sequence[3, 2]],
        0, 2.5
      ),
      exp(model@logparams[model@param_sequence[4, 1]:model@param_sequence[4, 2]]), 1
    )
  }
  invisible(model)
})

#' compute_residuals
#'
#' Compute residuals based on beta parameters
#'
#' @param model The model to compute the residual of
#' @param Y the training dependent variable matrix
#' @param Y.hat the predicted training dependent variable matrix based on
#' beta hat
#'
#' @return a model with computed residuals
#' @noRd
setMethod("compute_residuals", "PrestoGPModel", function(model, Y, Y.hat,
  family) {
  res <- as.double(Y - Y.hat)
  if (family == "gaussian") {
    model@res <- res
  } else {
    model@res <- sign(res) * sqrt(-2 * (Y * log(Y.hat) + (1 - Y) *
          log(1 - Y.hat)))
  }
  model@vecchia_approx$zord <- model@res[model@vecchia_approx$ord]
  invisible(model)
})
