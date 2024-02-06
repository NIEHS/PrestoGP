setOldClass("cv.glmnet")

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
#' @slot beta A matrix containing the regression coefficients.
#' @slot lambda_1se_idx Stores the index of the optimal tuning parameter for
#' the glmnet model. See \code{\link[glmnet]{glmnet}}.
#' @slot vecchia_approx The output of the Vecchia specify function. See
#' \code{\link[GPvecchia]{vecchia_specify}} and \code{\link{vecchia_Mspecify}}.
#' @slot X_train A matrix containing the original predictors. This will be a
#' "super matrix" for multivariate models. See
#' \code{\link[psych]{superMatrix}}.
#' @slot Y_train A column matrix containing the original response values.
#' @slot X_tilde The matrix of transformed predictors.
#' @slot y_tilde The column matrix containing the transformed response values.
#' @slot res numeric.
#' @slot locs_train A list containing the location coordinates. Each element
#' of the list corresponds to a different outcome. (The list will have length
#' 1 for univariate models.)
#' @slot res A numeric vector of the residuals.
#' @slot linear_model The glmnet model. See \code{\link[glmnet]{glmnet}} and
#' \code{\link[glmnet]{cv.glmnet}}.
#' @slot converged Did the model fitting process converge (boolean)?
#' @slot LL_Vecchia_krig The value of the negative log likelihood function
#' after optimization.
#' @slot error Penalized model error. See References for details.
#' @slot n_eighbors Number of neighbors to condition on for the Vecchia
#' approximation. Ignored for full models.
#' @slot min_m Minimum permissible number of neighbors.
#' @slot alpha Parameter alpha for glmnet. See \code{\link[glmnet]{glmnet}}.
#' @slot scaling The indices of the scale parameters. See
#' \code{link{prestogp_fit}}.
#' @slot nscale The number of scale parameters in the model.
#' @slot apanasovich Should the Apanasovich covariance model be used? See
#' References.
#' @slot param_sequence Records the indices of the various Matern parameters.
#' See \code{\link{create.param.sequence}}.
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
    lambda_1se_idx = "numeric", # "list", #the index of the best model
    vecchia_approx = "list", # the output of vecchia_specify
    y_tilde = "dgeMatrix", # iid transformed dependent variables matrix
    X_tilde = "dgeMatrix", # iid transformed independent variables matrix
    res = "numeric", # residuals
    linear_model = "cv.glmnet", # the linear model
    X_train = "matrix", # the original independent variable matrix
    Y_train = "matrix", # the original dependent variable matrix
    locs_train = "list", # the location / temporal matrix
    converged = "logical", # a logical variable that is true if the model fitting process has converged
    LL_Vecchia_krig = "numeric", # the value of the negative log likelihood function after optimization
    error = "numeric", # negative log likelihood + SCAD penalty likelihood
    n_neighbors = "numeric", # the number of neighbors to condition on for the Vecchia approximation
    min_m = "numeric", # the minimum m required by the specific model type (full vs Vecchia)
    alpha = "numeric", # the alpha ratio of ridge to lasso penalty
    scaling = "numeric", # the indices of the scale parameters,
    nscale = "numeric", # the number of scale parameters
    apanasovich = "logical", # should the Apanasovich model be used
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
  .Object@linear_model <- structure(list(), class = "cv.glmnet")
  .Object@alpha <- 1 # 0.5
  .Object <- callNextMethod()
  validObject(.Object)
  .Object
})

setGeneric("show_theta", function(object, Y_names)
    standardGeneric("show_theta"))
setGeneric(
  "prestogp_fit",
  function(model, Y, X, locs, scaling = NULL, apanasovich = FALSE,
    covparams = NULL, beta.hat = NULL, tol = 0.999999, max_iters = 100,
    verbose = FALSE, optim.method = "Nelder-Mead",
    optim.control = list(trace = 0, reltol = 1e-3, maxit = 5000),
    family = c("gaussian", "binomial"), nfolds = 10, foldid = NULL,
    parallel = FALSE) {
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
  function(model, X = "matrix", locs = "matrix", m = "numeric", ordering.pred = c("obspred", "general"),
    pred.cond = c("independent", "general"), return.values = c("mean", "meanvar")) {
    standardGeneric("prestogp_predict")
  }
)
setGeneric("calc_covparams", function(model, locs, Y, covparams) standardGeneric("calc_covparams"))
setGeneric("specify", function(model, ...) standardGeneric("specify"))
setGeneric("compute_residuals", function(model, Y, Y.hat, family) standardGeneric("compute_residuals"))
setGeneric("transform_data", function(model, Y, X) standardGeneric("transform_data"))
setGeneric("estimate_theta", function(model, locs, optim.control, method) standardGeneric("estimate_theta"))
setGeneric("estimate_betas", function(model, family, nfolds, foldid, parallel) standardGeneric("estimate_betas"))
setGeneric("compute_error", function(model, y, X) standardGeneric("compute_error"))
setGeneric("scale_locs", function(model, locs) standardGeneric("scale_locs"))
setGeneric("theta_names", function(model) standardGeneric("theta_names"))
setGeneric("transform_covariance_parameters", function(model) standardGeneric("transform_covariance_parameters"))
setGeneric("check_input", function(model, Y, X, locs) standardGeneric("check_input"))
setGeneric("check_input_pred", function(model, X, locs) standardGeneric("check_input_pred"))

#' show
#'
#' Print a summary of the model and its parameters
#'
#' @param object the PrestoGP model object
setMethod(
  "show", "PrestoGPModel", # TODO consider exporting this
  function(object) {
    cat("PrestoGP Model\n") # TODO print out type of model
    cat("Negative Log-Likelihood: ", object@error, "\n")
    cat("Covariance Parameters:\n")
    Y_names <- colnames(object@Y_train)
    if (is.null(Y_names)) {
      Y_names <- unlist(lapply(seq_len(ncol(object@Y_train)), function(x) {
        paste("Outcome", x)
      }))
    }
    show_theta(object, Y_names)

    y_hat <- matrix(predict(object@linear_model, newx = object@X_train), nrow = nrow(object@X_train), ncol = ncol(object@Y_train))
    mse <- crossprod((object@Y_train - y_hat)) / (nrow(object@Y_train) - colSums(object@beta))
    cat("\nTraining MSE: ", diag(mse), "\n")
    X <- cbind(1, object@X_train)
    covm <- MASS::ginv(t(X) %*% X)
    cat("Non-zero Beta Parameters:\n")
    # TODO compare to zero within a tolerance
    # nnz_betas <- lapply(object@beta, 2, function(x){which(x != 0.0)})
    nnz_betas <- list()
    for (col in seq_len(ncol(object@Y_train))) {
      nnz_betas <- append(nnz_betas, list(which(object@beta[, col] != 0.0)))
    }
    X_names <- colnames(object@X_train)
    if (is.null(X_names)) {
      X_names <- unlist(lapply(seq_len(ncol(object@X_train)), function(x) {
        paste("Ind. Variable", x)
      }))
    }
    X_names <- append("Intercept", X_names)
    for (i in seq_len(ncol(object@Y_train))) {
      cat(Y_names[i], " Parameters:\n")
      beta_summary <- data.frame(matrix(ncol = 4, nrow = 0, dimnames = list(NULL, c("Parameter", "Estimate", "Standard Error", "Walds P-value"))))
      # for(nnz in nnz_betas[i]){
      for (j in seq_along(nnz_betas[[i]])) {
        nnz <- nnz_betas[[i]][[j]]
        walds <- wald.test(covm * mse[i, i], object@beta[, i], Terms = nnz)
        std_err <- sqrt(diag(covm) * mse[i, i])
        walds_p <- walds$result$chi2[3]
        beta_summary[nrow(beta_summary) + 1, ] <- list(X_names[nnz], object@beta[nnz, i], std_err[nnz], walds_p)
      }
      print(beta_summary, row.names = FALSE)
      cat("\n")
    }
    invisible(object)
  }
)

#' show_theta
#'
#' Print the covariance parameters in a table
#'
#' @param object the PrestoGP model object
#' @param Y_names the names of the different outcome variables (may just be numbers if not provided in training input)
setMethod(
  "show_theta", "PrestoGPModel",
  function(object, Y_names) {
    theta_name_arr <- theta_names(object)
    theta_summary <- data.frame(matrix(ncol = ncol(object@Y_train) + 1, nrow = length(theta_name_arr), dimnames = list(NULL, c("Parameter", Y_names))))
    for (i in seq_along(theta_name_arr)) {
      theta_row <- object@covparams[((i - 1) * ncol(object@Y_train) + 1):(i * ncol(object@Y_train))]
      for (j in seq_len(ncol(object@Y_train))) {
        theta_summary[i, j + 1] <- theta_row[j]
      }
    }
    for (j in seq_along(theta_name_arr)) {
      theta_summary[j, 1] <- theta_name_arr[j]
    }
    print(theta_summary, row.names = FALSE)
    # TODO show Rho matrix if there are 2 or more outcomes
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
#' @param scaling A vector of consecutive positive integers that is used to
#' specify which columns of locs should have the same scaling parameter. For
#' example, in a spatiotemporal model with two spatial measures and a time
#' measure, the value of scaling would be c(1, 1, 2). The length of scaling
#' must match the number of columns of locs. If it is not specified, all
#' columns of locs will have a common scale parameter.
#' @param apanasovich Should the multivariate Matern model described in
#' Apanasovich et al. (2012) be used? Defaults to TRUE if there is only one
#' scale parameter for each outcome and FALSE otherwise.
#' @param covparams The initial covariance parameters estimate (optional).
#' @param beta.hat The initial beta parameters estimates (optional).
#' @param tol The model is considered converged when error is not less than
#' tol*previous_error (optional). Defaults to 0.999999.
#' @param max_iters Maximum number of iterations for the model fitting
#' procedure. Defaults to 100.
#' @param verbose If TRUE, additional information about model fit will be
#' printed. Defaults to FALSE.
#' @param optim.method Optimization method to be used for the maximum
#' likelihood estimation that is passed to optim. Defaults to "Nelder-Mead".
#' See \code{\link[stats]{optim}}.
#' @param optim.control Control parameter that is passed to optim. See
#' \code{\link[stats]{optim}}.
#' @param family Family parameter for the glmnet model. Currently only
#' "gaussian" and "binomial" are supported. Defaults to "gaussian". See
#' \code{\link[glmnet]{glmnet}}.
#' @param nfolds Number of cross-validation folds for cv.glmnet. Defaults to
#' 10. See \code{\link[glmnet]{cv.glmnet}}.
#' @param foldid Optional vector of values between 1 and "nfolds" specifying
#' what fold each observation should be assigned to in the cv.glmnet
#' cross-validation. See \code{\link[glmnet]{cv.glmnet}}.
#' @param parallel Should cv.glmnet use parallel "foreach" to fit each fold?
#' Defaults to FALSE. See \code{\link[glmnet]{cv.glmnet}}.
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
  function(model, Y, X, locs, scaling = NULL, apanasovich = NULL,
    covparams = NULL, beta.hat = NULL, tol = 0.999999,
    max_iters = 100, verbose = FALSE, optim.method = "Nelder-Mead",
    optim.control = list(trace = 0, reltol = 1e-3, maxit = 5000),
    family = c("gaussian", "binomial"),
    nfolds = 10, foldid = NULL, parallel = FALSE) {
    model <- check_input(model, Y, X, locs)
    family <- match.arg(family)
    if (!is.null(beta.hat)) {
      if (!is.vector(beta.hat) | !is.numeric(beta.hat)) {
        stop("beta.hat parameter must be a numeric vector")
      }
      if (length(beta.hat) != (ncol(model@X_train) + 1)) {
        stop("Length of beta.hat must match the number of predictors")
      }
      beta.hat <- as.matrix(beta.hat)
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
    if (is.null(apanasovich)) {
      if (nscale == 1) {
        apanasovich <- TRUE
      } else {
        apanasovich <- FALSE
      }
    }
    if (apanasovich & nscale > 1) {
      stop("Apanasovich models require a common scale parameter")
    }
    model@scaling <- scaling
    model@nscale <- nscale
    model@apanasovich <- apanasovich
    if (model@apanasovich & (length(model@locs_train) > 1)) {
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
      beta0.glmnet <- cv.glmnet(model@X_train, model@Y_train,
        parallel = parallel,
        foldid = foldid
      )
      beta.hat <- as.matrix(predict(beta0.glmnet,
          type = "coefficients", s = beta0.glmnet$lambda.1se))
    }
    Y.hat <- beta.hat[1, 1] + model@X_train %*% beta.hat[-1, ]
    if (family == "binomial") {
      Y.hat <- exp(Y.hat) / (1 + exp(Y.hat))
    }

    # Begining algorithm (Algorithm 1 from Messier and Katzfuss 2020)
    model@converged <- FALSE
    prev.error <- 1e10
    iter <- 1
    if (verbose) {
      cat("\n")
    }
    while (!model@converged && (iter < max_iters)) {
      model <- compute_residuals(model, model@Y_train, Y.hat, family)
      res_matrix <- matrix(model@res, nrow = nrow(model@Y_train), ncol = ncol(model@Y_train))
      if (verbose) {
        cat("MSE: ", colMeans(res_matrix^2), "\n")
      }
      model <- estimate_theta(model, locs, optim.control, optim.method)
      # transform data to iid
      if (!model@apanasovich) {
        model <- specify(model)
      }
      model <- transform_data(model, model@Y_train, model@X_train)
      model <- estimate_betas(model, family, nfolds, foldid, parallel)
      min.error <- compute_error(model)
      ### Check min-error against the previous error and tolerance
      if (min.error < prev.error * tol) {
        prev.error <- min.error
        model@error <- prev.error
        beta.hat <- sparseToDenseBeta(model@linear_model)
        model@beta <- beta.hat
        Y.hat <- as.matrix(
          predict(
            model@linear_model,
            newx = model@X_train,
            s = "lambda.1se",
            type = "response")
        )
        covparams.iter <- model@covparams
        Vecchia.SCAD.iter <- model@linear_model
      } else {
        model@converged <- TRUE
        model@beta <- beta.hat
        model@covparams <- covparams.iter
        model@linear_model <- Vecchia.SCAD.iter
        model@error <- prev.error
      }
      if (verbose) {
        cat("\nIteration: ", iter, "\n")
        show(model)
      }
      iter <- iter + 1
    }
    return(model)
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
  foldid, parallel) {
  if (ncol(model@Y_train) > 1) {
    model@linear_model <- cv.glmnet(
      as.matrix(model@X_tilde),
      as.matrix(model@y_tilde),
      family = "mgaussian",
      alpha = model@alpha,
      nfolds = nfolds,
      foldid = foldid,
      parallel = parallel
    )
  } else {
      model@linear_model <- cv.glmnet(as.matrix(model@X_tilde),
        as.matrix(model@y_tilde), alpha = model@alpha, family = family,
        fnolds = nfolds, foldid = foldid, parallel = parallel)
  }
  idmin <- which(model@linear_model$lambda == model@linear_model$lambda.min)
  semin <- model@linear_model$cvm[idmin] + model@linear_model$cvsd[idmin]
  lambda_1se <- max(model@linear_model$lambda[model@linear_model$cvm <= semin])
  model@lambda_1se_idx <- which(model@linear_model$lambda == lambda_1se)
  invisible(model)
})

#' sparseToDenseBeta
#'
#' Convert the sparse beta coefficients matrix from glmnet to a dense matrix
#'
#' @param linear_model the glmnet model
#'
#' @return A dense matrix
#' @noRd
sparseToDenseBeta <- function(linear_model) {
  coefs <- coef(linear_model)
  if (!is.list(coefs)) {
    coefs <- list(coefs)
  }
  beta_construct <- matrix(data = 0, nrow = coefs[[1]]@Dim[1], ncol = length(coefs))
  # coefs[[1]]@Dim[1]+2s because dgCMatrix is 0 offset, and we want to include intercept
  for (i in seq_along(coefs)) {
    for (j in seq_along(coefs[[i]]@i)) {
      k <- coefs[[i]]@i[j]
      # beta_construct[k+1,i] <- coefs[[i]]@x[j]
      beta_construct[k + 1, i] <- coefs[[i]]@x[j]
    }
  }
  # show(beta_construct)
  beta <- matrix(beta_construct, nrow = coefs[[1]]@Dim[1], ncol = length(coefs))
  beta
}

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
  beta.iter <- sparseToDenseBeta(model@linear_model)
  # beta.iter <- as.numeric(coef(model@linear_model))
  # beta.iter <- do.call(cbind, coef(model@linear_model))

  ### Lambdas
  lambda.iter <- model@linear_model$lambda[model@lambda_1se_idx]

  ### Get SCAD penalty values
  # LL.vecchia.beta <- SCAD_Penalty_Loglike(beta.iter,lambda.iter)

  ### Compute log-likelihood
  # error <- model@LL_Vecchia_krig + LL.vecchia.beta[model@lambda_1se_idx[[1]]]
  error <- model@LL_Vecchia_krig + glmnet_penalty(beta.iter, lambda.iter, model@alpha)

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
  pseq <- create.param.sequence(P, model@nscale)
  if (is.null(covparams)) {
    col.vars <- rep(NA, P)
    D.sample.bar <- rep(NA, model@nscale * P)
    for (i in 1:P) {
      col.vars[i] <- var(Y[[i]])
      N <- length(Y[[i]])
      # TODO find a better way to compute initial spatial range
      for (j in 1:model@nscale) {
        d.sample <- sample(1:N, max(2, ceiling(N / 50)), replace = FALSE)
        D.sample <- rdist(locs[[i]][d.sample, model@scaling == j])
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
  if (model@apanasovich) {
    return(locs)
  } else {
    locs.out <- locs
    for (i in seq_along(locs)) {
      for (j in 1:model@nscale) {
        locs.out[[i]][, model@scaling == j] <-
          locs[[i]][, model@scaling == j] /
            model@covparams[model@param_sequence[2, 1] + model@nscale * (i - 1) + j - 1]
      }
    }
    return(locs.out)
  }
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
