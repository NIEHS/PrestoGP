# Look at https://github.com/cran/gpclib/blob/master/R/Rgpc.R
library(methods)

setOldClass("cv.glmnet")


#' PrestoGPModel
#'
#' @slot covparams numeric.
#' @slot beta numeric.
#' @slot lambda_1se_idx numeric.
#' @slot vecchia_approx list.
#' @slot y_tilde numeric.
#' @slot X_tilde dgeMatrix.
#' @slot res numeric.
#' @slot Vecchia_SCAD_fit cv.ncvreg.
#'
#' @export PrestoGPModel
#' @import GPvecchia Matrix fields ncvreg readxl scoringRules MASS glmnet
#' @importFrom stats optim predict var
#' @importFrom aod wald.test
#' @importFrom dplyr %>%
#'
#' @examples
#' @noRd
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
    logparams = "numeric"
  ) # transformed version of the Matern parameters
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

setGeneric("show_theta", function(object, Y_names) standardGeneric("show_theta"))
setGeneric(
  "prestogp_fit",
  function(model, Y, X, locs, scaling = NULL, apanasovich = FALSE,
           covparams = NULL, beta.hat = NULL, tol = 0.999999, max_iters = 100, verbose = FALSE,
           optim.method = "Nelder-Mead", optim.control = list(trace = 0, reltol = 1e-3, maxit = 5000),
           parallel = FALSE, foldid = NULL) {
    standardGeneric("prestogp_fit")
  }
)
setGeneric(
  "prestogp_predict",
  function(model, X = "matrix", locs = "matrix", m = "numeric", ordering.pred = c("obspred", "general"),
           pred.cond = c("independent", "general"), return.values = c("mean", "meanvar")) {
    standardGeneric("prestogp_predict")
  }
)
setGeneric("calc_covparams", function(model, locs, Y, covparams) standardGeneric("calc_covparams"))
setGeneric("specify", function(model, ...) standardGeneric("specify"))
setGeneric("compute_residuals", function(model, Y, Y.hat) standardGeneric("compute_residuals"))
setGeneric("transform_data", function(model, Y, X) standardGeneric("transform_data"))
setGeneric("estimate_theta", function(model, locs, optim.control, method) standardGeneric("estimate_theta"))
setGeneric("estimate_betas", function(model, parallel, foldid) standardGeneric("estimate_betas"))
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
#' This method fits any PrestoGP model given a matrix of locations, a matrix of independent variables, and a matrix of dependent variables.
#'
#' @param model The model object being fit.
#' @param Y A matrix containing training values for the dependent variable.
#' @param X A matrix containing training values for the independent varaibles.
#' @param locs A matrix containing the training spatial coordinates and times.
#' @param covparams The initial covariance parameters to use (optional).
#' @param beta.hat The initial beta parameters to use (optional).
#' @param tol The model is considered converged when error isn't less than tol*previous_error (optional).
#' @param m The number of neighboring datapoints to condition on in the likelihood function (optional).
#' @param verbose If TRUE, additional information about model fit will be printed.
#'
#' @return An object containing model parameters for spatiotemporal prediction.
#' @export
#'
#' @examples
#'
#' US_NO2_Data <- read.csv("data/US_ST_NO2_Data.csv")
#' lat <- US_NO2_Data$Latitude # Latitude
#' lon <- US_NO2_Data$Longitude # Longitude
#' logNO2 <- log(US_NO2_Data$Y) # ozone data
#' logNO2 <- logNO2[1:1000, ]
#' time <- US_NO2_Data$YearFrac # time in fractional years
#' N <- length(logNO2)
#' locs <- cbind(lon, lat, time) # Coordinates in R3 (x,y,t)
#' locs <- locs[1:1000, ]
#' X.Design <- US_NO2_Data[, 7:145] # Design matrix
#' X.Design <- scale(X.Design)
#' X <- X.Design[1:1000, ]
#' model <- PrestoGPSpatiotemporalModel()
#' model <- prestogp_fit(model, logNO2, X, locs)
#' ...
setMethod(
  "prestogp_fit", "PrestoGPModel",
  function(model, Y, X, locs, scaling = NULL, apanasovich = NULL,
           covparams = NULL, beta.hat = NULL, tol = 0.999999,
           max_iters = 100, verbose = FALSE, optim.method = "Nelder-Mead",
           optim.control = list(trace = 0, reltol = 1e-3, maxit = 5000),
           parallel = FALSE, foldid = NULL) {
    model <- check_input(model, Y, X, locs)
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
        type = "coefficients",
        s = beta0.glmnet$lambda.1se
      ))
    }
    Y.hat <- beta.hat[1, 1] + model@X_train %*% beta.hat[-1, ]

    # Begining algorithm (Algorithm 1 from Messier and Katzfuss 2020)
    model@converged <- FALSE
    prev.error <- 1e10
    iter <- 1
    if (verbose) {
      cat("\n")
    }
    while (!model@converged && (iter < max_iters)) {
      model <- compute_residuals(model, model@Y_train, Y.hat)
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
      model <- estimate_betas(model, parallel)
      min.error <- compute_error(model)
      ### Check min-error against the previous error and tolerance
      if (min.error < prev.error * tol) {
        prev.error <- min.error
        model@error <- prev.error
        beta.hat <- sparseToDenseBeta(model@linear_model)
        model@beta <- beta.hat
        # Y.hat <- as.matrix(predict(model@linear_model,newx = X, s=model@linear_model$lambda[model@lambda_1se_idx]))
        Y.hat <- as.matrix(predict(model@linear_model, newx = model@X_train, s = "lambda.1se"))
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
    invisible(model)
  }
)

#' estimate_betas
#'
#' Estimate the beta coefficients for a model (not called by user)
#'
#' @param model the model to estimate coeffients for
#'
#' @return A model with updated coefficients
setMethod("estimate_betas", "PrestoGPModel", function(model, parallel, foldid) {
  if (ncol(model@Y_train) > 1) {
    model@linear_model <- cv.glmnet(
      as.matrix(model@X_tilde),
      as.matrix(model@y_tilde),
      family = "mgaussian",
      alpha = model@alpha,
      parallel = parallel,
      foldid = foldid
    )
  } else {
    model@linear_model <- cv.glmnet(as.matrix(model@X_tilde), as.matrix(model@y_tilde), alpha = model@alpha, parallel = parallel, foldid = foldid)
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
#' @param Y.hat the predicted training dependent variable matrix based on beta hat
#'
#' @return a model with computed residuals
setMethod("compute_residuals", "PrestoGPModel", function(model, Y, Y.hat) {
  model@res <- as.double(Y - Y.hat)
  model@vecchia_approx$zord <- model@res[model@vecchia_approx$ord]
  invisible(model)
})
