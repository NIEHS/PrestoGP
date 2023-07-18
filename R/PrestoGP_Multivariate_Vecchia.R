#' Title
#'
#' @slot model PrestoGPModel.
#'
#' @export MultivariateSpatiotemporalModel
#'
#' @examples
#' @include PrestoGP_Model.R
#' @noRd
MultivariateSpatiotemporalModel <- setClass("MultivariateSpatiotemporalModel",
                                    contains = "PrestoGPModel",
                                    slots = c(
                                      model = "PrestoGPModel",
                                      param_sequence = "matrix",
                                      logparams = "numeric"
                                    ))

validityMultivariateSpatiotemporalModel<-function(object){
  TRUE
}
setValidity("MultivariateSpatiotemporalModel",
            validityMultivariateSpatiotemporalModel)

setMethod("initialize", "MultivariateSpatiotemporalModel", function(.Object, ...) {
  .Object@n_neighbors <- 0#25
  .Object@min_m <- 0
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
#' model <- SpatiotemporalModel()
#' model <- prestogp_fit(model, logNO2, X, locs)
#' prediction <- prestogp_predict(model, X.test, locs.test)
#' Vec.mean <- prediction[[1]]
#' Vec.sds <- prediction[[2]]
setMethod("prestogp_predict", "MultivariateSpatiotemporalModel", function(model, X, locs, m=NULL) {
  #validate parameters
  if(!is.matrix(X)){
    stop("X parameter must be a matrix.")
  }
  if(!is.matrix(locs)){
    stop("The locs parameter must be a matrix.")
  }
  if(nrow(X) != nrow(locs)){
    stop("The number of locations must match the number of X observations.")
  }
  if(ncol(locs) != 3){
    stop("The locs parameter must have 3 columns.")
  }
  if(is.null(m)){ #m defaults to the value used for training
    m <- model@m
  }
  stopifnot((m > 0)) #FIXME m is not required by full model

  # Vecchia prediction at new locations
  Vecchia.Pred <- predict(model@Vecchia_SCAD_fit, X = X, which = model@lambda_1se_idx)
  # Vecchia trend prediction at observed data
  Vecchia.hat <- predict(model@Vecchia_SCAD_fit, X = model@X_train, which = model@lambda_1se_idx)

  # Test set prediction
  res = model@Y_train - Vecchia.hat

  locs.train.scaled = scale_locs(model, model@locs_train)
  locs.scaled = scale_locs(model, locs)
  vec.approx.test = vecchia_specify(locs.train.scaled, m, locs.pred=locs.scaled)

  ## carry out prediction
  pred = vecchia_prediction(res, vec.approx.test, c(model@covparams[1], 1, 0.5), model@covparams[4])

  #prediction function can return both mean and sds
  # returns a list with elements mu.pred,mu.obs,var.pred,var.obs,V.ord
  Vec.mean = pred$mu.pred + Vecchia.Pred #residual + mean trend
  #option to include or exclude theta below
  Vec.sds = sqrt(pred$var.pred + model@covparams[4]) #standard deviation

  return(list("means" = Vec.mean, "standard deviations" = Vec.sds))
})

setMethod("calc_covparams", "MultivariateSpatiotemporalModel", function(model, locs, Y) {
  #cor.matrix <- cor(Y)
  #P <- ncol(Y)
  #col.vars <- apply(Y, 2, var)
  #N <- length(Y)
    P <- length(Y)
    col.vars <- rep(NA, P)
    D.sample.bar <- rep(NA, 2*P)
    for (i in 1:P) {
        col.vars[i] <- var(Y[[i]])
        N <- length(Y[[i]])
        #TODO find a better way to compute initial spatial range
        d.sample <- sample(1:N,max(2, ceiling(N/50)),replace = FALSE)
        D.sample = rdist(locs[[i]][d.sample,1:2])
        D.sample.bar[2*i-1] <- mean(D.sample)/4
        d.sample <- sample(1:N,max(2, ceiling(N/50)),replace = FALSE)
        D.sample = rdist(locs[[i]][d.sample,1:2])
        D.sample.bar[2*i] <- mean(D.sample)/4
    }
  model@logparams <- create.initial.values.flex(c(0.9*col.vars), #marginal variance
                                               D.sample.bar, #range
                                               rep(0.5,P), #smoothness
                                               c(.1*col.vars), #nuggets
                                               rep(0, choose(P,2)),
                                               P)
  model@param_sequence <- create.param.sequence(P, 2)
  model <- transform_covariance_parameters(model)
  invisible(model)
})

setMethod("specify", "MultivariateSpatiotemporalModel", function(model, locs, m) {
  locs.scaled = scale_locs(model, locs)
  model@vecchia_approx=vecchia_Mspecify(locs.scaled,m)
  olocs.scaled <- model@vecchia_approx$locsord
  for (i in 1:length(locs)) {
      olocs.scaled[model@vecchia_approx$ondx==i,1:2] <-
          olocs.scaled[model@vecchia_approx$ondx==i,1:2] *
          model@covparams[model@param_sequence[2,1]+2*i-2]
      olocs.scaled[model@vecchia_approx$ondx==i,3] <-
          olocs.scaled[model@vecchia_approx$ondx==i,3] *
          model@covparams[model@param_sequence[2,1]+2*i-1]
  }
  model@vecchia_approx$locsord <- olocs.scaled
  invisible(model)
})

setMethod("scale_locs", "MultivariateSpatiotemporalModel", function(model, locs) {
    locs.out <- locs
    for (i in 1:length(locs)) {
        locs.out[[i]][,1:2] <- locs[[i]][,1:2] /
            model@covparams[model@param_sequence[2,1]+2*i-2]
        locs.out[[i]][,3] <- locs[[i]][,3] /
            model@covparams[model@param_sequence[2,1]+2*i-1]
    }
    return(locs.out)
})

setMethod("compute_residuals", "MultivariateSpatiotemporalModel", function(model, Y, Y.hat) {
  model@res = as.double(Y-Y.hat)
  model@vecchia_approx$zord = model@res[model@vecchia_approx$ord]
  invisible(model)
})

setMethod("estimate_theta", "MultivariateSpatiotemporalModel", function(model, locs, optim.control, method) {
  P <- length(locs)
#  locs_list <- list()
#  y <- list()
#  for(i in 1:P){
#    locs_list[[i]] <- locs
#    y[[i]] <- model@Y_train[,i]
#  }
#  show(model@covparams)
#  show(model@param_sequence)
#  show(model@logparams)
  vecchia.result<- optim(par = model@logparams,
                  fn = mvnegloglik_ST,
                  vecchia.approx=model@vecchia_approx,
                  y = model@res,
                  P = P,
                  param.seq = model@param_sequence,
                  method = method,
                  control=optim.control)

  model@LL_Vecchia_krig <- vecchia.result$value
  model@logparams <- vecchia.result$par
  model <- transform_covariance_parameters(model)
  invisible(model)
})

setMethod("transform_covariance_parameters", "MultivariateSpatiotemporalModel", function(model) {
  P <- length(model@Y_train)
  if(P > 1){
      model@covparams <- c(exp(model@logparams[1:model@param_sequence[2,2]]),
                           gtools::inv.logit(model@logparams[model@param_sequence[3,1]:
                                                   model@param_sequence[3,2]],
                                         0,2.5),
                           exp(model@logparams[model@param_sequence[4,1]:
                                               model@param_sequence[4,2]]),
                           tanh(model@logparams[model@param_sequence[5,1]:
                                                model@param_sequence[5,2]]))
  } else {
      model@covparams <- c(exp(model@logparams[1:model@param_sequence[2,2]]),
                           gtools::inv.logit(model@logparams[model@param_sequence[3,1]:
                                                   model@param_sequence[3,2]],
                                         0,2.5),
                           exp(model@logparams[model@param_sequence[4,1]:
                                               model@param_sequence[4,2]]),1)
  }
  invisible(model)
})

setMethod("transform_data", "MultivariateSpatiotemporalModel", function(model, Y, X) {
    vecchia.approx <- model@vecchia_approx
    params <- model@covparams
    param.seq <- model@param_sequence
    locs.scaled <- vecchia.approx$locsord
    for (i in 1:vecchia.approx$P) {
        locs.scaled[vecchia.approx$ondx==i,1:2] <-
            locs.scaled[vecchia.approx$ondx==i,1:2] /
            params[param.seq[2,1]+2*i-2]
        locs.scaled[vecchia.approx$ondx==i,3] <-
            locs.scaled[vecchia.approx$ondx==i,3] /
            params[param.seq[2,1]+2*i-1]
    }
    vecchia.approx$locsord <- locs.scaled

   transformed.data=transform_miid(cbind(Y,as.matrix(X)),
                                   vecchia.approx = vecchia.approx,
                                   c(params[1:param.seq[1,2]],
                                   rep(1, param.seq[2,2]-param.seq[2,1]+1),
                                   params[param.seq[3,1]:param.seq[5,2]]))
   xcols <- ncol(model@X_train)
   ycols <- ncol(model@Y_train)
   tcols <- ncol(transformed.data)
   model@y_tilde <- Matrix(transformed.data[,1:ncol(model@Y_train)])
   model@X_tilde <- Matrix(transformed.data[,(ncol(model@Y_train)+1):ncol(transformed.data)], sparse=FALSE)
   invisible(model)
})

setMethod("theta_names", "MultivariateSpatiotemporalModel", function(model) {
  c("Marginal Variance", "Range", "Smoothness", "Nugget")
})
