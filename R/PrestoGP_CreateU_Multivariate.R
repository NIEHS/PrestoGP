#' Extract specific Matern parameters from a parameter sequence
#'
#' This function is used to obtain specific Matern parameters (e.g.,
#' range or smoothness) from the covparams slot of a PrestoGPModel object.
#'
#' @param P Number of outcome variables
#' @param ns Number of scale parameters
#'
#' @details This function is intended for advanced users who want to specify
#' the input Matern parameters for functions such as
#' \code{\link{vecchia_Mlikelihood}} or \code{\link{createUMultivariate}}.
#' To extract the Matern parameters from a fitted PrestoGP model, it is
#' strongly recommended to use \code{link{get_theta}} instead.
#'
#' @return A matrix with five rows and two columns as described below:
#' \describe{
#' \item{Row 1:}{Starting and ending indices for the sigma parameter(s)}
#' \item{Row 2:}{Starting and ending indices for the scale parameter(s)}
#' \item{Row 3:}{Starting and ending indices for the smoothness parameter(s)}
#' \item{Row 4:}{Starting and ending indices for the nugget(s)}
#' \item{Row 5:}{Starting and ending indices for the correlation parameter(s)}
#' }
#'
#' @seealso \code{\link{PrestoGPModel-class}}
#'
#' @references
#' \itemize{
#' \item Apanasovich, T.V., Genton, M.G. and Sun, Y. "A valid Matérn class of
#' cross-covariance functions for multivariate random fields with any number
#' of components", Journal of the American Statistical Association (2012)
#' 107(497):180-193.
#' \item Genton, M.G. "Classes of kernels for machine learning: a statistics
#' perspective", The Journal of Machine Learning Research (2001) 2:299-312.
#' }
#'
#' @export
#' @examples
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
#'
#' pseq <- create_param_sequence(1, 2)
#' soil2.params <- soil.vm2@covparams
#' # sigma
#' soil2.params[pseq[1,1]:pseq[1,2]]
#' # scale parameters
#' soil2.params[pseq[2,1]:pseq[2,2]]
#' # smoothness parameter
#' soil2.params[pseq[3,1]:pseq[3,2]]
#' # nugget
#' soil2.params[pseq[4,1]:pseq[4,2]]
#'
#' # Multivariate model
#' ym <- list()
#' ym[[1]] <- soil250[,4] # predict sand/silt portion of the sample
#' ym[[2]] <- soil250[,5]
#' ym[[3]] <- soil250[,6]
#' Xm <- list()
#' Xm[[1]] <- Xm[[2]] <- Xm[[3]] <- as.matrix(soil250[,7:22])
#' locsm <- list()
#' locsm[[1]] <- locsm[[2]] <- locsm[[3]] <- as.matrix(soil250[,1:3])
#'
#' soil.mvm <-  new("MultivariateVecchiaModel", n_neighbors = 10)
#' soil.mvm <- prestogp_fit(soil.mvm, ym, Xm, locsm)
#'
#' pseq <- create_param_sequence(3, 2)
#' soil.params <- soil.mvm@covparams
#' # sigmas
#' soil.params[pseq[1,1]:pseq[1,2]]
#' # scale parameters
#' scale.seq <- pseq[2,1]:pseq[2,2]
#' # scale parameter for location, outcome 1
#' soil.params[scale.seq[1]]
#' # scale parameter for elevation, outcome 1
#' soil.params[scale.seq[2]]
#' # scale parameter for location, outcome 2
#' soil.params[scale.seq[3]]
#' # scale parameter for elevation, outcome 2
#' soil.params[scale.seq[4]]
#' # scale parameter for location, outcome 3
#' soil.params[scale.seq[5]]
#' # scale parameter for elevation, outcome 3
#' soil.params[scale.seq[6]]
#' # smoothness parameters
#' soil.params[pseq[3,1]:pseq[3,2]]
#' # nuggets
#' soil.params[pseq[4,1]:pseq[4,2]]
#' # correlation
#' soil.corr <- diag(2) / 2
#' soil.corr[upper.tri(soil.corr)] <- soil.params[pseq[5,1]:pseq[5,2]]
#' soil.corr <- soil.corr + t(soil.corr)
create_param_sequence <- function(P, ns = 1) {
  nk <- choose(P, 2)
  if (nk == 0) {
    nk <- 1 # univariate case
  }

  param.sequence.begin <- c(1, P + 1, seq(P * (ns + 1) + 1, length = 3, by = P))
  param.sequence.end <- c(P, ns * P, P, P, nk) %>% cumsum()
  param.sequence <- cbind(param.sequence.begin, param.sequence.end)

  param.sequence
}

#' Maximum minimum distance ordering
#'
#' Returns the indices of an exact maximum-minimum distance ordering. This is
#' similar to the \code{\link[GPvecchia]{order_maxmin_exact}} function. The
#' main difference is that it allows the user to specify a distance function.
#'
#' @param locs A matrix with one row per location and any number of columns
#' (x, y, time, etc.).
#' @param dist.func Any distance function with a signature of
#' dist(query_location, locations_matrix). Defaults to Euclidean distance.
#'
#' @details For Euclidean distance, this function will return the same
#' results as \code{\link[GPvecchia]{order_maxmin_exact}}, but will be much
#' slower for large data sets. \code{\link[GPvecchia]{order_maxmin_exact}}
#' should be used instead of this function when the distance function is
#' Euclidean.
#'
#' @return A vector of indices giving the ordering. Element \emph{i} of this
#' vector is the index of the \emph{i}th location.
#'
#' @seealso \code{\link[GPvecchia]{order_maxmin_exact}}
#'
#' @references
#' \itemize{
#' \item Katzfuss, M., and Guinness, J. "A general framework for Vecchia
#' approximations of Gaussian processes", Statistical Science (2021)
#' 36(1):124-141.
#' \item Guiness, J. "Permutation methods for sharpening Gaussian process
#' approximations", Technometrics (2018) 60(4):415-429.
#' }
#'
#' @export
#' @examples
#' data(weather)
#' locs <- weather[,3:4]
#' max_min_ordering(locs)
max_min_ordering <- function(locs, dist.func = NULL) {
  if (is.null(dist.func)) {
    dist.func <- fields::rdist
  }
  center <- matrix(colMeans(locs), ncol = ncol(locs))
  # find the point closest to the mean of all points
  dists <- dist.func(center, locs)
  first <- which.min(dists)
  unsolved <- seq_len(nrow(locs))
  unsolved <- unsolved[-first]
  order <- c(first)

  while (length(order) < nrow(locs)) {
    max_min <- 0
    max_min_i <- unsolved[1]
    in_order <- locs[order[seq_along(order)], ]
    dim(in_order) <- c(length(order), ncol(locs))
    for (i in unsolved) {
      loc_i <- locs[i, ]
      dim(loc_i) <- c(1, ncol(locs))
      dists <- dist.func(loc_i, in_order)
      candiate_dist <- min(dists) # min distance from loc(i) to already ordered points
      if (candiate_dist > max_min) {
        max_min <- candiate_dist # loc(i) has a larger minimum distance from already ordered points
        max_min_i <- i
      }
    }
    order <- c(order, max_min_i)
    unsolved <- unsolved[-which(unsolved == max_min_i)] # mark the max-min loc as solved
  }
  order
}

#' knn_indices
#'
#' Find the index of K nearest neighbors within a set of locations for a given query and distance function
#'
#' @param ordered_locs A matrix with one row per location, where locations are ordered using max_min ordering.
#' @param query Find the nearest neighbors to this location
#' @param n_neighbors The number of neighbors to find (K)
#' @param dist_func Any distance function with a signature of dist(query_location, locations_matrix)
#'
#' @return A vector containing the indices of the neighbors
#' @noRd
knn_indices <- function(ordered_locs, query, n_neighbors, dist_func, dist_func_code) {
  if (dist_func_code == "custom") {
    dists <- dist_func(query, ordered_locs)
    dists_order <- order(dists)
    nearest_neighbors <- dists_order[1:n_neighbors]
    junk <- list(
      "indices" = nearest_neighbors,
      "distances" = dists[nearest_neighbors]
    )
  } else {
    cur.nn <- nn2(ordered_locs, query, n_neighbors)
    junk <- list("indices" = cur.nn$nn.idx, "distances" = cur.nn$nn.dists)
  }
  junk
}

#' sparseNN
#'
#' Find the index of K nearest neighbors within a set of locations for a given query and distance function
#'
#' @param ordered_locs A matrix with one row per location, where locations are ordered using max_min ordering.
#' @param n_neighbors The number of neighbors to find (K) for each location
#' @param dist_func Any distance function with a signature of dist(query_location, locations_matrix)
#'
#' @return A list containing two matrices, each with one row per location:
#' an indices matrix with the indices of nearest neighbors for each location, and a distance matrix with the associated distances
#' @noRd
sparseNN <- function(ordered_locs, n_neighbors, dist_func, dist_func_code, ordered_locs_pred = NULL) {
  #  ee <- min(apply(ordered_locs, 2, stats::sd))
  #  n <- nrow(ordered_locs)
  #  ordered_locs <- ordered_locs + matrix(
  #    ee * 1e-04 *
  #      stats::rnorm(n * ncol(ordered_locs)),
  #    n, ncol(ordered_locs)
  #  )
  ordered_locs <- noise_locs(ordered_locs)
  indices_matrix <- matrix(
    data = NA, nrow = nrow(ordered_locs),
    ncol = n_neighbors
  )
  distances_matrix <- matrix(
    data = NA, nrow = nrow(ordered_locs),
    ncol = n_neighbors
  )
  for (row in 1:n_neighbors) {
    # for the locations from 1 to n_neighbors, use the entire locs list to find the neighbors
    nn <- knn_indices(
      ordered_locs[1:(n_neighbors + 1), , drop = FALSE][-row, ,
        drop = FALSE
      ],
      ordered_locs[row, , drop = FALSE], n_neighbors,
      dist_func, dist_func_code
    )
    indices_matrix[row, 1:n_neighbors] <- nn$indices[1:n_neighbors]
    distances_matrix[row, 1:n_neighbors] <- nn$distances[1:n_neighbors]
  }
  for (row in (n_neighbors + 1):nrow(ordered_locs)) {
    # get the m nearest neighbors from the locs before this one in the max-min order
    nn <- knn_indices(
      ordered_locs[1:(row - 1), , drop = FALSE],
      ordered_locs[row, , drop = FALSE], n_neighbors,
      dist_func, dist_func_code
    )
    indices_matrix[row, 1:n_neighbors] <- nn$indices[1:n_neighbors]
    distances_matrix[row, 1:n_neighbors] <- nn$distances[1:n_neighbors]
  }
  if (!is.null(ordered_locs_pred)) {
    indices_matrix_pred <- matrix(
      data = NA, nrow = nrow(ordered_locs_pred),
      ncol = n_neighbors
    )
    distances_matrix_pred <- matrix(
      data = NA, nrow = nrow(ordered_locs_pred),
      ncol = n_neighbors
    )
    for (row in seq_len(nrow(ordered_locs_pred))) {
      nn <- knn_indices(
        ordered_locs,
        ordered_locs_pred[row, , drop = FALSE], n_neighbors,
        dist_func, dist_func_code
      )
      indices_matrix_pred[row, 1:n_neighbors] <- nn$indices[1:n_neighbors]
      distances_matrix_pred[row, 1:n_neighbors] <- nn$distances[1:n_neighbors]
    }
    indices_matrix <- rbind(indices_matrix, indices_matrix_pred)
    distances_matrix <- rbind(distances_matrix, distances_matrix_pred)
  }
  list("indices" = indices_matrix, "distances" = distances_matrix)
}

# This partitions the vector q(i) of nearest neighbors into the subsets
# q_y(i) and q_z(i) as described in the second full paragraph on page
# 17 of Kyle's technical report.
calc.q <- function(nn.obj, firstind.pred) {
  m <- ncol(nn.obj)
  n <- nrow(nn.obj)

  q.y <- list(length = n)
  q.z <- list(length = n)
  q.y[[1]] <- NULL
  q.z[[1]] <- NULL

  for (i in 2:(m + 1)) {
    q.y[[i]] <- 1:(i - 1)
    q.z[[i]] <- NULL
  }

  for (i in (m + 1):n) {
    cur.q <- nn.obj[i, ]
    best.k <- cur.q[1]
    best.qy <- intersect(q.y[[best.k]], cur.q)
    for (j in 2:m) {
      cur.k <- cur.q[j]
      cur.qy <- intersect(q.y[[cur.k]], cur.q)
      if (length(cur.qy) > length(best.qy) && cur.k < firstind.pred) {
        best.k <- cur.k
        best.qy <- cur.qy
      }
    }
    q.y[[i]] <- union(best.k, best.qy)
    q.z[[i]] <- setdiff(cur.q, q.y[[i]])
    latent.pred <- q.z[[i]][q.z[[i]] >= firstind.pred]
    if (length(latent.pred) > 0) {
      q.y[[i]] <- union(q.y[[i]], latent.pred)
      q.z[[i]] <- setdiff(q.z[[i]], latent.pred)
    }
  }
  list(q.y = q.y, q.z = q.z)
}

#' Specify a multivariate Vecchia approximation
#'
#' Specifies a multivariate Vecchia approximation for later use in likelihood
#' evaluation or prediction. This function does not depend on parameter values,
#' and only has to be run once before repeated likelihood evaluations. This
#' function is a multivariate version of
#' \code{\link[GPvecchia]{vecchia_specify}}.
#'
#' @param locs.list List of observed locations. Each each element should be a
#' matrix containing the locs for the corresponding outcome variable.
#' @param m Number of nearby points to condition on.
#' @param locs.list.pred List of locations at which to make predictions. Each
#' element should be a matrix containing the locs for the corresponding outcome
#' variable.
#' @param dist.func Any distance function with a signature of
#' dist(query_location, locations_matrix). Defaults to Euclidean distance.
#' @param ordering.pred Should "obspred" or "general" ordering be used for
#' prediction? See \code{\link[GPvecchia]{vecchia_specify}}. Defaults to
#' "obspred".
#' @param pred.cond Should prediction conditioning be "general" or
#' "independent"? See \code{\link[GPvecchia]{vecchia_specify}}. Defaults to
#' "independent".
#'
#' @details This function should produce identical results to
#' \code{\link[GPvecchia]{vecchia_specify}} for univariate problems, although
#' it has fewer options. We recommend that
#' \code{\link[GPvecchia]{vecchia_specify}} be used in the univariate case.
#'
#' @return An object that specifies the multivariate Vecchia approximation for
#' later use in likelihood evaluation or prediction.
#'
#' @seealso \code{\link[GPvecchia]{vecchia_specify}}
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
vecchia_Mspecify <- function(locs.list, m, locs.list.pred = NULL,
  dist.func = NULL,
  ordering.pred = c("obspred", "general"),
  pred.cond = c("independent", "general")) {
  ordering.pred <- match.arg(ordering.pred)
  pred.cond <- match.arg(pred.cond)

  dist.func.code <- "custom"
  if (is.null(dist.func)) {
    dist.func <- fields::rdist
    dist.func.code <- "rdist"
  }

  P <- length(locs.list)

  locs <- NULL
  ndx <- NULL
  for (i in 1:P) {
    if (!is.matrix(locs.list[[i]])) {
      stop("Each element of locs.list must be a matrix")
    }
    locs <- rbind(locs, locs.list[[i]])
    ndx <- c(ndx, rep(i, nrow(locs.list[[i]])))
  }
  n <- nrow(locs)

  if (!is.null(locs.list.pred)) {
    if (length(locs.list.pred) != P) {
      stop("locs.list and locs.list.pred must have the same length")
    }
    locs.pred <- NULL
    ndx.pred <- NULL
    for (i in 1:P) {
      if (!is.null(locs.list.pred[[i]])) {
        if (!is.matrix(locs.list.pred[[i]])) {
          stop("Each element of locs.list.pred must be a matrix or NULL")
        }
        locs.pred <- rbind(locs.pred, locs.list.pred[[i]])
        ndx.pred <- c(ndx.pred, rep(i, nrow(locs.list.pred[[i]])))
      }
    }
    locs.all <- rbind(locs, locs.pred)
    ndx.all <- c(ndx, ndx.pred)
    observed.obspred <- c(rep(TRUE, n), rep(FALSE, nrow(locs.pred)))
  } else {
    locs.all <- locs
    ndx.all <- ndx
    observed.obspred <- rep(TRUE, n)
    ordering.pred <- "general"
  }

  if (dist.func.code == "custom") {
    if (!is.null(locs.list.pred)) {
      stop("Only Euclidean distance currently supported for prediction")
    }
    loc.order <- max_min_ordering(locs.all, dist.func)
    loc.order <- c(unique(loc.order), setdiff(1:n, loc.order))
  } else {
    if (is.null(locs.list.pred) || ordering.pred == "general") {
      loc.order <- GPvecchia::order_maxmin_exact(locs.all)
      # I am not sure why the next two lines are here. I added them because
      # similar code exists in the GPvecchia package. But I don't know why
      # they did this.
      cutoff <- min(n, 9)
      loc.order <- c(
        loc.order[1], loc.order[-seq(1, cutoff)],
        loc.order[2:cutoff]
      )
      ord <- loc.order
      ord.z <- loc.order[loc.order <= n]
    } else {
      loc.order <- GPvecchia::order_maxmin_exact_obs_pred(
        locs,
        locs.pred
      )
      ord.z <- loc.order$ord
      ord <- c(ord.z, loc.order$ord_pred + nrow(locs))
    }
    olocs <- locs.all[ord, , drop = FALSE]
    ondx <- ndx.all[ord]
    obs <- observed.obspred[ord]
  }

  # Note that the corresponding function in the GPvecchia package
  # is non-deterministic, so there may be some slight differences
  # between the output of this function and the output of createU
  # in the GPvecchia package.
  if (is.null(locs.list.pred) || pred.cond == "general") {
    nn.mat <- sparseNN(olocs, m, dist.func, dist.func.code)
  } else {
    nn.mat <- sparseNN(
      olocs[1:n, , drop = FALSE], m, dist.func, dist.func.code,
      olocs[-(1:n), , drop = FALSE]
    )
  }
  last.obs <- max(which(obs))
  q.list <- calc.q(nn.mat$indices, last.obs + 1)

  list(
    locsord = olocs, obs = obs, ord = ord, ord.z = ord.z,
    ord.pred = ordering.pred, cond.yz = "SGV", conditioning = "NN",
    P = P, ondx = ondx, dist.func = dist.func,
    dist.func.code = dist.func.code, q.list = q.list,
    n.neighbors = m
  )
}

#' Create the sparse triangular matrix U for multivariate Vecchia models
#'
#' This creates the sparse triangular matrix U for multivariate Vecchia
#' models. This matrix can be used to estimate the likelihood or transform
#' the data to be iid. This function is a multivariate version of
#' \code{\link[GPvecchia]{createU}}.
#'
#' @param vec.approx Object returned by \code{\link{vecchia_Mspecify}}.
#' @param params Vector of covariance parameters. See
#' \code{\link{create_param_sequence}} or the examples below for details
#' about the format of this vector.
#' @param cov_func The function used to compute the covariance between two
#' observations. Defaults to a Matern model.
#'
#' @details This function will be much slower if a non-default cov_func is
#' specified. More importantly, there is no guarantee that the resulting
#' covariance matrices will be positive definite. We recommend using the
#' default (Matern) covariance function unless you know exactly what you are
#' doing. See Apanasovich et al. (2012) for a description of how the
#' cross-covariances are computed.
#'
#' @return A list containing the sparse upper trianguler U, plus additional
#' objects required for other functions.
#'
#' @seealso \code{\link[GPvecchia]{createU}}, \code{\link{vecchia_Mspecify}},
#' \code{\link{create_param_sequence}}
#'
#' @references
#' \itemize{
#' \item Apanasovich, T.V., Genton, M.G. and Sun, Y. "A valid Matérn class of
#' cross-covariance functions for multivariate random fields with any number
#' of components", Journal of the American Statistical Association (2012)
#' 107(497):180-193.
#' \item Katzfuss, M., and Guinness, J. "A general framework for Vecchia
#' approximations of Gaussian processes", Statistical Science (2021)
#' 36(1):124-141.
#' }
#'
#' @useDynLib PrestoGP
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
#' soil.u <- createUMultivariate(soil.va, params)
createUMultivariate <- function(vec.approx, params, cov_func = NULL) {
  if (is.null(cov_func)) {
    cov_func <- fields::Matern
  }

  dist_func <- vec.approx$dist.func
  P <- vec.approx$P
  q.list <- vec.approx$q.list
  olocs <- vec.approx$locsord
  n <- nrow(olocs)
  ondx <- vec.approx$ondx

  param.seq <- create_param_sequence(P)
  param.sequence.begin <- param.seq[, 1]
  param.sequence.end <- param.seq[, 2]

  sig2 <- params[param.sequence.begin[1]:param.sequence.end[1]]
  rangep <- params[param.sequence.begin[2]:param.sequence.end[2]]
  smoothness <- params[param.sequence.begin[3]:param.sequence.end[3]]
  nugget <- params[param.sequence.begin[4]:param.sequence.end[4]]

  rho <- params[param.sequence.begin[5]:param.sequence.end[5]]
  rho.mat <- matrix(0, nrow = P, ncol = P)
  rho.mat[upper.tri(rho.mat, diag = FALSE)] <- rho
  rho.mat <- rho.mat + t(rho.mat)
  diag(rho.mat) <- 1

  if (vec.approx$dist.func.code == "rdist") {
    uvec <- rep(NA, 7)

    uvec[1] <- sig2[ondx[1]]^(-1 / 2)
    uvec[3] <- nugget[ondx[1]]^(-1 / 2)
    uvec[2] <- -1 * uvec[3]
    uvec[7] <- nugget[ondx[2]]^(-1 / 2)
    uvec[6] <- -1 * uvec[7]

    vii <- smoothness[ondx[1]]
    vjj <- smoothness[ondx[2]]
    vij <- (vii + vjj) / 2
    aii <- 1 / rangep[ondx[1]]
    ajj <- 1 / rangep[ondx[2]]
    aij <- sqrt((aii^2 + ajj^2) / 2)
    K1 <- rho.mat[ondx[1], ondx[2]] * sqrt(sig2[ondx[1]]) * sqrt(sig2[ondx[2]]) *
      aii^vii * ajj^vjj * gamma(vij) / (aij^(2 * vij) * sqrt(gamma(vii) * gamma(vjj))) *
      cov_func(dist_func(olocs[1, , drop = FALSE], olocs[2, , drop = FALSE], ),
        smoothness = vij, alpha = aij
      )
    K2 <- sig2[ondx[1]]
    bi <- K1 / K2
    ri <- sig2[ondx[2]] - bi * K1

    uvec[5] <- ri^(-1 / 2)
    uvec[4] <- -1 * bi * ri^(-1 / 2)

    vijs <- outer(smoothness, smoothness, "+") / 2
    aijs <- outer(rangep, rangep, function(x, y) sqrt((1 / x^2 + 1 / y^2) / 2))
    gammas <- outer(smoothness, smoothness, function(x, y) gamma((x + y) / 2) / sqrt(gamma(x) * gamma(y)))
    expprod <- outer(rangep^(-smoothness), rangep^(-smoothness), "*")
    sigs <- outer(sig2, sig2, function(x, y) sqrt(x) * sqrt(y))

    full_const <- sigs * gammas * expprod * rho.mat / (aijs^(2 * vijs))

    if (sum(is.na(full_const)) > 0) {
      browser()
    }

    m <- length(q.list$q.y[[n]]) + length(q.list$q.z[[n]])
    q.list$q.y <- lapply(q.list$q.y, function(x) c(x, rep(NA, m - length(x))))
    q.list$q.z <- lapply(q.list$q.z, function(x) c(x, rep(NA, m - length(x))))
    cur.qys <- do.call(cbind, q.list$q.y)
    cur.qzs <- do.call(cbind, q.list$q.z)
    # browser()
    U <- createU_helper_mat(olocs, ondx, cur.qys, cur.qzs, vijs, aijs, full_const, nugget, sig2, uvec)
    # U <- sparseMatrix(i=U1[1,], j=U1[2,], x = U1[3,], triangular = TRUE)
  }

  if (vec.approx$dist.func.code == "custom") {
    # U <- matrix(0, nrow=2*n, ncol=2*n)
    # U <- sparseMatrix(i = 1, j = 1, x = sig2[ondx[1]]^(-1/2), dims = c(2*n, 2*n), triangular = T)
    U1 <- matrix(ncol = 3, nrow = 7)
    U1[1, ] <- c(1, 1, sig2[ondx[1]]^(-1 / 2))
    U1[2, ] <- c(2, 2, nugget[ondx[1]]^(-1 / 2))
    U1[3, ] <- c(1, 2, -U1[2, 3])
    U1[4, ] <- c(4, 4, nugget[ondx[2]]^(-1 / 2))
    U1[5, ] <- c(3, 4, -U1[4, 3])

    # U[1,1] <- sig2[ondx[1]]^(-1/2)
    # U[2,2] <- nugget[ondx[1]]^(-1/2)
    # U[1,2] <- -1*U[2,2]
    # U[4,4] <- nugget[ondx[2]]^(-1/2)
    # U[3,4] <- -1*U[4,4]
    vii <- smoothness[ondx[1]]
    vjj <- smoothness[ondx[2]]
    vij <- (vii + vjj) / 2
    aii <- 1 / rangep[ondx[1]]
    ajj <- 1 / rangep[ondx[2]]
    aij <- sqrt((aii^2 + ajj^2) / 2)
    K1 <- rho.mat[ondx[1], ondx[2]] * sqrt(sig2[ondx[1]]) * sqrt(sig2[ondx[2]]) *
      aii^vii * ajj^vjj * gamma(vij) / (aij^(2 * vij) * sqrt(gamma(vii) * gamma(vjj))) *
      cov_func(dist_func(olocs[1, , drop = FALSE], olocs[2, , drop = FALSE], ),
        smoothness = vij, alpha = aij
      )
    K2 <- sig2[ondx[1]]
    bi <- K1 / K2
    ri <- sig2[ondx[2]] - bi * K1
    U1[6, ] <- c(3, 3, ri^(-1 / 2))
    U1[7, ] <- c(1, 3, -1 * bi * ri^(-1 / 2))
    # U[3,3] <- ri^(-1/2)
    # U[1,3] <- -1*bi*ri^(-1/2)
    i <- NULL # lintr requirement
    U2 <- foreach(i = 3:n, .combine = rbind) %dopar% {
      # U[2*i,2*i] <- nugget[ondx[i]]^(-1/2)
      # U[2*i-1,2*i] <- -1*U[2*i,2*i]
      cur.qy <- q.list$q.y[[i]]
      cur.qz <- q.list$q.z[[i]]
      cur.q <- c(cur.qy, cur.qz)
      nq <- length(cur.q)
      cur.U <- matrix(nrow = 3 + nq, ncol = 3)
      # Computing K1 and K2 will be slow because I don't think I can
      # vectorize these calculations unless there is a common smoothness
      # parameter, which generally will not be true in the multivariate
      # case. Rewriting these calculations in C++ is probably going
      # to be the simplest way to speed this up. Currently it is going
      # to be too slow for large data sets.
      K1 <- rep(NA, nq)
      for (j in 1:nq) {
        vii <- smoothness[ondx[i]]
        vjj <- smoothness[ondx[cur.q[j]]]
        vij <- (vii + vjj) / 2
        aii <- 1 / rangep[ondx[i]]
        ajj <- 1 / rangep[ondx[cur.q[j]]]
        aij <- sqrt((aii^2 + ajj^2) / 2)
        # The funky multiplier before the covariance function is necessary
        # is necessary to ensure that the final covariance matrix is
        # positive definite. See equation (9) in Apanasovich (2011).
        K1[j] <- rho.mat[ondx[i], ondx[cur.q[j]]] *
          sqrt(sig2[ondx[i]]) * sqrt(sig2[ondx[cur.q[j]]]) *
          aii^vii * ajj^vjj * gamma(vij) / (aij^(2 * vij) * sqrt(gamma(vii) * gamma(vjj))) *
          cov_func(
            dist_func(
              olocs[i, , drop = FALSE],
              olocs[cur.q[j], , drop = FALSE]
            ),
            smoothness = vij, alpha = aij
          )
      }
      K2 <- matrix(nrow = nq, ncol = nq)
      for (j in 1:nq) {
        for (k in j:nq) {
          vii <- smoothness[ondx[cur.q[j]]]
          vjj <- smoothness[ondx[cur.q[k]]]
          vij <- (vii + vjj) / 2
          aii <- 1 / rangep[ondx[cur.q[j]]]
          ajj <- 1 / rangep[ondx[cur.q[k]]]
          aij <- sqrt((aii^2 + ajj^2) / 2)
          K2[j, k] <- rho.mat[ondx[cur.q[j]], ondx[cur.q[k]]] *
            sqrt(sig2[ondx[cur.q[j]]]) * sqrt(sig2[ondx[cur.q[k]]]) *
            aii^vii * ajj^vjj * gamma(vij) /
            (aij^(2 * vij) * sqrt(gamma(vii) * gamma(vjj))) *
            cov_func(
              dist_func(
                olocs[cur.q[j], , drop = FALSE],
                olocs[cur.q[k], , drop = FALSE]
              ),
              smoothness = vij, alpha = aij
            )
          if (j != k) {
            K2[k, j] <- K2[j, k]
          }
        }
      }
      K2 <- K2 + diag(c(rep(0, length(cur.qy)), nugget[ondx[cur.qz]]),
        nrow = nq, ncol = nq
      )
      bi <- t(solve(K2, K1))
      ri <- sig2[ondx[i]] - bi %*% K1
      cur.br <- -1 * as.vector(bi) * c(ri^(-1 / 2))
      cur.U[1, ] <- c(2 * i - 1, 2 * i - 1, ri^(-1 / 2))
      cur.U[2, ] <- c(2 * i, 2 * i, nugget[ondx[i]]^(-1 / 2))
      cur.U[3, ] <- c(2 * i - 1, 2 * i, -cur.U[2, 3])
      temp_len.y <- length(cur.qy)
      # U[(ind+3), ] <- c(2*cur.qy-1,2*i-1, cur.br[1:length(cur.qy)])
      # ind <- ind + 4
      cur.U[4:(4 + temp_len.y - 1), 1] <- 2 * cur.qy - 1
      # show("NO")
      cur.U[4:(4 + temp_len.y - 1), 2] <- 2 * i - 1
      cur.U[4:(4 + temp_len.y - 1), 3] <- cur.br[1:temp_len.y]
      # show("YES")
      # U[2*i-1,2*i-1] <- ri^(-1/2)
      # U[2*cur.qy-1,2*i-1] <- cur.br[1:length(cur.qy)]
      temp_len.z <- length(cur.qz)
      if (temp_len.z > 0) {
        # U[2*cur.qz,2*i-1] <- cur.br[-(1:length(cur.qy))]
        cur.U[(4 + temp_len.y):(3 + nq), 2] <- 2 * i - 1
        cur.U[(4 + temp_len.y):(3 + nq), 1] <- 2 * cur.qz
        cur.U[(4 + temp_len.y):(3 + nq), 3] <- cur.br[-(1:temp_len.y)]
      }
      cur.U
    }
    U <- rbind(U1, U2)
    # browser()
    U <- sparseMatrix(
      i = U[, 1], j = U[, 2], x = U[, 3], dims = c(2 * n, 2 * n),
      triangular = TRUE
    )
  }
  # I think the code below only works for obspred ordering. This probably
  # needs to be fixed when (if) we support other orderings.
  if (sum(!vec.approx$obs) > 0) {
    if (vec.approx$ord.pred != "obspred") {
      stop("Currently only obspred ordering is supported")
    } else {
      drop.seq <- seq(from = (2 * sum(vec.approx$obs) + 2), to = ncol(U), by = 2)
      U <- U[-drop.seq, -drop.seq]
    }
  }
  latent <- rep(TRUE, nrow(U))
  latent[seq(from = 2, to = (2 * sum(vec.approx$obs)), by = 2)] <- FALSE
  list(
    U = U, latent = latent, ord = vec.approx$ord, obs = vec.approx$obs,
    ord.pred = vec.approx$ord.pred, ord.z = vec.approx$ord.z,
    cond.yz = vec.approx$cond.yz, ic0 = FALSE
  )
}
