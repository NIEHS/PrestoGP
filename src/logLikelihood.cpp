#include "utilFuncs.h"

// Transform the log Matern parameters back to the original
// [[Rcpp::export]]
arma::vec unlog_params(const arma::vec& logparams,
                        const arma::mat& param_seq,
                        int P) {

  // --- Segment 1: exp(logparams[1:param.seq[2, 2]]) ---
  // R: logparams[1 : param.seq(2,2)]  -> C++: logparams[0 : param.seq(1,1)-1]
  arma::uword end1 = param_seq(1, 1) - 1;
  arma::vec part1 = arma::exp(logparams.subvec(0, end1));

  // --- Segment 2: inv.logit(logparams[param.seq[3,1]:param.seq[3,2]], 0, 2.5) ---
  arma::uword start2 = param_seq(2, 0) - 1;
  arma::uword end2   = param_seq(2, 1) - 1;
  arma::vec seg2  = logparams.subvec(start2, end2);
  arma::vec part2 = 2.5 / (1.0 + arma::exp(-seg2));   // min = 0, max = 2.5

  // --- Segment 3: exp(logparams[param.seq[4,1]:param.seq[4,2]]) ---
  arma::uword start3 = param_seq(3, 0) - 1;
  arma::uword end3   = param_seq(3, 1) - 1;
  arma::vec part3 = arma::exp(logparams.subvec(start3, end3));

  arma::vec params = arma::join_cols(arma::join_cols(part1, part2), part3);

  // --- Segment 4: conditional on P ---
  if (P > 1) {
    // tanh(logparams[param.seq[5,1]:param.seq[5,2]])
    arma::uword start4 = param_seq(4, 0) - 1;
    arma::uword end4   = param_seq(4, 1) - 1;
    arma::vec part4 = arma::tanh(logparams.subvec(start4, end4));
    params = arma::join_cols(params, part4);
  } else {
    arma::vec one(1);
    one(0) = 1.0;
    params = arma::join_cols(params, one);
  }

  return params;
}

// ---------------------------------------------------------------------------
// Helper: solve L * x = b where L is sparse LOWER triangular (Matrix::solve
// with system = "L"). Implemented as column-oriented forward substitution
// directly on Armadillo's CSC storage. This avoids depending on SuperLU
// (arma::spsolve requires ARMA_USE_SUPERLU + linking the actual SuperLU
// library, which causes "undefined symbol" load errors if not linked
// correctly) and is also faster than a general sparse LU factorization
// since it exploits the triangular structure directly.
//
// Assumes L is square and has no nonzero entries strictly above the
// diagonal (true for a Cholesky-type factor).
// ---------------------------------------------------------------------------
static arma::vec sp_lower_tri_solve(const arma::sp_mat& L, const arma::vec& b) {
  const arma::uword n = L.n_rows;
  arma::vec x = b;                       // working vector, updated in place
  arma::vec d = arma::vec(L.diag());     // cache diagonal for O(1) lookups

  for (arma::uword j = 0; j < n; ++j) {
    x(j) /= d(j);
    for (arma::sp_mat::const_col_iterator it = L.begin_col(j); it != L.end_col(j); ++it) {
      const arma::uword i = it.row();
      if (i > j) {
        x(i) -= (*it) * x(j);
      }
    }
  }

  return x;
}

// ---------------------------------------------------------------------------
// Helper: extract a subset of rows (0-based indices) from a sparse matrix.
//
// Armadillo's support for non-contiguous sparse row/col indexing (e.g.
// A.rows(uvec)) is inconsistent across versions, so we build the submatrix
// directly from the non-zero entries. This is version-safe and O(nnz).
// ---------------------------------------------------------------------------
static arma::sp_mat sp_select_rows(const arma::sp_mat& A, const arma::uvec& row_idx) {
  const arma::uword n_new_rows = row_idx.n_elem;
  const arma::uword n_cols     = A.n_cols;

  // old row index -> new row index (-1 if not selected)
  arma::ivec row_map(A.n_rows);
  row_map.fill(-1);
  for (arma::uword i = 0; i < n_new_rows; ++i) {
    row_map(row_idx(i)) = static_cast<int>(i);
  }

  std::vector<arma::uword> new_rows;
  std::vector<arma::uword> new_cols;
  std::vector<double>      new_vals;
  new_rows.reserve(A.n_nonzero);
  new_cols.reserve(A.n_nonzero);
  new_vals.reserve(A.n_nonzero);

  for (arma::sp_mat::const_iterator it = A.begin(); it != A.end(); ++it) {
    const int mapped = row_map(it.row());
    if (mapped >= 0) {
      new_rows.push_back(static_cast<arma::uword>(mapped));
      new_cols.push_back(it.col());
      new_vals.push_back(*it);
    }
  }

  arma::umat locations(2, new_vals.size());
  for (size_t k = 0; k < new_vals.size(); ++k) {
    locations(0, k) = new_rows[k];
    locations(1, k) = new_cols[k];
  }
  arma::vec values(new_vals.data(), new_vals.size());

  return arma::sp_mat(locations, values, n_new_rows, n_cols);
}

// ---------------------------------------------------------------------------
// vecchia_likelihood_U (C++ version)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
double vecchia_likelihood_U_cpp(const arma::vec& z, List U_obj) {

  // --- unpack U.obj ----------------------------------------------------
  arma::sp_mat  U       = as<arma::sp_mat>(U_obj["U"]);   // requires dgCMatrix -> sp_mat converter
  LogicalVector latentR = U_obj["latent"];
  IntegerVector ordzR   = U_obj["ord.z"];

  const arma::uword n_total = latentR.size();

  arma::uvec latent(n_total);
  for (arma::uword i = 0; i < n_total; ++i) {
    latent(i) = latentR[i] ? 1u : 0u;
  }

  const arma::uword n_obs = ordzR.size();
  arma::uvec ord_z(n_obs);
  for (arma::uword i = 0; i < n_obs; ++i) {
    ord_z(i) = static_cast<arma::uword>(ordzR[i] - 1); // R is 1-based -> 0-based
  }

  // --- zord <- z[U.obj$ord.z] -------------------------------------------
  arma::vec zord = z.elem(ord_z);

  // --- row indices for latent / observed (!latent) ------------------------
  arma::uvec obs_idx    = arma::find(latent == 0);
  arma::uvec latent_idx = arma::find(latent == 1);

  // --- const <- sum(!latent) * log(2 * pi) --------------------------------
  const double const_term = static_cast<double>(obs_idx.n_elem) * std::log(2.0 * M_PI);

  // --- z1 <- crossprod(U[!latent, ], zord) = t(U_obs) %*% zord ------------
  arma::sp_mat U_obs = sp_select_rows(U, obs_idx);
  arma::vec z1 = arma::vec(U_obs.t() * zord);

  // --- quadform.num <- sum(z1^2) ------------------------------------------
  const double quadform_num = arma::accu(arma::square(z1));

  // --- logdet.num <- -2 * sum(log(diag(U))) -------------------------------
  arma::vec diagU = arma::vec(U.diag());
  const double logdet_num = -2.0 * arma::accu(arma::log(diagU));

  double quadform_denom = 0.0;
  double logdet_denom   = 0.0;

  if (latent_idx.n_elem > 0) {
    // U.y <- U[latent, ]
    arma::sp_mat U_y = sp_select_rows(U, latent_idx);

    // z2 <- as.numeric(U.y %*% z1)
    arma::vec z2 = arma::vec(U_y * z1);

    // V.ord <- U2V_cpp(U.obj)
    arma::sp_mat V_ord = U2V_cpp(U_obj);

    // z3 <- solve(V.ord, rev(z2), system = "L")
    arma::vec z2_rev = arma::flipud(z2);
    arma::vec z3 = sp_lower_tri_solve(V_ord, z2_rev); // forward substitution, no SuperLU needed

    // quadform.denom <- sum(z3^2)
    quadform_denom = arma::accu(arma::square(z3));

    // logdet.denom <- -2 * sum(log(diag(V.ord)))
    arma::vec diagV = arma::vec(V_ord.diag());
    logdet_denom = -2.0 * arma::accu(arma::log(diagV));
  }

  const double neg2loglik = logdet_num - logdet_denom + quadform_num - quadform_denom + const_term;
  const double loglik = -neg2loglik / 2.0;

  return loglik;
}

//' Evaluation of the multivariate Vecchia likelihood
//'
//' This function is used to evaluate the multivariate Vecchia likelihood.
//'
//' @param z The observed data.
//' @param vecchia_approx A Vecchia object returned by
//' \code{\link{vecchia_Mspecify}}.
//' @param covparams Vector of covariance parameters. See
//' \code{\link{create_param_sequence}} or the examples below for details
//' about the format of this vector.
//'
//' @return The log likelihood implied by the multivariate Vecchia
//' approximation.
//'
//' @seealso \code{\link[GPvecchia]{vecchia_likelihood}},
//' \code{\link{vecchia_Mspecify}}, \code{\link{create_param_sequence}}
//'
//' @references
//' \itemize{
//' \item Katzfuss, M., and Guinness, J. "A general framework for Vecchia
//' approximations of Gaussian processes", Statistical Science (2021)
//' 36(1):124-141.
//' }
//'
//' @export
//' @examples
//' data(soil)
//' soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
//' locs <- as.matrix(soil[,1:2])
//' locsm <- list()
//' locsm[[1]] <- locsm[[2]] <- locs
//' soil.va <- vecchia_Mspecify(locsm, m=10)
//'
//' pseq <- create_param_sequence(2)
//' # Initialize the vector of covariance parameters
//' params <- rep(NA, pseq[5,2])
//' # Sigma parameters:
//' params[pseq[1,1]:pseq[1,2]] <- c(100, 80)
//' # Scale parameters:
//' params[pseq[2,1]:pseq[2,2]] <- c(60, 50)
//' # Smoothness parameters:
//' params[pseq[3,1]:pseq[3,2]] <- c(0.5, 0.5)
//' # Nuggets:
//' params[pseq[4,1]:pseq[4,2]] <- c(30, 30)
//' # Correlation:
//' params[pseq[5,1]:pseq[5,2]] <- -0.9
//'
//' vecchia_Mlikelihood(rnorm(nrow(locs)), soil.va, params)
// [[Rcpp::export]]
double vecchia_Mlikelihood(arma::vec z, Rcpp::List vecchia_approx,
			   arma::vec covparams) {

    // Build the U object (equivalent to U.obj <- createUMultivariate(...))
    List U_obj = createUMultivariate(vecchia_approx, covparams);

    double junk;

    // Equivalent to R's try(..., silent = TRUE): if
    // vecchia_likelihood_U_cpp throws (e.g. because of a numerically
    // singular matrix), catch it and fall back to -Inf instead of
    // letting the exception propagate.
    try {
        junk = vecchia_likelihood_U_cpp(z, U_obj);
    } catch (...) {
        junk = -std::numeric_limits<double>::infinity();
    }

    return junk;
}

//  Input-
//  logparams: A numeric vector of length (4*P)+(4*choose(P,2)).
//             To construct these parameters we unlist a list of the 7 covariance
//             categories- in order: (1) marginal variances, (2) Marginal ranges,
//             (3) Marginal smoothness, (4) Nuggets, and
//             (5) cross-covariance correlation. These seven parameters are to be
//             created in a list. The variance, range, smoothness, and nugget
//             have P terms, and the correlations have choose(P,2)
//             terms. Use unlist() to create the vector of parameters.
//
// locs: list of the location coordinates, each outcome in y is a separate cell
//       in the list
// y    :  multivariate outcome, each out outcome in a separate entry in a list
// param_seq: The vector of parameter index sequences created by the function
//            create_param_sequence - Used to identify the beginning and end
//            index locations of each parameter.

// P <- length(y)
// transform the postively constrained parameters from log-space to normal-space
// [[Rcpp::export]]
double mvnegloglik(arma::vec logparams, Rcpp::List vecchia_approx,
		   arma::vec y, arma::mat param_seq, int P) {
  arma::vec params = unlog_params(logparams, param_seq, P);
  return -1 * vecchia_Mlikelihood(y, vecchia_approx, params);
}

// NOTE on indexing: param_seq, scaling, and ondx all carry values that
// were originally meant as 1-based R indices/labels. R's matrix/vector
// *positions* (row i, column j) become 0-based in Armadillo (row i-1,
// column j-1), but the *values stored inside* param_seq that are used
// to slice `params` are still 1-based R indices, so 1 must be
// subtracted from them whenever they're used to index an arma::vec.

// [[Rcpp::export]]
double mvnegloglik_ST(arma::vec logparams, List vecchia_approx, arma::vec y,
		      arma::mat param_seq, int P, arma::ivec scaling,
		      int nscale) {

    // params <- unlog_params(logparams, param.seq, P)
    arma::vec params = unlog_params(logparams, param_seq, P);

    // locs.scaled <- vecchia.approx$locsord
    arma::mat locs_scaled = as<arma::mat>(vecchia_approx["locsord"]);
    arma::ivec ondx        = as<arma::ivec>(vecchia_approx["ondx"]);

    for (int i = 1; i <= P; i++) {
        arma::uvec rows_idx = arma::find(ondx == i);
        for (int j = 1; j <= nscale; j++) {
            arma::uvec cols_idx = arma::find(scaling == j);

            // R: params[param.seq[2, 1] + nscale * (i - 1) + j - 1]
            // param.seq[2, 1] (R, 1-based) -> param_seq(1, 0) (arma, 0-based)
            int base_1based   = param_seq(1, 0);
            int idx_1based    = base_1based + nscale * (i - 1) + j - 1;
            double divisor    = params(idx_1based - 1); // -1 to convert to 0-based

            if (rows_idx.n_elem > 0 && cols_idx.n_elem > 0) {
                locs_scaled.submat(rows_idx, cols_idx) /= divisor;
            }
        }
    }

    // vecchia.approx$locsord <- locs.scaled
    // (clone so the original R list/object passed in is left untouched)
    List vecchia_approx_scaled = clone(vecchia_approx);
    vecchia_approx_scaled["locsord"] = locs_scaled;

    // Build:
    //   c(params[1:param.seq[1, 2]],
    //     rep(1, param.seq[2, 2] - param.seq[2, 1] + 1),
    //     params[param.seq[3, 1]:param.seq[5, 2]])

    // params[1:param.seq[1, 2]]  ->  subvec(0, param_seq(0,1) - 1)
    arma::vec part1 = params.subvec(0, param_seq(0, 1) - 1);

    // rep(1, param.seq[2, 2] - param.seq[2, 1] + 1)
    int ones_len = param_seq(1, 1) - param_seq(1, 0) + 1;
    arma::vec part2(ones_len, arma::fill::ones);

    // params[param.seq[3, 1]:param.seq[5, 2]] -> subvec(param_seq(2,0)-1, param_seq(4,1)-1)
    arma::vec part3 = params.subvec(param_seq(2, 0) - 1, param_seq(4, 1) - 1);

    arma::vec covparams = arma::join_cols(arma::join_cols(part1, part2), part3);

    // -1 * vecchia_Mlikelihood(y, vecchia.approx, covparams)
    return -1.0 * vecchia_Mlikelihood(y, vecchia_approx_scaled, covparams);
}
