#include "utilFuncs.h"

double matern1(const double &x, const double &smooth, const double &alpha) {
  if (x == 0.0) {
    return 1.0;
  } else if (smooth == 0.5) {
    return exp(-x * alpha);
  } else if (smooth == 1.5) {
    double normcon = (pow(alpha, 3) / sqrt(M_PI)) *
      (boost::math::tgamma(2.0) / boost::math::tgamma(1.5));
    return normcon * 0.5 * M_PI * pow(alpha, -3) * exp(-x * alpha) *
           (1 + x * alpha);
  } else if (smooth == 2.5) {
    double normcon = (pow(alpha, 5) / sqrt(M_PI)) *
      (boost::math::tgamma(3.0) / boost::math::tgamma(2.5));
    double scaled_x = x * alpha;
    return normcon * 0.125 * M_PI * pow(alpha, -5) * exp(-scaled_x) *
      (3 + 3 * scaled_x + scaled_x * scaled_x); 
  } else {
    double normcon = pow(2.0, 1.0 - smooth) / boost::math::tgamma(smooth);
    double scaled_x = x * alpha;
    return normcon * pow(scaled_x, smooth) *
      boost::math::cyl_bessel_k(smooth, scaled_x);
  }
}

// computes Euclidean distance between 2 vectors ([slightly slower] helper
// function for rdist in C++)
double dist1(const arma::rowvec &a, const arma::rowvec &b) {
  return arma::norm(a - b, 2);
}

arma::vec na_omit_c(arma::vec x) { return x(arma::find_finite(x)); }

// ---------------------------------------------------------------------------
// Helper: infinity norm of a matrix  (max absolute row sum),
//         matching R's  norm(M, "I")
// ---------------------------------------------------------------------------
static double norm_inf(const mat& M) {
  if (M.n_elem == 0) return 0.0;
  double val = 0.0;
  for (uword i = 0; i < M.n_rows; ++i) {
    double row_sum = accu(abs(M.row(i)));
    if (row_sum > val) val = row_sum;
  }
  return val;
}

// ---------------------------------------------------------------------------
// Helper: symmetrize  (X + X') / 2
// ---------------------------------------------------------------------------
static inline mat sym(const mat& X) {
  return 0.5 * (X + X.t());
}

//' Extract specific Matern parameters from a parameter sequence
//'
//' This function is used to obtain specific Matern parameters (e.g.,
//' range or smoothness) from the covparams slot of a PrestoGPModel object.
//'
//' @param P Number of outcome variables
//' @param ns Number of scale parameters
//'
//' @details This function is intended for advanced users who want to specify
//' the input Matern parameters for functions such as
//' \code{\link{vecchia_Mlikelihood}} or \code{\link{createUMultivariate}}.
//' To extract the Matern parameters from a fitted PrestoGP model, it is
//' strongly recommended to use \code{link{get_theta}} instead.
//'
//' @return A matrix with five rows and two columns as described below:
//' \describe{
//' \item{Row 1:}{Starting and ending indices for the sigma parameter(s)}
//' \item{Row 2:}{Starting and ending indices for the scale parameter(s)}
//' \item{Row 3:}{Starting and ending indices for the smoothness parameter(s)}
//' \item{Row 4:}{Starting and ending indices for the nugget(s)}
//' \item{Row 5:}{Starting and ending indices for the correlation parameter(s)}
//' }
//'
//' @seealso \code{\link{PrestoGPModel-class}}
//'
//' @references
//' \itemize{
//' \item Apanasovich, T.V., Genton, M.G. and Sun, Y. "A valid Matérn class of
//' cross-covariance functions for multivariate random fields with any number
//' of components", Journal of the American Statistical Association (2012)
//' 107(497):180-193.
//' \item Genton, M.G. "Classes of kernels for machine learning: a statistics
//' perspective", The Journal of Machine Learning Research (2001) 2:299-312.
//' }
//'
//' @export
//' @examples
//' # Space/elevation model
//' data(soil250, package="geoR")
//' y2 <- soil250[,7]               # predict pH level
//' X2 <- as.matrix(soil250[,c(4:6,8:22)])
//' # columns 1+2 are location coordinates; column 3 is elevation
//' locs2 <- as.matrix(soil250[,1:3])
//'
//' soil.vm2 <- new("VecchiaModel", n_neighbors = 10)
//' # fit separate scale parameters for location and elevation
//' soil.vm2 <- prestogp_fit(soil.vm2, y2, X2, locs2, scaling = c(1, 1, 2))
//'
//' pseq <- create_param_sequence(1, 2)
//' soil2.params <- soil.vm2@covparams
//' # sigma
//' soil2.params[pseq[1,1]:pseq[1,2]]
//' # scale parameters
//' soil2.params[pseq[2,1]:pseq[2,2]]
//' # smoothness parameter
//' soil2.params[pseq[3,1]:pseq[3,2]]
//' # nugget
//' soil2.params[pseq[4,1]:pseq[4,2]]
//'
//' # Multivariate model
//' ym <- list()
//' ym[[1]] <- soil250[,4] # predict sand/silt portion of the sample
//' ym[[2]] <- soil250[,5]
//' ym[[3]] <- soil250[,6]
//' Xm <- list()
//' Xm[[1]] <- Xm[[2]] <- Xm[[3]] <- as.matrix(soil250[,7:22])
//' locsm <- list()
//' locsm[[1]] <- locsm[[2]] <- locsm[[3]] <- as.matrix(soil250[,1:3])
//'
//' soil.mvm <-  new("MultivariateVecchiaModel", n_neighbors = 10)
//' soil.mvm <- prestogp_fit(soil.mvm, ym, Xm, locsm)
//'
//' pseq <- create_param_sequence(3, 2)
//' soil.params <- soil.mvm@covparams
//' # sigmas
//' soil.params[pseq[1,1]:pseq[1,2]]
//' # scale parameters
//' scale.seq <- pseq[2,1]:pseq[2,2]
//' # scale parameter for location, outcome 1
//' soil.params[scale.seq[1]]
//' # scale parameter for elevation, outcome 1
//' soil.params[scale.seq[2]]
//' # scale parameter for location, outcome 2
//' soil.params[scale.seq[3]]
//' # scale parameter for elevation, outcome 2
//' soil.params[scale.seq[4]]
//' # scale parameter for location, outcome 3
//' soil.params[scale.seq[5]]
//' # scale parameter for elevation, outcome 3
//' soil.params[scale.seq[6]]
//' # smoothness parameters
//' soil.params[pseq[3,1]:pseq[3,2]]
//' # nuggets
//' soil.params[pseq[4,1]:pseq[4,2]]
//' # correlation
//' soil.corr <- diag(2) / 2
//' soil.corr[upper.tri(soil.corr)] <- soil.params[pseq[5,1]:pseq[5,2]]
//' soil.corr <- soil.corr + t(soil.corr)
// [[Rcpp::export]]
arma::mat create_param_sequence(const double P, const double ns = 1) {

  double nk = P * (P - 1) / 2;
  if (nk == 0) {
    nk = 1;
  }

  // param.sequence.begin <- c(1, P + 1, seq(P * (ns + 1) + 1, length = 3, by = P))
  // seq(..., length=3, by=P) produces 3 values starting at P*(ns+1)+1, stepping by P
  double seq_start = P * (ns + 1) + 1;
  arma::vec param_sequence_begin = {
    1,
    P + 1,
    seq_start,
    seq_start + P,
    seq_start + 2 * P
  };

  // param.sequence.end <- c(P, ns * P, P, P, nk) |> cumsum()
  arma::vec end_raw = {P, ns * P, P, P, nk};
  arma::vec param_sequence_end(5);
  double running = 0;
  for (int i = 0; i < 5; i++) {
    running += end_raw[i];
    param_sequence_end[i] = running;
  }

  // cbind into a 5x2 matrix
  arma::mat param_sequence(5, 2);
  param_sequence.col(0) = param_sequence_begin;
  param_sequence.col(1) = param_sequence_end;

  return param_sequence;
}

// ---------------------------------------------------------------------------
// nearPD_cpp
//
// Arguments mirror the R function; conv.norm.type is fixed to "I" (infinity).
// Returns a named List with the same fields as the R "nearPD" S3 object:
//   mat         — the nearest positive-definite matrix (NumericMatrix)
//   eigenvalues — eigenvalues of mat (NumericVector)
//   corr        — whether corr = TRUE was requested (bool)
//   normF       — Frobenius norm of (original X − result)
//   iterations  — number of iterations performed
//   rel.tol     — final relative tolerance achieved
//   converged   — whether convergence criterion was satisfied
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List nearPD_cpp(
    arma::mat      x,                    // n × n approximately PD matrix
    bool           corr       = false,   // TRUE → force diagonal 1 (correlation)
    bool           keepDiag   = false,   // TRUE → preserve original diagonal
    bool           do2eigen   = true,    // TRUE → posdefify eigen step after main loop
    bool           doSym      = false,   // TRUE → symmetrise inside loop
    bool           doDykstra  = true,    // TRUE → use Dykstra's correction
    bool           only_values = false,  // TRUE → return eigenvalues only
    double         eig_tol    = 1e-6,   // relative positivity threshold
    double         conv_tol   = 1e-7,   // convergence tolerance
    double         posd_tol   = 1e-8,   // tolerance for the posdefify step
    int            maxit      = 100,     // maximum iterations
    bool           trace      = false    // print iteration info to Rcout
) {
    // -----------------------------------------------------------------------
    // 0.  Copy input, optionally symmetrize (caller should have done this;
    //     we mirror R's ensureSymmetry logic by always symmetrizing lightly).
    // -----------------------------------------------------------------------
  mat X = x;
  const uword n = X.n_rows;

  if (X.n_cols != n)
    stop("Input matrix must be square.");

  // Store original diagonal when keepDiag = TRUE
  vec diagX0(n);
  if (keepDiag)
    diagX0 = X.diag();

  // -----------------------------------------------------------------------
  // 1.  Initialise Dykstra correction matrix D_S (n × n zeros)
  // -----------------------------------------------------------------------
  mat D_S(n, n, fill::zeros);

  // -----------------------------------------------------------------------
  // 2.  Main Higham (2002) alternating-projection loop
  // -----------------------------------------------------------------------
  int  iter      = 0;
  bool converged = false;
  double conv    = datum::inf;

  while (iter < maxit && !converged) {

    mat Y = X;

    // R ← Y − D_S  (Dykstra's running correction)
    mat R = doDykstra ? (Y - D_S) : Y;

    // -------------------------------------------------------------------
    // 2a. Project onto PSD: eigendecomposition, threshold small λ, rebuild
    // -------------------------------------------------------------------
    vec  d;
    mat  Q;
    // symmetric = true → faster, more numerically stable
    if (!eig_sym(d, Q, R)) {
      stop("Eigendecomposition failed.");
    }

    // arma::eig_sym returns eigenvalues in *ascending* order;
    // R's eigen() returns them in *descending* order.
    // We need the largest eigenvalue to set the relative threshold.
    double d_max = d(n - 1);   // largest eigenvalue (ascending order)

    // Mask: keep only eigenvalues > eig_tol * d_max
    // (mirrors R:  p <- d > eig.tol * d[1])
    uvec p = find(d > eig_tol * d_max);
    if (p.is_empty())
      stop("Matrix seems negative semi-definite.");

    // X = Q_p * diag(d_p) * Q_p'
    mat  Qp  = Q.cols(p);          // n × |p|
    vec  dp  = d(p);               // |p| positive eigenvalues
    // Efficient:  (Q_p .* d_p') * Q_p'
    X = (Qp.each_row() % dp.t()) * Qp.t();

    // -------------------------------------------------------------------
    // 2b. Update Dykstra correction:  D_S ← X − R
    // -------------------------------------------------------------------
    if (doDykstra)
      D_S = X - R;

    // -------------------------------------------------------------------
    // 2c. Optional symmetrisation inside loop
    // -------------------------------------------------------------------
    if (doSym)
      X = sym(X);

    // -------------------------------------------------------------------
    // 2d. Project onto correlation / fixed-diagonal constraint
    // -------------------------------------------------------------------
    if (corr)
      X.diag().fill(1.0);
    else if (keepDiag)
      X.diag() = diagX0;

    // -------------------------------------------------------------------
    // 2e. Convergence check (infinity norm ratio, matching R default)
    // -------------------------------------------------------------------
    double normY  = norm_inf(Y);
    double normDiff = norm_inf(Y - X);
    conv = (normY > 0.0) ? (normDiff / normY) : normDiff;

    ++iter;

    if (trace)
      Rcout << "iter " << iter
	    << " : #{p}=" << p.n_elem
	    << ", ||Y-X||/||Y|| = " << conv << "\n";

    converged = (conv <= conv_tol);
  }

  if (!converged)
    warning("'nearPD()' did not converge in %d iterations", iter);

  // -----------------------------------------------------------------------
  // 3.  Optional posdefify eigen step  (mirrors sfsmisc::posdefify in R)
  // -----------------------------------------------------------------------
  vec d_final;
  mat Q_final;

  if (do2eigen || only_values) {
    if (!eig_sym(d_final, Q_final, X))
      stop("Eigendecomposition failed in posdefify step.");

    double Eps = posd_tol * std::abs(d_final(n - 1));  // d_final ascending
    
    if (d_final(0) < Eps) {                            // smallest eigenvalue
      // Clamp all eigenvalues below Eps
      d_final.elem(find(d_final < Eps)).fill(Eps);

      if (!only_values) {
	vec o_diag = X.diag();
	// X ← Q * diag(d_final) * Q'
	X = (Q_final.each_row() % d_final.t()) * Q_final.t();
	// Rescale to restore diagonal: D = sqrt(o_diag / diag(X))
	vec D_scale = sqrt(clamp(o_diag, Eps, datum::inf) / X.diag());
	// X ← D * X * D  (element-wise outer product of D_scale)
	X = X % (D_scale * D_scale.t());
      }
    }

    if (only_values) {
      // Return eigenvalues in descending order (as R does)
      return List::create(
			  Named("eigenvalues") = wrap(reverse(d_final))
			  );
    }

    // Enforce diagonal constraints again after posdefify
    if (corr)
      X.diag().fill(1.0);
    else if (keepDiag)
      X.diag() = diagX0;

  } else {
    // We still need eigenvalues for the return value
    if (!eig_sym(d_final, Q_final, X))
      stop("Final eigendecomposition failed.");
  }

  // -----------------------------------------------------------------------
  // 4.  Compute Frobenius norm of (original − result)
  // -----------------------------------------------------------------------
  mat X0 = x;
  double normF = norm(X0 - X, "fro");

  // -----------------------------------------------------------------------
  // 5.  Build return list  (eigenvalues in descending order, as R does)
  // -----------------------------------------------------------------------
  return List::create(
        Named("mat")         = wrap(X),
        Named("eigenvalues") = wrap(reverse(d_final)),
        Named("corr")        = corr,
        Named("normF")       = normF,
        Named("iterations")  = iter,
        Named("rel.tol")     = conv,
        Named("converged")   = converged
    );
}

arma::mat nearPD_wrap(arma::mat x) {
  return nearPD_cpp(x)["mat"];
}

// =====================================================================
// Sparse-matrix helpers
// =====================================================================

// Reverse row order AND column order of a sparse matrix in place of
// R's revMat(): mat[rev(seq_len(nrow)), rev(seq_len(ncol))].
// Implemented by remapping the (row, col) index of every stored
// nonzero -- O(nnz), no dense intermediate.
static arma::sp_mat revMat_sp(const arma::sp_mat& X) {
  if (X.n_rows == 0 || X.n_cols == 0) return X;

  const arma::uword nr = X.n_rows, nc = X.n_cols;
  std::vector<arma::uword> rows_out, cols_out;
  std::vector<double> vals_out;
  rows_out.reserve(X.n_nonzero);
  cols_out.reserve(X.n_nonzero);
  vals_out.reserve(X.n_nonzero);

  for (arma::sp_mat::const_iterator it = X.begin(); it != X.end(); ++it) {
    rows_out.push_back(nr - 1 - it.row());
    cols_out.push_back(nc - 1 - it.col());
    vals_out.push_back(*it);
  }

  arma::umat locations(2, rows_out.size());
  for (size_t i = 0; i < rows_out.size(); ++i) {
    locations(0, i) = rows_out[i];
    locations(1, i) = cols_out[i];
  }
  arma::vec values(vals_out.data(), vals_out.size());
  return arma::sp_mat(locations, values, nr, nc);
}

// Select an arbitrary (possibly non-contiguous) set of ROWS from a
// sparse matrix, keeping all columns and preserving the order given
// by row_idx (0-based). O(n_rows) for the index map + O(nnz) scan.
// (For CONTIGUOUS ranges, prefer X.rows(a,b) / X.submat(...), which
// Armadillo supports natively and more efficiently for sp_mat.)
static arma::sp_mat sp_select_rows(const arma::sp_mat& X, const arma::uvec& row_idx) {
  const arma::uword n_new = row_idx.n_elem;
  arma::ivec row_map(X.n_rows);
  row_map.fill(-1);
  for (arma::uword k = 0; k < n_new; ++k) row_map(row_idx(k)) = static_cast<int>(k);

  std::vector<arma::uword> rows_out, cols_out;
  std::vector<double> vals_out;
  rows_out.reserve(X.n_nonzero);
  cols_out.reserve(X.n_nonzero);
  vals_out.reserve(X.n_nonzero);

  for (arma::sp_mat::const_iterator it = X.begin(); it != X.end(); ++it) {
    int nr = row_map(it.row());
    if (nr >= 0) {
      rows_out.push_back(static_cast<arma::uword>(nr));
      cols_out.push_back(it.col());
      vals_out.push_back(*it);
    }
  }

  arma::umat locations(2, rows_out.size());
  for (size_t i = 0; i < rows_out.size(); ++i) {
    locations(0, i) = rows_out[i];
    locations(1, i) = cols_out[i];
  }
  arma::vec values(vals_out.data(), vals_out.size());
  return arma::sp_mat(locations, values, n_new, X.n_cols);
}

// Same idea, but selecting COLUMNS (keeping all rows).
static arma::sp_mat sp_select_cols(const arma::sp_mat& X, const arma::uvec& col_idx) {
  const arma::uword n_new = col_idx.n_elem;
  arma::ivec col_map(X.n_cols);
  col_map.fill(-1);
  for (arma::uword k = 0; k < n_new; ++k) col_map(col_idx(k)) = static_cast<int>(k);

  std::vector<arma::uword> rows_out, cols_out;
  std::vector<double> vals_out;
  rows_out.reserve(X.n_nonzero);
  cols_out.reserve(X.n_nonzero);
  vals_out.reserve(X.n_nonzero);

  for (arma::sp_mat::const_iterator it = X.begin(); it != X.end(); ++it) {
    int nc = col_map(it.col());
    if (nc >= 0) {
      rows_out.push_back(it.row());
      cols_out.push_back(static_cast<arma::uword>(nc));
      vals_out.push_back(*it);
    }
  }

  arma::umat locations(2, rows_out.size());
  for (size_t i = 0; i < rows_out.size(); ++i) {
    locations(0, i) = rows_out[i];
    locations(1, i) = cols_out[i];
  }
  arma::vec values(vals_out.data(), vals_out.size());
  return arma::sp_mat(locations, values, X.n_rows, n_new);
}

// Read U.obj$U as a sparse matrix, whether it arrives as an S4
// (dgCMatrix etc, handled directly by RcppArmadillo's as<arma::sp_mat>)
// or, defensively, as an ordinary dense R matrix.
static arma::sp_mat as_sparse_mat(SEXP x) {
  if (Rf_isS4(x)) {
    return as<arma::sp_mat>(x);
  }
  arma::mat Xd = as<arma::mat>(x);
  return arma::sp_mat(Xd);
}

// ---------------------------------------------------------------------
// Conversions between arma::sp_mat and Eigen::SparseMatrix<double>,
// needed because Armadillo has no sparse Cholesky of its own.
// ---------------------------------------------------------------------
static Eigen::SparseMatrix<double> arma_sp_to_eigen(const arma::sp_mat& X) {
  std::vector<Eigen::Triplet<double>> trip;
  trip.reserve(X.n_nonzero);
  for (arma::sp_mat::const_iterator it = X.begin(); it != X.end(); ++it) {
    trip.emplace_back(static_cast<int>(it.row()), static_cast<int>(it.col()), *it);
  }
  Eigen::SparseMatrix<double> Y(X.n_rows, X.n_cols);
  Y.setFromTriplets(trip.begin(), trip.end());
  return Y;
}

static arma::sp_mat eigen_to_arma_sp(const Eigen::SparseMatrix<double>& X) {
  const int nnz = X.nonZeros();
  arma::umat locations(2, nnz);
  arma::vec values(nnz);
  int idx = 0;
  for (int k = 0; k < X.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(X, k); it; ++it) {
      locations(0, idx) = it.row();
      locations(1, idx) = it.col();
      values(idx) = it.value();
      ++idx;
    }
  }
  return arma::sp_mat(locations, values, X.rows(), X.cols());
}

// ---------------------------------------------------------------------
// Sparse Cholesky with fallback to nearPD_cpp, replicating:
//   res <- try(V <- t(chol(X)), silent = TRUE)
//   if (inherits(res, "try-error")) {
//     X <- nearPD_cpp(X)$mat
//     V <- t(chol(X))
//   }
// but performed with Eigen::SimplicialLLT (sparse) as the primary
// attempt, since arma::chol only accepts dense matrices.
//
// NaturalOrdering<int> disables Eigen's fill-reducing permutation, so
// the returned factor's rows/cols line up 1-to-1 with X's -- this is
// required for correctness elsewhere in U2V (e.g. the zero-block
// stitching in the obspred branch).
//
// NOTE: the nearPD_cpp fallback densifies X. This assumes X (a W or A
// block, not the full U matrix) is small enough to do so -- true in
// the typical Vecchia setting where these blocks correspond to
// conditioning/observed sets, but worth keeping in mind if your
// blocks can also become enormous.
// ---------------------------------------------------------------------
static arma::sp_mat chol_with_fallback_sparse(const arma::sp_mat& X) {
  Eigen::SparseMatrix<double> X_eig = arma_sp_to_eigen(X);

  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower,
                        Eigen::NaturalOrdering<int>> llt;
  llt.compute(X_eig);

  if (llt.info() == Eigen::Success) {
    Eigen::SparseMatrix<double> L = llt.matrixL();
    return eigen_to_arma_sp(L);
  }

  // --- fallback: X is (numerically) singular / not PD ---------------
  arma::mat X_dense(X);                       // densify just this block
  List pd = nearPD_cpp(X_dense);
  arma::mat X_fixed = as<arma::mat>(pd["mat"]);

  arma::mat R;
  bool ok = false;
  try { ok = arma::chol(R, X_fixed); } catch (...) { ok = false; }
  if (!ok) {
    Rcpp::stop("U2V_cpp: Cholesky decomposition failed even after nearPD_cpp correction.");
  }
  arma::mat L_dense = R.t();
  return arma::sp_mat(L_dense);               // re-sparsify before returning
}

// =====================================================================
// RcppArmadillo/RcppEigen translation of GPvecchia's U2V(), sparse
// throughout: input U, all intermediates, and the returned V.ord are
// arma::sp_mat (wrap()'d back to R as a dgCMatrix -- there is no
// Armadillo "triangularMatrix" class, so it comes back as a standard
// sparse matrix, as requested).
//
// Expects U_obj to be an R list with elements:
//   U         - n x n sparse matrix (e.g. dgCMatrix), the Vecchia U factor
//   latent    - logical vector of length n
//   cond.yz   - single string ("zy" or something else)
//   ord.pred  - single string ("obspred" or something else)
// =====================================================================
// [[Rcpp::export]]
arma::sp_mat U2V_cpp(List U_obj) {

  arma::sp_mat U        = as_sparse_mat(U_obj["U"]);
  LogicalVector latentR = U_obj["latent"];
  std::string cond_yz   = as<std::string>(U_obj["cond.yz"]);
  std::string ord_pred  = as<std::string>(U_obj["ord.pred"]);

  arma::uvec latent     = as<arma::uvec>(latentR);        // 0/1 vector
  arma::uvec latent_idx = arma::find(latent == 1);

  // U.y <- U.obj$U[U.obj$latent, ]   (non-contiguous row selection)
  arma::sp_mat U_y = sp_select_rows(U, latent_idx);

  arma::sp_mat V_ord;

  if (cond_yz == "zy") {

    // -----------------------------------------------------------
    // V.ord <- revMat(U.y[, U.obj$latent, drop = FALSE])
    // -----------------------------------------------------------
    arma::sp_mat sub = sp_select_cols(U_y, latent_idx);
    V_ord = revMat_sp(sub);

  } else if (ord_pred != "obspred") {

    // -----------------------------------------------------------
    // W <- tcrossprod(U.y); W.rev <- revMat(W)
    // V.ord <- t(chol(W.rev))  (nearPD_cpp fallback if singular)
    // -----------------------------------------------------------
    arma::sp_mat W     = U_y * U_y.t();     // sparse x sparse -> sparse
    arma::sp_mat W_rev = revMat_sp(W);
    V_ord = chol_with_fallback_sparse(W_rev);

  } else {  // ord.pred == "obspred"

    arma::uword n = latent.n_elem;

    // last.obs <- max(which(!U.obj$latent))   (0-based here)
    arma::uvec obs_idx = arma::find(latent == 0);
    if (obs_idx.n_elem == 0) {
      Rcpp::stop("U2V_cpp: no observed locations found (all latent == TRUE).");
    }
    arma::uword last_obs0 = obs_idx.max();

    // latents.before <- sum(U.obj$latent[1:last.obs])
    arma::uword latents_before = arma::accu(latent.subvec(0, last_obs0));

    // latents.after <- sum(U.obj$latent[-(1:last.obs)])
    arma::uword latents_after = 0;
    if (last_obs0 + 1 < n) {
      latents_after = arma::accu(latent.subvec(last_obs0 + 1, n - 1));
    }

    // -----------------------------------------------------------
    // V.pr <- revMat(U.y[, (last.obs + 1):ncol(U.y), drop = FALSE])
    // contiguous column range -> native sp_mat slicing
    // -----------------------------------------------------------
    arma::sp_mat V_pr;
    if (last_obs0 + 1 <= n - 1) {
      V_pr = revMat_sp(U_y.cols(last_obs0 + 1, n - 1));
    } else {
      V_pr = arma::sp_mat(U_y.n_rows, 0);   // no prediction columns
    }

    // -----------------------------------------------------------
    // U.oo <- U.y[1:latents.before, 1:last.obs]   (contiguous)
    // A <- tcrossprod(U.oo); A.rev <- revMat(A)
    // V.oor <- t(chol(A.rev))  (nearPD_cpp fallback if singular)
    // -----------------------------------------------------------
    arma::sp_mat U_oo;
    if (latents_before > 0) {
      U_oo = U_y.submat(0, 0, latents_before - 1, last_obs0);
    } else {
      U_oo = arma::sp_mat(0, last_obs0 + 1);
    }

    arma::sp_mat A     = U_oo * U_oo.t();
    arma::sp_mat A_rev = revMat_sp(A);
    arma::sp_mat V_oor = chol_with_fallback_sparse(A_rev);

    // -----------------------------------------------------------
    // zeromat.sparse <- sparseMatrix(c(), c(), dims=c(latents.after, latents.before))
    // V.or <- rbind(zeromat.sparse, V.oor)
    // V.ord <- as(cbind(V.pr, V.or), "triangularMatrix")
    //   -> returned as a plain sparse matrix (dgCMatrix) instead,
    //      since Armadillo has no triangularMatrix class
    // -----------------------------------------------------------
    arma::sp_mat zeromat(latents_after, latents_before);   // all-zero, sparse
    arma::sp_mat V_or = arma::join_cols(zeromat, V_oor);   // rbind
    V_ord              = arma::join_rows(V_pr, V_or);      // cbind
  }

  return V_ord;   // -> dgCMatrix on the R side
}
