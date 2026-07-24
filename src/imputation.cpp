#include "utilFuncs.h"

arma::mat create_param_sequence(const double P, const double ns = 1);

// [[Rcpp::export]]
arma::mat MMatern_cov(const arma::mat &locs,
                      const arma::ivec &y_ndx,
                      const arma::vec &covparams,
                      const double P) {

  // Get parameter sequences (assumed to return a 5x2 integer matrix,
  // where rows correspond to sigma, rangep, smoothness, nuggets, rho
  // and columns correspond to start/end indices)
  arma::mat param_seq = create_param_sequence(P);

  // Extract parameters using 0-based indices from param_seq
  arma::vec sigma      = covparams.subvec(param_seq(0,0) - 1,
					  param_seq(0,1) - 1);
  arma::vec rangep     = covparams.subvec(param_seq(1,0) - 1,
					  param_seq(1,1) - 1);
  arma::vec smoothness = covparams.subvec(param_seq(2,0) - 1,
					  param_seq(2,1) - 1);
  arma::vec nuggets    = covparams.subvec(param_seq(3,0) - 1,
					  param_seq(3,1) - 1);
  arma::vec rho        = covparams.subvec(param_seq(4,0) - 1,
					  param_seq(4,1) - 1);

  // Build rho matrix (P x P): upper triangle filled with rho values,
  // then symmetrized, diagonal set to 1
  arma::mat rho_mat(P, P, arma::fill::zeros);
  int rho_idx = 0;
  for (int j = 1; j < P; j++) {       // column
    for (int i = 0; i < j; i++) {     // row (upper triangle, col-major)
      rho_mat(i, j) = rho(rho_idx++);
    }
  }
  rho_mat = rho_mat + rho_mat.t();
  rho_mat.diag().ones();

  int n = locs.n_rows;
  arma::mat Sigma_hat(n, n, arma::fill::zeros);

  for (int i = 0; i < P; i++) {
    for (int j = i; j < P; j++) {

      // Find row indices where y_ndx == i+1 and y_ndx == j+1 (1-based)
      arma::uvec which_i = arma::find(y_ndx == (i + 1));
      arma::uvec which_j = arma::find(y_ndx == (j + 1));

      if (which_i.n_elem == 0 || which_j.n_elem == 0) continue;

      // Compute cross-covariance parameters
      double smooth_ii = smoothness(i);
      double smooth_jj = smoothness(j);
      double smooth_ij = (smooth_ii + smooth_jj) / 2.0;

      double alpha_ii  = 1.0 / rangep(i);
      double alpha_jj  = 1.0 / rangep(j);
      double alpha_ij  = std::sqrt((alpha_ii * alpha_ii + alpha_jj * alpha_jj) / 2.0);

      // Scaling constant (mirrors the R formula)
      double scaling = rho_mat(i, j)
        * std::sqrt(sigma(i)) * std::sqrt(sigma(j))
        * std::pow(alpha_ii, smooth_ii)
        * std::pow(alpha_jj, smooth_jj)
        * boost::math::tgamma(smooth_ij)
        / (std::pow(alpha_ij, 2.0 * smooth_ij)
           * std::sqrt(boost::math::tgamma(smooth_ii)
                       * boost::math::tgamma(smooth_jj)));

      // Fill the (which_i x which_j) block element-by-element,
      // computing pairwise Euclidean distances on-the-fly (mimics rdist)
      int ni = which_i.n_elem;
      int nj = which_j.n_elem;

      for (int ii = 0; ii < ni; ii++) {
        for (int jj = 0; jj < nj; jj++) {
          double d = dist1(locs.row(which_i(ii)), locs.row(which_j(jj)));
          double cov_val = scaling * matern1(d, smooth_ij, alpha_ij);
          Sigma_hat(which_i(ii), which_j(jj)) = cov_val;
        }
      }

      if (i != j) {
        // Symmetrize off-diagonal block
        for (int ii = 0; ii < ni; ii++) {
          for (int jj = 0; jj < nj; jj++) {
            Sigma_hat(which_j(jj), which_i(ii)) =
              Sigma_hat(which_i(ii), which_j(jj));
          }
        }
      } else {
        // Add nugget to diagonal of the i==j block
        for (int ii = 0; ii < ni; ii++) {
          Sigma_hat(which_i(ii), which_i(ii)) += nuggets(i);
        }
      }
    }
  }

  return Sigma_hat;
}
