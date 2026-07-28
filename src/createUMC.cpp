#include "utilFuncs.h"

arma::mat create_param_sequence(const double P, const double ns = 1);

Rcpp::List createU_helper_mat(const arma::mat &olocs, const arma::vec &ondx,
			      const arma::mat &curqys,
			      const arma::mat &curqzs, const arma::mat &vijs,
			      const arma::mat &aijs,
			      const arma::mat &full_const,
			      const arma::vec &nugget, const arma::vec &sig2,
			      const arma::vec &U_beginning,
			      const int n_cores) {
  int n = arma::as_scalar(ondx.n_elem);
  int m = arma::as_scalar(curqys.n_rows);
  int n_inds = 2 * n * (m + 3);
  n_inds -= sum(arma::linspace(0, m, m + 1));
  // arma::mat feeder(3, n_inds);
  arma::umat ndx(2, n_inds);
  arma::vec vals(n_inds);
  // vals.rows(0, 6) = U_beginning;
  // feeder.cols(0, 6) = {{1, 1, 2, 1, 3, 3, 4},
  //                     {1, 2, 2, 3, 3, 4, 4}};
  ndx.cols(0, 6) = {{0, 0, 1, 0, 2, 2, 3}, {0, 1, 1, 2, 2, 3, 3}};
  vals(arma::span(0, 6)) = U_beginning;
  // feeder.cols(0,6) = {{0, 0, 1, 0, 2, 2, 3},
  //                     {0, 1, 1, 2, 2, 3, 3}, {1,1,1,1,1,1,1}};
  // feeder(2, arma::span(0,6)) = U_beginning.t();
  int ind = 7;
  // int ind = 0;
  if (n_cores > -1) {
    omp_set_num_threads(n_cores);
  }
  bool stable = true;
#ifdef _OPENMP
#pragma omp parallel for reduction(&& : stable)
#endif
  for (arma::uword i = 3; i < (ondx.n_elem + 1); i++) {
    double temppi = ondx(i - 1);
    arma::vec cy = na_omit_c(curqys.col(i - 1));
    arma::vec cz = (curqzs.col(i - 1));
    arma::vec cq = na_omit_c(arma::join_cols(cy, cz));
    arma::uword nq = cq.n_elem;
    arma::vec k1(nq);
    arma::mat k2(nq, nq);
    // Rcpp::Rcout << i << " ";
    // Rcpp::Rcout << cq;
    for (arma::uword j = 0; j < nq; j++) {
      double temppj = ondx(cq(j) - 1);
      k1(j) =
          full_const(temppi - 1, temppj - 1) *
          matern1(dist1(olocs.row(i - 1), olocs.row(cq(j) - 1)),
                  vijs(temppi - 1, temppj - 1), aijs(temppi - 1, temppj - 1));
      for (arma::uword k = j; k < nq; k++) {
        double temppk = ondx(cq(k) - 1);
        k2(j, k) =
            full_const(temppj - 1, temppk - 1) *
            matern1(dist1(olocs.row(cq(j) - 1), olocs.row(cq(k) - 1)),
                    vijs(temppj - 1, temppk - 1), aijs(temppj - 1, temppk - 1));
      }
      // // add nugget if neighbor is not latent
      if (j >= cy.n_elem) {
        k2(j, j) += nugget(temppj - 1);
      }
    }
    // make symmetric
    k2 = arma::symmatu(k2);
    // Rcpp::Rcout << k2;
    arma::vec bi(nq);
    bool status = arma::solve(bi, k2, k1);
    if (!status) {
      stable = false;
      k2 = nearPD_wrap(k2);
      status = arma::solve(bi, k2, k1);
    }
    double epsilon = 1e-9;
    while (!status) {
      epsilon *= 2;
      k2.diag() += epsilon;
      status = arma::solve(bi, k2, k1);
    }
    double ritemp = arma::as_scalar(sig2(temppi - 1) - bi.t() * k1);
    // Note: This section is causing crashes
    //    if (ritemp <= 0) {
    //      k2 = nearPD_wrap(k2);
    //      arma::solve(bi, k2, k1);
    //      ritemp = arma::as_scalar(sig2(temppi - 1) - bi.t() * k1);
    //    }
    epsilon = 1e-9;
    while (ritemp <= 0) {
      stable = false;
      epsilon *= 2;
      k2.diag() += epsilon;
      arma::solve(bi, k2, k1);
      ritemp = arma::as_scalar(sig2(temppi - 1) - bi.t() * k1);
      // Rcpp::Rcout << ritemp << "\n";
    }
    double ri = sqrt(ritemp);
    bi *= (-1.0 / ri);
    // arma::vec bi = arma::solve(k2, k1);
    // uhat(2*i - 1, 2*i - 1) = pow(nugget(temppi - 1), -0.5);
    // uhat(2*i - 2, 2*i - 1) = -1.0 * uhat(2*i - 1, 2*i - 1);
    int siz = cq.n_elem + 3;
    arma::umat ndx_private(2, siz);
    arma::vec vals_private(siz);
    // uhat(2*i - 2, 2*i - 2) = 1.0 / ri;
    for (arma::uword j = 0; j < cq.n_elem; j++) {
      arma::uword xind = 2 * cq(j) - 2;
      if (j >= cy.n_elem) {
        xind += 1;
      }
      // feeder.col(ind) = {(double)xind,(double) 2*i - 2, bi(j)};
      // ndx.col(ind) = {xind, 2 * i - 2};
      // vals(ind) = bi(j);
      ndx_private.col(j) = {xind, 2 * i - 2};
      vals_private(j) = bi(j);
      // vals[ind] = bi(j);
      // ind++;
      // uhat(xind, 2.0*i - 2.0) = bi(j);
    }
    // vals[ind] = 1.0 / ri;
    double feed_temp =
        pow(nugget(temppi - 1),
            -0.5); // used twice so temp storage to avoid recomputation
    // feeder.col(ind) = {(double)2*i - 2,(double)2*i - 2, 1.0 / ri};
    // feeder.col(ind + 2) = {(double)2*i - 1,(double) 2*i - 1, feed_temp};
    // ndx.col(ind) = {2 * i - 2, 2 * i - 2};
    // vals(ind) = 1.0 / ri;
    // ndx.col(ind + 2) = {2 * i - 1, 2 * i - 1};
    // vals(ind + 2) = feed_temp;
    ndx_private.col(nq) = {2 * i - 2, 2 * i - 2};
    vals_private(nq) = 1.0 / ri;
    ndx_private.col(nq + 2) = {2 * i - 1, 2 * i - 1};
    vals_private(nq + 2) = feed_temp;
    // vals[ind + 2] = pow(nugget(temppi - 1), -0.5);
    // vals[ind + 1] = -1.0 * vals[ind + 2];
    // feeder.col(ind + 1) = {(double)2*i - 2,(double) 2*i - 1, -1.0 *
    // feed_temp};
    // ndx.col(ind + 1) = {2 * i - 2, 2 * i - 1};
    // vals(ind + 1) = -1.0 * feed_temp;
    // ind += 3;
    ndx_private.col(nq + 1) = {2 * i - 2, 2 * i - 1};
    vals_private(nq + 1) = -1.0 * feed_temp;
    int local_ind;
 #pragma omp atomic capture
    {local_ind = ind; ind += siz;}
    ndx.cols(local_ind, local_ind + siz - 1) = ndx_private;
    vals.subvec(local_ind, local_ind + siz - 1) = vals_private;
  }
  arma::sp_mat uhat(ndx, vals);
  // 0 indexing -> 1 indexing for R
  // feeder.row(0) += 1.0;
  // feeder.row(1) += 1.0;
  // return feeder;
  // return uhat;
  return List::create(Named("U") = uhat, Named("stable") = stable);
}

// ---------------------------------------------------------------------
// Helper: convert a (possibly NULL) element of an R list of index vectors
// into an arma::vec. Used to rebuild q.list$q.y / q.list$q.z.
// ---------------------------------------------------------------------
//static arma::vec qelem_to_vec(SEXP s) {
//  if (s == R_NilValue) {
//    return arma::vec();
//  }
//  Rcpp::NumericVector nv = Rcpp::as<Rcpp::NumericVector>(s);
//  arma::vec out(nv.begin(), nv.size(), false, true);
//  return arma::vec(out); // copy, since nv may go out of scope
//}

//' Create the sparse triangular matrix U for multivariate Vecchia models
//'
//' This creates the sparse triangular matrix U for multivariate Vecchia
//' models. This matrix can be used to estimate the likelihood or transform
//' the data to be iid. This function is a multivariate version of
//' \code{\link[GPvecchia]{createU}}.
//'
//' @param vec_approx Object returned by \code{\link{vecchia_Mspecify}}.
//' @param params Vector of covariance parameters. See
//' \code{\link{create_param_sequence}} or the examples below for details
//' about the format of this vector.
//'
//' @return A list containing the sparse upper trianguler U, plus additional
//' objects required for other functions.
//'
//' @seealso \code{\link[GPvecchia]{createU}}, \code{\link{vecchia_Mspecify}},
//' \code{\link{create_param_sequence}}
//'
//' @references
//' \itemize{
//' \item Apanasovich, T.V., Genton, M.G. and Sun, Y. "A valid Matérn class of
//' cross-covariance functions for multivariate random fields with any number
//' of components", Journal of the American Statistical Association (2012)
//' 107(497):180-193.
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
//' soil.u <- createUMultivariate(soil.va, params)
// [[Rcpp::export]]
Rcpp::List createUMultivariate(Rcpp::List vec_approx, arma::vec params) {

  // -------------------- unpack vec.approx --------------------
  int P = Rcpp::as<int>(vec_approx["P"]);

  arma::mat qys = Rcpp::as<arma::mat>(vec_approx["qy.mat"]);
  arma::mat qzs = Rcpp::as<arma::mat>(vec_approx["qz.mat"]);

  //  Rcpp::List q_list = vec_approx["q.list"];
  //  Rcpp::List q_y = q_list["q.y"];
  //  Rcpp::List q_z = q_list["q.z"];

  arma::mat olocs = Rcpp::as<arma::mat>(vec_approx["locsord"]);
  int n = olocs.n_rows;

  // ondx may come over as integer or double; coerce defensively.
  arma::vec ondx = Rcpp::as<arma::vec>(vec_approx["ondx"]);

  int ncores = Rcpp::as<int>(vec_approx["n.cores"]);

  Rcpp::LogicalVector obs_lv = vec_approx["obs"];
  arma::uvec obs(obs_lv.size());
  for (int i = 0; i < obs_lv.size(); i++) obs(i) = obs_lv[i] ? 1u : 0u;

  std::string ord_pred = Rcpp::as<std::string>(vec_approx["ord.pred"]);

  // -------------------- parameter unpacking --------------------
  // 1-indexed [begin,end] rows: 1=sig2, 2=rangep, 3=smoothness, 4=nugget, 5=rho
  arma::mat param_seq = create_param_sequence(P);

  arma::vec sig2       = params.subvec(param_seq(0, 0) - 1, param_seq(0, 1) - 1);
  arma::vec rangep     = params.subvec(param_seq(1, 0) - 1, param_seq(1, 1) - 1);
  arma::vec smoothness = params.subvec(param_seq(2, 0) - 1, param_seq(2, 1) - 1);
  arma::vec nugget     = params.subvec(param_seq(3, 0) - 1, param_seq(3, 1) - 1);
  arma::vec rho        = params.subvec(param_seq(4, 0) - 1, param_seq(4, 1) - 1);

  // rho.mat: fill strict upper triangle in R's column-major upper.tri()
  // order (col j: rows 0..j-1), then symmetrize, then set diagonal to 1.
  arma::mat rho_mat(P, P, arma::fill::zeros);
  {
    int idx = 0;
    for (int j = 0; j < P; j++) {
      for (int i = 0; i < j; i++) {
        rho_mat(i, j) = rho(idx++);
      }
    }
  }
  rho_mat = rho_mat + rho_mat.t();
  rho_mat.diag().ones();

  // -------------------- uvec (7 special first-two-location entries) ---
  arma::vec uvec(7, arma::fill::zeros);

  int ondx1 = ondx(0) - 1; // R's ondx[1] -> 0-indexed
  int ondx2 = ondx(1) - 1; // R's ondx[2] -> 0-indexed

  uvec(0) = std::pow(sig2(ondx1), -0.5);
  uvec(2) = std::pow(nugget(ondx1), -0.5);
  uvec(1) = -1.0 * uvec(2);
  uvec(6) = std::pow(nugget(ondx2), -0.5);
  uvec(5) = -1.0 * uvec(6);

  double vii = smoothness(ondx1);
  double vjj = smoothness(ondx2);
  double vij = (vii + vjj) / 2.0;
  double aii = 1.0 / rangep(ondx1);
  double ajj = 1.0 / rangep(ondx2);
  double aij = std::sqrt((aii * aii + ajj * ajj) / 2.0);

  // R: Matern(rdist(olocs[1,], olocs[2,]), smoothness = vij, alpha = aij)
  double dist12 = dist1(olocs.row(0), olocs.row(1)); // rows 1,2 in R -> 0,1
  double mat12  = matern1(dist12, vij, aij);

  double K1 = rho_mat(ondx1, ondx2) * std::sqrt(sig2(ondx1)) * std::sqrt(sig2(ondx2)) *
              std::pow(aii, vii) * std::pow(ajj, vjj) *
              boost::math::tgamma(vij) /
              (std::pow(aij, 2.0 * vij) *
               std::sqrt(boost::math::tgamma(vii) * boost::math::tgamma(vjj))) *
              mat12;

  double K2 = sig2(ondx1);
  double bi = K1 / K2;
  double ri = sig2(ondx2) - bi * K1;

  uvec(4) = std::pow(ri, -0.5);
  uvec(3) = -1.0 * bi * std::pow(ri, -0.5);

  // -------------------- P x P outer-product-style matrices ------------
  arma::mat vijs(P, P), aijs(P, P), gammas(P, P), expprod(P, P), sigs(P, P);
  for (int i = 0; i < P; i++) {
    for (int j = 0; j < P; j++) {
      vijs(i, j) = (smoothness(i) + smoothness(j)) / 2.0;
      aijs(i, j) = std::sqrt((1.0 / (rangep(i) * rangep(i)) +
                              1.0 / (rangep(j) * rangep(j))) / 2.0);
      gammas(i, j) = boost::math::tgamma((smoothness(i) + smoothness(j)) / 2.0) /
                     std::sqrt(boost::math::tgamma(smoothness(i)) *
                               boost::math::tgamma(smoothness(j)));
      expprod(i, j) = std::pow(rangep(i), -smoothness(i)) *
                      std::pow(rangep(j), -smoothness(j));
      sigs(i, j) = std::sqrt(sig2(i)) * std::sqrt(sig2(j));
    }
  }

  arma::mat full_const(P, P);
  for (int i = 0; i < P; i++) {
    for (int j = 0; j < P; j++) {
      full_const(i, j) = sigs(i, j) * gammas(i, j) * expprod(i, j) * rho_mat(i, j) /
                         std::pow(aijs(i, j), 2.0 * vijs(i, j));
    }
  }

  // -------------------- rebuild cur.qys / cur.qzs as NaN-padded mats ---
  //  int m;
  //  {
  //    arma::vec last_y = qelem_to_vec(q_y[n - 1]);
  //    arma::vec last_z = qelem_to_vec(q_z[n - 1]);
  //    m = (int)last_y.n_elem + (int)last_z.n_elem;
  //  }

  //  arma::mat qys(m, n);
  //  arma::mat qzs(m, n);
  //  qys.fill(arma::datum::nan);
  //  qzs.fill(arma::datum::nan);

  //  for (int i = 0; i < n; i++) {
  //    arma::vec yv = qelem_to_vec(q_y[i]);
  //    if (yv.n_elem > 0) {
  //      qys.submat(0, i, yv.n_elem - 1, i) = yv;
  //    }
  //    arma::vec zv = qelem_to_vec(q_z[i]);
  //    if (zv.n_elem > 0) {
  //      qzs.submat(0, i, zv.n_elem - 1, i) = zv;
  //    }
  //  }

  // -------------------- call the (already-ported) U builder -----------
  Rcpp::List UL = createU_helper_mat(olocs, ondx, qys, qzs, vijs, aijs,
                                      full_const, nugget, sig2, uvec, ncores);
  arma::sp_mat U = UL["U"];

  // -------------------- drop prediction-only columns/rows -------------
  int nobs = (int)arma::accu(obs);
  int n_pred = n - nobs;

  if (n_pred > 0) {
    if (ord_pred != "obspred") {
      Rcpp::stop("Currently only obspred ordering is supported");
    }

    int ncolU = (int)U.n_cols;

    // Reproduce R: drop.seq <- seq(from = 2*nobs + 2, to = ncolU, by = 2)
    // (1-indexed); build a 0-indexed keep mask.
    arma::uvec keep_mask(ncolU, arma::fill::ones);
    for (int idx1 = 2 * nobs + 2; idx1 <= ncolU; idx1 += 2) {
      keep_mask(idx1 - 1) = 0;
    }
    arma::uvec keep_idx = arma::find(keep_mask);

    arma::ivec newpos(ncolU);
    newpos.fill(-1);
    for (arma::uword k = 0; k < keep_idx.n_elem; k++) {
      newpos(keep_idx(k)) = (int)k;
    }

    // Rebuild U restricted to kept rows/cols via triplet reconstruction
    // (arma::sp_mat has no direct arbitrary-index submat like R's `-drop.seq`).
    std::vector<arma::uword> rr, cc;
    std::vector<double> vv;
    rr.reserve(U.n_nonzero);
    cc.reserve(U.n_nonzero);
    vv.reserve(U.n_nonzero);

    for (arma::sp_mat::const_iterator it = U.begin(); it != U.end(); ++it) {
      arma::uword r = it.row();
      arma::uword c = it.col();
      if (newpos(r) >= 0 && newpos(c) >= 0) {
        rr.push_back((arma::uword)newpos(r));
        cc.push_back((arma::uword)newpos(c));
        vv.push_back(*it);
      }
    }

    arma::umat loc_mat(2, rr.size());
    for (size_t k = 0; k < rr.size(); k++) {
      loc_mat(0, k) = rr[k];
      loc_mat(1, k) = cc[k];
    }
    arma::vec vals(vv);

    arma::sp_mat U_new(loc_mat, vals, keep_idx.n_elem, keep_idx.n_elem);
    U = U_new;
  }

  // -------------------- latent indicator --------------------
  int n_final = (int)U.n_rows;
  Rcpp::LogicalVector latent(n_final, true);
  for (int idx1 = 2; idx1 <= 2 * nobs; idx1 += 2) {
    latent[idx1 - 1] = false;
  }

  // -------------------- assemble return list --------------------
  return Rcpp::List::create(
    Rcpp::Named("U") = U,
    Rcpp::Named("latent") = latent,
    Rcpp::Named("ord") = vec_approx["ord"],
    Rcpp::Named("obs") = vec_approx["obs"],
    Rcpp::Named("ord.pred") = vec_approx["ord.pred"],
    Rcpp::Named("ord.z") = vec_approx["ord.z"],
    Rcpp::Named("cond.yz") = vec_approx["cond.yz"],
    Rcpp::Named("ic0") = false,
    Rcpp::Named("stable") = UL["stable"]
  );
}
