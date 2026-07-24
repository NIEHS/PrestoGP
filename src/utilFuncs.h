#ifndef UTIL_FUNCS_H
#define UTIL_FUNCS_H

#define ARMA_64BIT_WORD 1
#define ARMA_WARN_LEVEL 1
#ifdef _OPENMP
#include <omp.h>
#endif
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]
#include <cmath>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <RcppEigen.h>
#include <Eigen/SparseCholesky>
#include <vector>
#include <limits>

using namespace Rcpp;
using namespace arma;

double matern1(const double &x, const double &smooth, const double &alpha);
double dist1(const arma::rowvec &a, const arma::rowvec &b);
arma::vec na_omit_c(arma::vec x);
arma::mat create_param_sequence(const double P, const double ns);
arma::mat nearPD_wrap(arma::mat x);
arma::sp_mat U2V_cpp(List U_obj);
Rcpp::List createUMultivariate(Rcpp::List vec_approx, arma::vec params);

#endif
