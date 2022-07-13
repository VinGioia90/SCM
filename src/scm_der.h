#ifndef _SCM_H
#define _SCM_H



#include <RcppArmadillo.h>
using namespace Rcpp;

double d1_mcd_eta(NumericVector& eta, NumericVector& y, NumericVector& res, IntegerVector& z, IntegerVector& w, IntegerMatrix& G);

double d2_mcd_eta(const arma::mat& eta, const arma::mat& y, arma::mat& res,  Rcpp::List&  idx_jk, arma::vec& z, arma::vec& w, arma::mat& G, arma::vec& t);

#endif
