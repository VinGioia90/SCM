#ifndef _SCM_H
#define _SCM_H



#include <RcppArmadillo.h>
using namespace Rcpp;

// mcd parametrisation
double d1_mcd_eta(const arma::mat& eta, const arma::mat& y, arma::mat& d1l, arma::vec& z, arma::vec& w, arma::mat& G);
double d2_mcd_eta(const arma::mat& eta, const arma::mat& y, arma::mat& d2l, arma::vec& d2l_v,
                  arma::vec& z, arma::vec& w, arma::mat& G, arma::vec& t,
                  Rcpp::List&  b1_eta, Rcpp::List&  b3, arma::vec& ib1_eta, arma::vec& ib3);
double d3_mcd_eta (const arma::mat& eta, const arma::mat& y, arma::mat& res,
                   arma::vec& z, arma::vec& w, arma::mat& G, arma::vec& t,
                   Rcpp::IntegerVector& h, Rcpp::IntegerVector& h2, Rcpp::IntegerVector& h3,
                   Rcpp::List& idx3_1, Rcpp::List& idx3_2,  Rcpp::List& idx3_3,
                   Rcpp::List& idx3_4, Rcpp::List& idx3_5, Rcpp::List& idx3_6, Rcpp::List& idx_neq0);

// logm parametrisation
double  d1_logm_eta(const arma::mat& eta, const arma::mat& y, arma::mat& d1l, arma::vec& z, arma::vec& w);
double  d2_logm_eta(const arma::mat& eta, const arma::mat& y, arma::mat& d2l, arma::vec& d2l_v,
                    arma::vec& z, arma::vec& w, Rcpp::List&  b1_eta, Rcpp::List&  b3, arma::vec& ib1_eta, arma::vec& ib3);



#endif
