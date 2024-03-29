#ifndef _SCM_H
#define _SCM_H

#include <RcppArmadillo.h>
using namespace Rcpp;
//NumericMatrix lt_inversion(NumericMatrix& X);
arma::mat lt_inversion(arma::mat& X);

//NumericMatrix mcd_LD(NumericVector& x, uint32_t& d);
arma::mat mcd_LD(arma::rowvec& x, uint32_t& d);

//NumericMatrix mcd_Sigma(NumericMatrix& L, NumericMatrix& D,  uint32_t& d);
arma::mat mcd_Sigma(arma::mat& L, arma::mat& D,  uint32_t& d);

arma::mat logM_Sigma(arma::rowvec& x,   uint32_t& d);
#endif
