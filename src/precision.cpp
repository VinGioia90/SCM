// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(name="precision")]]
Rcpp::NumericMatrix precision(Rcpp::NumericMatrix& X){
 using namespace Rcpp;
 arma::mat res;
 res = arma::inv_sympd(as<arma::mat>(X));

 return(wrap(res));
}
