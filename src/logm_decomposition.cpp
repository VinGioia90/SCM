// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(name="logm_decomposition")]]
Rcpp::NumericMatrix logm_decomposition(Rcpp::NumericMatrix& X){
 using namespace Rcpp;
 
 arma::mat eigvec;
 arma::vec eigval;

 arma::mat res;
 
 arma::eig_sym(eigval,eigvec,as<arma::mat>(X));
 
 res = eigvec * diagmat(log(eigval))* eigvec.t(); 

 return(wrap(res));
}
