// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(name="logm_decomposition")]]
arma::mat logm_decomposition(arma::mat& X){
 using namespace arma;

 mat eigvec;
 vec eigval;
 mat out;

 eig_sym(eigval, eigvec, X);
 out = eigvec * diagmat(log(eigval)) * eigvec.t();

 return(out);
}
