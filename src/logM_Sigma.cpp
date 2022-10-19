// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(name="logM_Sigma")]]
arma::mat logM_Sigma(arma::rowvec& x, uint32_t& d){
 using namespace arma;

 mat logSigma(d, d, arma::fill::zeros);
 mat Sigma(d, d, arma::fill::zeros);

 uint32_t j;
 uint32_t k;

 logSigma(0, 0) = x(d);
 uint32_t count = d;
 for(j = 1; j < d; j++){
  logSigma(j, j) = x(j + d);
  for(k = 0; k <  j; k++){
   logSigma(j, k) = logSigma(k, j) = x(count + d);
   count += 1;
  }
 }

 Sigma = expmat_sym(logSigma);
 return(Sigma);
}
