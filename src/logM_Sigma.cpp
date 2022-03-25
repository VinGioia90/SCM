// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(name="logM_Sigma")]]
Rcpp::NumericMatrix logM_Sigma(Rcpp::NumericVector& x, uint32_t& d){
 using namespace Rcpp;

 arma::mat logSigma(d, d, arma::fill::zeros);
 arma::mat Sigma(d, d, arma::fill::zeros);

 uint32_t j;
 uint32_t k;

 arma::vec x2 = as<arma::vec>(x);


 logSigma(0,0) = x2(d);
 uint32_t count = d;

 for(j = 1; j < d; j++){
  logSigma(j,j) = x2(j+d);
  for(k = 0; k <  j; k++){
   logSigma(j,k) = logSigma(k,j) = x2(count + d);
   count += 1;
  }
 }

 Sigma = arma::expmat_sym(logSigma);


 return(wrap(Sigma));
}
