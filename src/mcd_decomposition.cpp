#include "scm.h"

// [[Rcpp::export(name="mcd_decomposition")]]
Rcpp::NumericMatrix mcd_decomposition(Rcpp::NumericMatrix& X){
 using namespace Rcpp;
 uint32_t d = X.cols();
 uint32_t j;
 uint32_t k;

 arma::mat C = arma::chol(as<arma::mat>(X));
 arma::vec D = C.diag();

 NumericMatrix L(d,d);
 NumericMatrix res(d,d);

 for(j = 0; j < d; j++){
  D(j) = C(j,j);
  L(j,j) = 1.0;
  if(j > 0){
   for(k = 0; k < j; k++){
    L(j,k) = C(k,j)/D(k);
   }
  }
 }

 res = lt_inversion(L);
 for(j = 0; j < d; j++){
  res(j,j) = log(pow(D(j),2));
  if(j > 0){
   for(k = 0; k < j; k++){
    res(k,j)= res(j,k);
   }
  }
 }
 return(res);
}
