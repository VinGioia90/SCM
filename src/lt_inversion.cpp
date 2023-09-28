// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(name="lt_inversion")]]
arma::mat lt_inversion(arma::mat& X){
 using namespace arma;

 uint32_t d = X.n_cols;//();
 uint32_t j;
 uint32_t k;
 uint32_t l;
 double s;
 mat out(d, d, fill::zeros);

 for(j = 0; j < d; j++){
  out(j, j) = 1.0;
  if(j > 0){
   for(k = 0; k < j; k++){
    s = 0.0;
    for(l = k; l < j; l++){
     s += X(j, l) * out(l, k);
    }
    out(j, k) = -s;
   }
  }
 }
 return(out);
}
