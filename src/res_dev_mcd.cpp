// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(name="res_dev_mcd")]]
double res_dev_mcd(arma::mat eta, arma::mat& y, arma::mat& res){
 using namespace arma;
 uint32_t n = y.n_rows;
 uint32_t d = y.n_cols;
 double aux_out = 0.0;

 uint32_t i;
 uint32_t j;
 uint32_t k;
 uint32_t count;

 for(i = 0; i < n; i++){
  res(i, 0) = exp(-0.5 * eta(i, d)) * (y(i, 0) - eta(i, 0));

  count = 2 * d;
  for(j = d + 1; j < 2 * d; j++){
   for(k = d; k < j; k++){
    aux_out = aux_out + (y(i, k - d) - eta(i, k - d)) * eta(i, count);
    count += 1;
   }
   res(i, j - d) =  exp(-0.5 * eta(i, j)) * (aux_out + y(i, j - d) - eta(i, j - d));
   aux_out = 0.0;
  }
 }
 return(1.0);
}
