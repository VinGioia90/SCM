// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(name="mcd_Sigma")]]
arma::mat mcd_Sigma(arma::mat& L, arma::mat& D, uint32_t& d){
  using namespace arma;

  mat P(d, d);
  mat Sigma(d, d);

  uint32_t j;
  uint32_t k;
  uint32_t l;

  for(j = 0; j < d; j++){
    P(j, j) = D(j, j);
    if ( j > 0 ) {
      for(k = 0; k < j; k++){
        P(j, k) = D(k, k) * L(j, k);
      }
    }
  }

  for(j = 0; j < d; j++){
    for(k = 0; k < (j + 1); k++){
      for(l = 0; l < (j + 1); l++){
        Sigma(k, j) = Sigma(k, j) + P(j, l) * P(k, l);
      }
    }
  }
 return(Sigma);
}
