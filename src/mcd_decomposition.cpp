#include "scm.h"

// [[Rcpp::export(name="mcd_decomposition")]]
arma::mat mcd_decomposition(arma::mat& X){
  using namespace arma;
  uint32_t d = X.n_cols;//cols();
  uint32_t j;
  uint32_t k;

  mat C = chol(X);
  vec D = C.diag();
  mat L(d, d, fill::zeros);
  mat out(d, d, fill::zeros);

  for(j = 0; j < d; j++){
    D(j) = C(j, j);
    L(j, j) = 1.0;
    if(j > 0){
      for(k = 0; k < j; k++){
        L(j, k) = C(k, j)/D(k);
      }
    }
  }

  out = lt_inversion(L);
  for(j = 0; j < d; j++){
    out(j, j) = log(pow(D(j), 2));
    if(j > 0){
      for(k = 0; k < j; k++){
        out(k, j)= out(j, k);
      }
    }
  }
  return(out);
}
