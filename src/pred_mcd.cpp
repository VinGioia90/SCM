#include "scm.h"

// [[Rcpp::export(name = "pred_mcd")]]
double pred_mcd(arma::mat& eta, arma::mat& pred,  uint32_t& d, uint32_t& cor_flag){
  using namespace arma;
  uint32_t n = eta.n_rows;
  uint32_t i;
  uint32_t j;
  uint32_t k;
  uint32_t count;

  //NumericMatrix aux_LD(d, d);
  //NumericMatrix L(d, d);
  //NumericMatrix D(d, d);
  //NumericMatrix Sigma(d, d);
  //NumericVector etai(d);

  mat aux_LD(d, d, fill::zeros);
  mat L(d, d, fill::zeros);
  mat D(d, d, fill::zeros);
  mat Sigma(d, d, fill::zeros);
  rowvec etai(d, fill::zeros);


  for(i = 0; i < n; i++){
    //etai = eta(i,_);
    etai = eta.row(i);
    aux_LD = mcd_LD(etai, d);
    L = aux_LD;
    for(j = 0; j < d; j++){
      D(j, j) = aux_LD(j, j);
      L(j, j) = 1.0;
    }
    Sigma = mcd_Sigma(L, D, d);

    count = d;
    for(j = 0; j < d; j++){
      pred(i, j) = eta(i, j);
      pred(i, j + d) = Sigma(j, j);
      if ( j > 0 ) {
        for(k = 0; k < j; k++){
          if ( cor_flag == 0 ) {
            pred(i, count + d) = Sigma(k, j);
          } else {
            pred(i, count + d) = Sigma(k, j)/(sqrt(Sigma(k, k)) * sqrt(Sigma(j, j)));
          }
          count += 1;
        }
      }
    }
  }
  return (1.0);
}
