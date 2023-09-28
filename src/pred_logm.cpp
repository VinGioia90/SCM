#include "scm.h"

// [[Rcpp::export(name="pred_logm")]]
double pred_logm(arma::mat& eta, arma::mat& pred,  uint32_t& d, uint32_t& cor_flag){
  using namespace arma;
  uint32_t n = eta.n_rows;
  uint32_t i;
  uint32_t j;
  uint32_t k;
  uint32_t count;

  arma::mat Sigma(d, d, fill::zeros);
  rowvec etai(d, fill::zeros);

  for(i = 0; i < n; i++){
    etai = eta.row(i);
    Sigma = logM_Sigma(etai, d);

    count = d;
    for(j = 0; j < d; j++){
      pred(i, j) = eta(i, j);
      pred(i, j + d) = Sigma(j, j);
      if(j > 0){
        for(k = 0; k < j; k++){
          if (cor_flag == 0) {
            pred(i, count + d) = Sigma(k, j);
          } else {
            pred(i, count + d) = Sigma(k, j)/(sqrt(Sigma(k, k)) * sqrt(Sigma(j, j)));
          }
          count = count + 1;
        }
      }
    }
  }
  return (1.0);
}
