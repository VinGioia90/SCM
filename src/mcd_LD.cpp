#include "scm.h"

// [[Rcpp::export(name="mcd_LD")]]
arma::mat mcd_LD(arma::rowvec& x, uint32_t& d){
  using namespace arma;

  //NumericMatrix T(d, d);
  //NumericMatrix Y(d, d);
  mat T(d, d, fill::zeros);
  mat Y(d, d, fill::zeros);

  uint32_t j;
  uint32_t k;
  uint32_t count = d;
  T(0,0) = 1.0;
  for(j = 1; j < d; j++){
    T(j, j) = 1.0;
    for(k = 0; k < j; k++){
      T(j, k) = x(count + d);
      count += 1;
    }
  }

  Y = lt_inversion(T);
  for(j = 0; j < d; j++){
    Y(j, j) = exp(0.5 * x(j + d));
   }
  return(Y);
}
