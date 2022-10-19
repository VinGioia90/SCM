// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(name="ll_mcd")]]
double ll_mcd(const arma::mat& eta, const arma::mat& y){
 using namespace arma;
 uint32_t d = y.n_cols;
 uint32_t n = y.n_rows;

 double out = 0.0;
 double aux1 = 0.0;

 double eij;
 uint32_t i;
 uint32_t j;
 uint32_t k;
 uint32_t c1;
 rowvec r(d, fill::zeros);

 for(i = 0; i < n; i++){
  r = y.row(i) - eta(i, span(0, d - 1));
  eij = eta.at(i, d);
  out +=  0.5 * eij + 0.5 * exp(-eij) * pow(r[0], 2);

  c1 = 2 * d;
  for(j = d + 1; j < 2 * d; j++){
   eij = eta.at(i, j);
   for(k = d; k < j; k++){
    aux1 += r[k - d] * eta.at(i, c1);
    c1 += 1;
   }
   out +=  0.5 * eij + 0.5 * exp(-eij) * pow(aux1 + r[j - d], 2);
   aux1 = 0.0;
  }
 }
 return -out;
}

