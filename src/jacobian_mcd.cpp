#include "scm.h"

// [[Rcpp::export(name="jacobian_mcd")]]
double jacobian_mcd(Rcpp::NumericMatrix& eta, Rcpp::NumericMatrix& res,  uint32_t& d,
                uint32_t& S_row, uint32_t& S_col,
                Rcpp::NumericVector rc_idx_s, Rcpp::NumericVector& rc_idx_t){
 using namespace Rcpp;
 uint32_t n = eta.rows();
 uint32_t i;
 uint32_t j;
 uint32_t k;
 int count;

 double aux_out1;
 double aux_out2;

 NumericVector etai(d);
 NumericMatrix aux_LD(d,d);
 NumericMatrix D(d,d);
 NumericMatrix L(d,d);

 for(i = 0; i < n; i++){
  etai = eta(i,_);
  aux_LD = mcd_LD(etai, d);
  L = aux_LD;
  for(j = 0; j < d; j++){
    D(j,j) = aux_LD(j,j);
    L(j,j) = 1.0;
    res(i, j + d) = L(S_row, j)*L(S_col, j)*pow(D(j,j),2);
  }
  count = 0;
  for(j = d; j < (d * (d + 1)/2); j++){
   aux_out1 = 0.0;
   aux_out2 = 0.0;
   for(k = 0; k < d; k++){
    aux_out1 = aux_out1 - L(rc_idx_t(count), k)*L(S_col, k)*pow(D(k,k),2);
    aux_out2 = aux_out2 - L(rc_idx_t(count), k)*L(S_row, k)*pow(D(k,k),2);
   }
   res(i, j + d) = L(S_row, rc_idx_s(count))*aux_out1 + L(S_col, rc_idx_s(count))*aux_out2;
   count += 1;
  }
 }
 return (1.0);
}
