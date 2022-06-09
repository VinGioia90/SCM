#include "scm.h"

// [[Rcpp::export(name="jacobian_mcd_rho")]]
double jacobian_mcd_rho(Rcpp::NumericMatrix& eta,  Rcpp::NumericMatrix& res, uint32_t& d,
                uint32_t& S_row, uint32_t& S_col,
                Rcpp::NumericVector rc_idx_s, Rcpp::NumericVector& rc_idx_t,  uint32_t& cor_flag){
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
 NumericVector aux_vec_res(d+d*(d+1)/2);
 NumericVector aux_vec_cov1(d+d*(d+1)/2);
 NumericVector aux_vec_cov2(d+d*(d+1)/2);
 NumericMatrix Sigma(d,d);


 for(i = 0; i < n; i++){
  etai = eta(i,_);
  aux_LD = mcd_LD(etai, d);
  L = aux_LD;
   for(j = 0; j < d; j++){
     D(j,j) = aux_LD(j,j);
     L(j,j) = 1.0;
   }
   Sigma = mcd_Sigma(L,D,d);
   
  
  for(j = 0; j < d; j++){
    if((cor_flag == 0) || (S_row == S_col)){
      res(i, j + d) = L(S_row, j)*L(S_col, j)*pow(D(j,j),2);
    }
    if((S_row != S_col) && (cor_flag == 1)){
      aux_vec_res(j) = 0.0;
      aux_vec_res(j + d) = L(S_row, j)*L(S_col, j)*pow(D(j,j),2);
      aux_vec_cov1(j) = 0.0;
      aux_vec_cov1(j + d) = L(S_row, j)*L(S_row, j)*pow(D(j,j),2);
      aux_vec_cov2(j) = 0.0;
      aux_vec_cov2(j + d) = L(S_col, j)*L(S_col, j)*pow(D(j,j),2);
      res(i, j + d) = pow(Sigma(S_row, S_row), -0.5) * aux_vec_res(j + d) * pow(Sigma(S_col, S_col), -0.5) - 
        0.5 * Sigma(S_col, S_row) * pow(Sigma(S_row,S_row), -1.5) *  pow(Sigma(S_col, S_col),-0.5) * aux_vec_cov1(j + d) - 
        0.5 * Sigma(S_col, S_row) * pow(Sigma(S_row, S_row), -0.5) *  pow(Sigma(S_col, S_col),-1.5) * aux_vec_cov2(j + d) ;

      }
  }

  count = 0;
  for(j = d; j < (d * (d + 1)/2); j++){
   aux_out1 = 0.0;
   aux_out2 = 0.0;
   for(k = 0; k < d; k++){
    aux_out1 = aux_out1 - L(rc_idx_t(count), k)*L(S_col, k)*pow(D(k,k),2);
    aux_out2 = aux_out2 - L(rc_idx_t(count), k)*L(S_row, k)*pow(D(k,k),2);
   }
   if((cor_flag == 0) || (S_row == S_col)){
     res(i, j + d) = L(S_row, rc_idx_s(count))*aux_out1 + L(S_col, rc_idx_s(count))*aux_out2;
   }
   if((S_row != S_col) && (cor_flag == 1)){
     aux_vec_res(j + d) = L(S_row, rc_idx_s(count))*aux_out1 + L(S_col, rc_idx_s(count))*aux_out2;
     
     aux_vec_cov1(j + d) = L(S_row, rc_idx_s(count))*aux_out2 + L(S_row, rc_idx_s(count))*aux_out2;
     aux_vec_cov2(j + d) = L(S_col, rc_idx_s(count))*aux_out1 + L(S_col, rc_idx_s(count))*aux_out1;
     
    res(i, j + d) = pow(Sigma(S_row, S_row), -0.5) * aux_vec_res(j + d) * pow(Sigma(S_col, S_col), -0.5) - 
    0.5 * Sigma(S_col, S_row) * pow(Sigma(S_row,S_row), -1.5) *  pow(Sigma(S_col, S_col),-0.5) * aux_vec_cov1(j + d) - 
    0.5 * Sigma(S_col, S_row) * pow(Sigma(S_row, S_row), -0.5) *  pow(Sigma(S_col, S_col),-1.5) * aux_vec_cov2(j + d) ;
   }
   count += 1;
  }
  
 }
 return (1.0);
}