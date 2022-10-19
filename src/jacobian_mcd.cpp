#include "scm.h"

// [[Rcpp::export(name="jacobian_mcd")]]
double jacobian_mcd(arma::mat& eta,  arma::mat& res, uint32_t& d,
                    uint32_t& S_r, uint32_t& S_c,
                    arma::vec rc_idx_s, arma::vec& rc_idx_t,  uint32_t& cor_flag){
  using namespace arma;
  uint32_t n = eta.n_rows;
  uint32_t no_eta = eta.n_cols;
  uint32_t i;
  uint32_t j;
  uint32_t k;
  uint32_t count;

  double aux_out1;
  double aux_out2;

  //NumericVector etai(d);
  //NumericMatrix aux_LD(d, d);
  //NumericMatrix D(d, d);
  //NumericMatrix L(d, d);
  //NumericVector aux_vec_res(no_eta);
  //NumericVector aux_vec_cov1(no_eta);
  //NumericVector aux_vec_cov2(no_eta);
  //NumericMatrix Sigma(d, d);

  rowvec etai(d);
  mat aux_LD(d, d);
  mat D(d, d);
  mat L(d, d);
  rowvec aux_vec_res(no_eta);
  rowvec aux_vec_cov1(no_eta);
  rowvec aux_vec_cov2(no_eta);
  mat Sigma(d, d);


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

    for(j = 0; j < d; j++){
      if ( (cor_flag == 0) || (S_r == S_c) ) {
        res(i, j + d) = L(S_r, j) * L(S_c, j) * pow(D(j, j), 2);
      }
      if ( (S_r != S_c) && (cor_flag == 1) ) {
        aux_vec_res(j) = 0.0;
        aux_vec_res(j + d) = L(S_r, j) * L(S_c, j) * pow(D(j, j), 2);
        aux_vec_cov1(j) = 0.0;
        aux_vec_cov1(j + d) = L(S_r, j) * L(S_r, j) * pow(D(j, j),2);
        aux_vec_cov2(j) = 0.0;
        aux_vec_cov2(j + d) = L(S_c, j) * L(S_c, j) * pow(D(j, j), 2);
        res(i, j + d) = pow(Sigma(S_r, S_r), -0.5) * aux_vec_res(j + d) * pow(Sigma(S_c, S_c), -0.5) -
                            0.5 * Sigma(S_c, S_r) * pow(Sigma(S_r, S_r), -1.5) *  pow(Sigma(S_c, S_c), -0.5) * aux_vec_cov1(j + d) -
                            0.5 * Sigma(S_c, S_r) * pow(Sigma(S_r, S_r), -0.5) *  pow(Sigma(S_c, S_c), -1.5) * aux_vec_cov2(j + d);
      }
    }

    count = 0;
    for(j = d; j < (d * (d + 1)/2); j++){
      aux_out1 = 0.0;
      aux_out2 = 0.0;
      for(k = 0; k < d; k++){
        aux_out1 = aux_out1 - L(rc_idx_t(count), k) * L(S_c, k) * pow(D(k, k), 2);
        aux_out2 = aux_out2 - L(rc_idx_t(count), k) * L(S_r, k) * pow(D(k, k), 2);
      }
      if ( (cor_flag == 0) || (S_r == S_c) ) {
        res(i, j + d) = L(S_r, rc_idx_s(count)) * aux_out1 + L(S_c, rc_idx_s(count)) * aux_out2;
      }
      if ( (S_r != S_c) && (cor_flag == 1) ) {
        aux_vec_res(j + d) = L(S_r, rc_idx_s(count)) * aux_out1 + L(S_c, rc_idx_s(count)) * aux_out2;
        aux_vec_cov1(j + d) = L(S_r, rc_idx_s(count)) * aux_out2 + L(S_r, rc_idx_s(count)) * aux_out2;
        aux_vec_cov2(j + d) = L(S_c, rc_idx_s(count)) * aux_out1 + L(S_c, rc_idx_s(count)) * aux_out1;

        res(i, j + d) = pow(Sigma(S_r, S_r), -0.5) * aux_vec_res(j + d) * pow(Sigma(S_c, S_c), -0.5) -
                            0.5 * Sigma(S_c, S_r) * pow(Sigma(S_r, S_r), -1.5) *  pow(Sigma(S_c, S_c), -0.5) * aux_vec_cov1(j + d) -
                            0.5 * Sigma(S_c, S_r) * pow(Sigma(S_r, S_r), -0.5) *  pow(Sigma(S_c, S_c), -1.5) * aux_vec_cov2(j + d);
      }
      count += 1;
    }
  }
 return (1.0);
}
