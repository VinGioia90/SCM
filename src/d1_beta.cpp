#include "scm_der.h"

// [[Rcpp::export(name="d1_beta")]]
double d1_beta(arma::mat& X, arma::mat& eta,  arma::mat& y, Rcpp::List& jj, uint32_t& K, arma::vec& lb, arma::mat& l1, arma::mat& l1_l,  arma::uvec& ig, arma::vec& z, arma::vec& w, arma::mat& Gm, uint32_t& param){
  using namespace arma;
  uint32_t size_g = ig.n_elem - 1;
  uint32_t igf;
  uint32_t igl;

  uint32_t i;
  uint32_t j;

  Rcpp::IntegerVector ijjj;
  uint32_t l_jjj;
  uint32_t ijjjf;
  uint32_t ijjjl;


  for(i = 0; i < size_g; i++){
    igf = ig[i] + 1;
    igl = ig[i + 1];

    if ( i < (size_g - 1) ) {
      //param=1 means mcd
      //param=2 means logm
        if ( param == 1 )  d1_mcd_eta(eta.rows(igf, igl), y.rows(igf, igl), l1, z, w, Gm);
        if ( param == 2 )  d1_logm_eta(eta.rows(igf, igl), y.rows(igf, igl), l1, z, w);
    } else {
        if ( param == 1 )  d1_mcd_eta(eta.rows(igf, igl), y.rows(igf, igl), l1_l, z, w, Gm);
        if ( param == 2 )  d1_logm_eta(eta.rows(igf, igl), y.rows(igf, igl), l1_l, z, w);
    }


      for(j = 0; j < K; j++) {
        ijjj = jj[j];
        l_jjj = ijjj.length() - 1;
        ijjjf = ijjj[0];
        ijjjl = ijjj[l_jjj];
        if ( i < (size_g - 1) ) {
          lb.subvec(ijjjf, ijjjl) += X(span(igf, igl), span(ijjjf, ijjjl)).t() * l1.col(j);
        } else {
          lb.subvec(ijjjf, ijjjl) += X(span(igf, igl), span(ijjjf, ijjjl)).t() * l1_l.col(j);
        }
      }
  }
  return(1.0);
}
