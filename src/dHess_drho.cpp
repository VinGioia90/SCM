#include "scm_der.h"

// [[Rcpp::export(name="dHess_drho")]]
double dHess_drho(arma::mat& X, arma::mat& eta,  arma::mat& y,  Rcpp::List& jj,  uint32_t& K, arma::mat& l3,  arma::mat& l3_l,
                       arma::uvec& ig, arma::mat& d1b,
                       arma::vec& a, arma::vec& a_l, arma::vec& d1H, arma::mat& fh,
                       arma::vec& z, arma::vec& w, arma::mat& G, arma::vec& t, Rcpp::List& idx_l3,  Rcpp::List& idx_neq0, Rcpp::List& idx_jkq){
  using namespace arma;
  uint32_t size_g = ig.n_elem - 1;

  Rcpp::IntegerVector h = idx_l3[0];
  Rcpp::IntegerVector h2 = idx_l3[1];
  Rcpp::IntegerVector h3 = idx_l3[2];
  Rcpp::List idx3_1 = idx_l3[3];
  Rcpp::List idx3_2 = idx_l3[4];
  Rcpp::List idx3_3 = idx_l3[5];
  Rcpp::List idx3_4 = idx_l3[6];
  Rcpp::List idx3_5 = idx_l3[7];
  Rcpp::List idx3_6 = idx_l3[8];

  Rcpp::List idxl3 = idx_jkq[0];
  Rcpp::List idxq = idx_jkq[1];
  Rcpp::IntegerVector idxj = idx_jkq[2];
  Rcpp::IntegerVector idxk = idx_jkq[3];

  Rcpp::IntegerVector idxl3j;
  Rcpp::IntegerVector idxqj;
  uint32_t idxjj;
  uint32_t idxjk;
  Rcpp::IntegerVector ijjj;
  Rcpp::IntegerVector ijjk;
  Rcpp::IntegerVector ijjq;

  uint32_t ljjj;
  uint32_t ijjjf;
  uint32_t ijjjl;
  uint32_t ljjk;
  uint32_t ijjkf;
  uint32_t ijjkl;
  uint32_t ljjq;
  uint32_t ijjqf;
  uint32_t ijjql;
  uint32_t il3jq;
  uint32_t idxqq;
  uint32_t igf;
  uint32_t igl;
  uint32_t mult;

  uint32_t i;
  uint32_t j;
  uint32_t q;

  d1H.zeros();

  for(i = 0; i < size_g; i++){
    igf = ig[i] + 1;
    igl = ig[i + 1];

    if ( i < (size_g - 1) ) {
      d3_mcd_eta(eta.rows(igf, igl), y.rows(igf, igl), l3, z, w,  G,  t, h, h2, h3, idx3_1, idx3_2, idx3_3, idx3_4, idx3_5, idx3_6, idx_neq0);
    } else {
      d3_mcd_eta(eta.rows(igf, igl), y.rows(igf, igl), l3_l, z, w,  G,  t, h, h2, h3, idx3_1, idx3_2, idx3_3, idx3_4, idx3_5, idx3_6, idx_neq0);
    }

    for(j = 0; j < idxl3.length(); j++){
      idxl3j = idxl3[j];
      idxqj = idxq[j];

      idxjk = idxk[j] - 1;
      idxjj = idxj[j] - 1;

      ijjj = jj[idxjj];
      ljjj =  ijjj.length() - 1;
      ijjjf = ijjj[0];
      ijjjl = ijjj[ljjj];

      ijjk = jj[idxjk];
      ljjk =  ijjk.length() - 1;
      ijjkf = ijjk[0];
      ijjkl = ijjk[ljjk];


      if ( i < (size_g - 1) ) {
        a = sum(((X(span(igf, igl), span(ijjjf,ijjjl)) * fh( span(ijjjf,ijjjl),  span(ijjkf,ijjkl))) % X(span(igf, igl), span(ijjkf,ijjkl))),1);
      } else {
        a_l = sum(((X(span(igf, igl), span(ijjjf,ijjjl)) * fh( span(ijjjf,ijjjl),  span(ijjkf,ijjkl))) % X(span(igf, igl), span(ijjkf,ijjkl))),1);
      }

      if ( idxj[j] == idxk[j] ){
        mult = 1;
      } else {
        mult = 2;
      }

      for(q = 0; q < idxl3j.length(); q++){
        il3jq = idxl3j[q] - 1;

        idxqq = idxqj[q] - 1;
        ijjq = jj[idxqq];
        ljjq = ijjq.length() - 1;
        ijjqf = ijjq[0];
        ijjql = ijjq[ljjq];

        if ( i < (size_g - 1) ) {
          d1H += mult * d1b.rows(span(ijjqf,ijjql)).t()*(X(span(igf, igl), span(ijjqf,ijjql)).t()  * (a % l3.col(il3jq)));
        }
        else {
          d1H += mult * d1b.rows(span(ijjqf,ijjql)).t()*(X(span(igf, igl), span(ijjqf,ijjql)).t()  * (a_l % l3_l.col(il3jq)));
        }
      }
    }
  }
  return(1.0);
}
