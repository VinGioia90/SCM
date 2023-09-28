#include "scm_der.h"

// [[Rcpp::export(name="d2_beta")]]
double d2_beta(arma::mat& X,  arma::mat& eta,  arma::mat& y,  Rcpp::List& jj,  uint32_t& K,
               arma::mat& lbb, arma::mat& l2, arma::vec& l2_v, arma::mat& l2_l, arma::vec& l2_v_l,
               arma::uvec& ig,  arma::vec& z, arma::vec& w, arma::mat& G, arma::vec& t,
               Rcpp::List&  b1_eta, Rcpp::List&  b1, Rcpp::List&  b2, Rcpp::List&  b3, arma::vec& ib1_eta, arma::vec& ib3,
               arma::vec& l2_el, arma::vec& l2_el2, uint32_t param){
  using namespace arma;
  uint32_t size_g = ig.n_elem - 1;

  uint32_t i;
  uint32_t j;
  uint32_t k;

  Rcpp::IntegerVector ijjj;
  Rcpp::IntegerVector ijjk;
  Rcpp::IntegerVector ik1;
  Rcpp::IntegerVector ik2;
  Rcpp::IntegerVector ik3;

  uint32_t ik1k;
  uint32_t ik2k;
  uint32_t ik3k;

  uint32_t ljjj;
  uint32_t ljjk;

  uint32_t igf;
  uint32_t igl;
  uint32_t ijjjf;
  uint32_t ijjjl;
  uint32_t ijjkf;
  uint32_t ijjkl;

  uint32_t c1;
  uint32_t c2;
  uint32_t c3;
  mat XXj;

  for(i = 0; i < size_g; i++){
    // Index for the observations block (igf: first row; igl: last row)
    igf = ig[i] + 1;
    igl = ig[i + 1];
    // set to zeros the cumulated second derivatives vector
    l2_v.zeros();
    l2_v_l.zeros();

    // Compute the 2nd derivatives w.r.t. eta
    if ( i < (size_g - 1) ) {
      if ( param == 1 ) d2_mcd_eta(eta.rows(igf, igl), y.rows(igf, igl), l2, l2_v, z, w, G, t, b1_eta, b3, ib1_eta, ib3);
      if ( param == 2 ) d2_logm_eta(eta.rows(igf, igl), y.rows(igf, igl), l2, l2_v, z, w, b1_eta, b3, ib1_eta, ib3);
    } else {
      if ( param == 1 ) d2_mcd_eta(eta.rows(igf, igl), y.rows(igf, igl), l2_l, l2_v_l, z, w, G, t, b1_eta, b3, ib1_eta, ib3);
      if ( param == 2 ) d2_logm_eta(eta.rows(igf, igl), y.rows(igf, igl), l2_l, l2_v_l, z, w, b1_eta, b3, ib1_eta, ib3);
    }

    // Counters
    c1 = 0;
    c2 = 0;
    c3 = 0;

    for(j = 0; j < K; j++){

      ijjj = jj[j];
      ik1 = b1[j];
      ik2 = b2[j];
      ik3 = b3[j];
      ljjj =  ijjj.length() - 1;
      ijjjf = ijjj[0];
      ijjjl = ijjj[ljjj];
      XXj = X(span(igf, igl), span(ijjjf, ijjjl));

      // Blocks full vs full
      if ( ik1.length() > 0 ) {
        for(k = 0; k < ik1.length(); k++){
          ik1k = ik1[k];
          ijjk = jj[ik1k];
          ljjk =  ijjk.length() - 1;
          ijjkf = ijjk[0];
          ijjkl = ijjk[ljjk];

          if ( i < (size_g - 1) ) {
            lbb(span(ijjjf, ijjjl), span(ijjkf, ijjkl)) += XXj.t() * (X(span(igf, igl), span(ijjkf, ijjkl)).each_col() % l2.col(l2_el[c1]));
          } else {
            lbb(span(ijjjf, ijjjl), span(ijjkf, ijjkl)) +=  XXj.t() * (X(span(igf, igl), span(ijjkf, ijjkl)).each_col() % l2_l.col(l2_el[c1]));
          }
          if ( (ik1k > j) && (i == (size_g - 1)) ) {
            lbb(span(ijjkf, ijjkl), span(ijjjf, ijjjl)) = lbb(span(ijjjf, ijjjl), span(ijjkf, ijjkl)).t();
          }
          c1 += 1;
        }
      }

      // Block full vs intercepts
      if ( (ik2.length() > 0) ) {
        for(k = 0; k < ik2.length(); k++){
          ik2k = ik2[k];
          ijjk = jj[ik2k];
          ljjk =  ijjk.length() - 1;
          ijjkf = ijjk[0];
          ijjkl = ijjk[ljjk];

          if ( i < (size_g - 1) ) {
            if ( ijjj.length() > 1 ) { // full vs intercepts
              lbb(span(ijjjf, ijjjl), span(ijjkf, ijjkl)) += XXj.t() * l2.col(l2_el2[c2]);
            } else { // intercepts vs full
              lbb(span(ijjjf, ijjjl), span(ijjkf, ijjkl)) += l2.col(l2_el2[c2]).t() * X(span(igf, igl), span(ijjkf, ijjkl));
            }
          } else {
            if ( ijjj.length() > 1 ) { // full vs intercepts
              lbb(span(ijjjf, ijjjl), span(ijjkf, ijjkl)) += XXj.t() * l2_l.col(l2_el2[c2]);
            } else { // intercepts vs full
              lbb(span(ijjjf, ijjjl), span(ijjkf, ijjkl)) += l2_l.col(l2_el2[c2]).t() * X(span(igf, igl), span(ijjkf, ijjkl));
            }
          }
          if ( (ik2k > j) && (i == (size_g - 1)) ) {
            lbb(span(ijjkf, ijjkl), span(ijjjf, ijjjl)) = lbb(span(ijjjf, ijjjl), span(ijjkf, ijjkl)).t();
          }
          c2 += 1;
        }
      }

      // Block intercepts vs intercepts
      if ( ik3.length() > 0 ) {
        for(k = 0; k < ik3.length(); k++){
          ik3k = ik3[k];
          ijjk = jj[ik3k];
          ijjkf = ijjk[0];

          if ( i < (size_g - 1) ) {
            lbb.at(ijjjf, ijjkf) += l2_v[c3];
          } else {
            lbb.at(ijjjf, ijjkf) += l2_v_l[c3];
          }

          if ( (ik3k > j) && (i == (size_g - 1)) ) {
            lbb.at(ijjkf, ijjjf) = lbb.at(ijjjf, ijjkf);
          }
          c3 += 1;
        }
      }
    }
  }
  return(1.0);
}

