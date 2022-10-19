#include "scm_der.h"

// [[Rcpp::export(name="dHess_drho")]]
double dHess_drho(arma::mat& X, arma::mat& eta,  arma::mat& y,  Rcpp::List& jj,  uint32_t& K, arma::mat& l3,  arma::mat& l3_l,
                  arma::uvec& ig, arma::mat& d1b,  arma::mat& d1eta, arma::mat& d1eta_l, arma::vec& V, arma::vec& V_l,
                  Rcpp::List& d1H,
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
  Rcpp::IntegerVector csj = idx_jkq[4];

  Rcpp::IntegerVector idxl3j;
  Rcpp::IntegerVector idxqj;
  uint32_t idxjj;
  uint32_t idxjk;
  Rcpp::IntegerVector ijjj;
  Rcpp::IntegerVector ijjk;

  Rcpp::NumericMatrix A;

  uint32_t no_eta = eta.n_cols;
  uint32_t m = d1b.n_cols;
  uint32_t ljjj;
  uint32_t ijjjf;
  uint32_t ijjjl;
  uint32_t ljjk;
  uint32_t ijjkf;
  uint32_t ijjkl;
  uint32_t il3jq;
  uint32_t idxqq;
  uint32_t igf;
  uint32_t igl;

  uint32_t i;
  uint32_t l;
  uint32_t j;
  uint32_t q;
  uint32_t c3;


  for(l = 0; l < m; l++){
    A = Rcpp::as<Rcpp::NumericMatrix>(d1H[l]);
    mat Bb(A.begin(),A.rows(), A.cols(), false);
    Bb.zeros();
  }

  for(i = 0; i < size_g; i++){
    igf = ig(i) + 1;
    igl = ig(i + 1);

    if ( i < (size_g - 1) ) {
      d3_mcd_eta(eta.rows(igf, igl), y.rows(igf, igl), l3, z, w,  G,  t, h, h2, h3, idx3_1, idx3_2, idx3_3, idx3_4, idx3_5, idx3_6, idx_neq0);
    } else {
      d3_mcd_eta(eta.rows(igf, igl), y.rows(igf, igl), l3_l, z, w,  G,  t, h, h2, h3, idx3_1, idx3_2, idx3_3, idx3_4, idx3_5, idx3_6, idx_neq0);
    }


    for(l = 0; l < m; l++){
      A = Rcpp::as<Rcpp::NumericMatrix>(d1H[l]);
      mat Aa(A.begin(),A.rows(), A.cols(), false);
      //Aa.zeros();
      for (j = 0; j < no_eta; j++) {
        ijjj = jj[j];
        ljjj =  ijjj.length() - 1;
        ijjjf = ijjj(0);
        ijjjl = ijjj(ljjj);
        if ( i < (size_g - 1) ) {
         d1eta.col(j) = X(span(igf, igl), span(ijjjf,ijjjl)) * d1b(span(ijjjf,ijjjl), l);
        } else {
         d1eta_l.col(j) = X(span(igf, igl), span(ijjjf,ijjjl)) * d1b(span(ijjjf,ijjjl), l);
        }
      }

      c3 = 0;

      for(j = 0; j < idxl3.length(); j++){
        V.zeros();
        V_l.zeros();

        idxl3j = idxl3[j];
        idxqj = idxq[j];
        idxjk = idxk(j)-1;

        idxjj = idxj(j)-1;
        ijjj = jj[idxjj];
        ljjj =  ijjj.length() - 1;
        ijjjf = ijjj(0);
        ijjjl = ijjj(ljjj);

        //if(j == csj(c3)){
        //  idxjj = idxj(j)-1;
        //  ijjj = jj[idxjj];
        //  ljjj =  ijjj.length() - 1;
        //  ijjjf = ijjj(0);
        //  ijjjl = ijjj(ljjj);
        //   XX = X(span(igf, igl), span(ijjjf,ijjjl));
        //   c3 += 1;
        //}


        ijjk = jj[idxjk];
        ljjk =  ijjk.length() - 1;
        ijjkf = ijjk(0);
        ijjkl = ijjk(ljjk);

        for(q = 0; q < idxl3j.length(); q++){
          il3jq = idxl3j(q)-1;
          idxqq = idxqj(q)-1;
          if ( i < (size_g - 1) ) {
           V += l3.col(il3jq) % d1eta.col(idxqq);
          }  else {
           V_l += l3_l.col(il3jq) % d1eta_l.col(idxqq);
          }
        }
        //if ( i < (size_g - 1) ) {
        //  Aa(span(ijjjf, ijjjl), span(ijjkf, ijjkl)) += XX.t() * (X(span(igf, igl), span(ijjkf,ijjkl)).each_col() % V);
        //} else {
        //  Aa(span(ijjjf, ijjjl), span(ijjkf, ijjkl)) += XX.t() * (X(span(igf, igl), span(ijjkf,ijjkl)).each_col() % V_l);
        //}
        //if((i == (size_g - 1)) && (idxjk > idxjj)) Aa(span(ijjkf, ijjkl), span(ijjjf, ijjjl)) = Aa(span(ijjjf, ijjjl), span(ijjkf, ijjkl)).t();

        if ( i < (size_g - 1) ) {
          Aa(span(ijjjf, ijjjl), span(ijjkf, ijjkl)) += X(span(igf, igl), span(ijjjf,ijjjl)).t() * (X(span(igf, igl), span(ijjkf,ijjkl)).each_col() % V);
        } else {
          Aa(span(ijjjf, ijjjl), span(ijjkf, ijjkl)) += X(span(igf, igl), span(ijjjf,ijjjl)).t() * (X(span(igf, igl), span(ijjkf,ijjkl)).each_col() % V_l);
        }
        if((i == (size_g - 1)) && (idxjk > idxjj)) Aa(span(ijjkf, ijjkl), span(ijjjf, ijjjl)) = Aa(span(ijjjf, ijjjl), span(ijjkf, ijjkl)).t();
      }
    }
  }
  return(1.0);
}

