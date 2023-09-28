// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(name="d2_mcd_eta")]]
double d2_mcd_eta(const arma::mat& eta, const arma::mat& y, arma::mat& d2l, arma::vec& d2l_v,   arma::vec& z, arma::vec& w, arma::mat& G, arma::vec& t, Rcpp::List&  b1_eta, Rcpp::List&  b3, arma::vec& ib1_eta, arma::vec& ib3){
  using namespace arma;
  uint32_t n = y.n_rows;
  uint32_t d = y.n_cols;
  uint32_t l = eta.n_cols;
  double ee_s;
  double ee;
  double eee;
  double sj;
  uint32_t wj;
  uint32_t k1;
  double swj;
  uint32_t zj;
  double ee_r;

  double aux1 = 0.0;

  mat out(l, l, fill::zeros);
  vec s(d, fill::zeros);
  rowvec r(d, fill::zeros);
  Rcpp::IntegerVector ik1;

  uint32_t c1 = 0;
  uint32_t i;
  uint32_t j;
  uint32_t k;
  uint32_t q;

  vec::iterator it_b1;
  vec::iterator it_b1_end;
  vec::iterator it_b3;
  vec::iterator it_b3_end;

  for(i = 0; i < n; i++){
    r = y.row(i) - eta(i,span(0, d - 1));

    // sum in brackets (computed one time and used below)
    s[0] = r[0];
    for(j = 1; j < d; j++){
      for(k = 0; k < j; k++){
        aux1 += r[k] * eta.at(i, G.at(j - 1, k));
      }
      s[j] = aux1 + r[j];
      aux1 = 0.0;
    }

    for(j = 0; j < d; j++){
      sj = s[j];
      ee = -exp(-eta.at(i, j + d));
      ee_s =  ee * sj;

      // Block (1,1) - mean vs mean
      out.at(j, j) = ee;

      for(k = j + 1; k < d; k++){
        eee = exp(-eta.at(i, k + d)) * eta.at(i, G.at(k - 1, j));
        // Block (1,1) - mean vs mean
        out.at(j, k) = -eee;
        for(q = k + 1; q < d; q++){
          aux1 +=  -exp(-eta.at(i, q + d)) * eta.at(i, G.at(q - 1, j)) * eta(i, G.at(q - 1,k));
        }

        // Block (1,1) - mean vs mean
        out.at(j, k) += aux1;
        aux1 = 0.0;

        // Block (1,2) - mean vs logD2
        out.at(j, k + d) = -eee * s[k] ;
        //Block (1,1) - mean vs mean
        out.at(j, j) -=  eee * eta.at(i, G.at(k - 1, j));
      }


      // Block (1,2) - mean vs logD2
      out.at(j, j + d) = ee_s;

      // Block (2,2) - logD2 vs logD2
      out.at(j + d, j + d) = 0.5 * ee_s * sj;
    }

    for(k = 0; k < d * (d - 1)/2; k++){
      k1 = k + 2 * d;
      wj = w[k];
      swj = s[wj];
      zj = z[k];
      ee = exp(-eta.at(i, wj + d));
      ee_r = ee * r[zj];
      //Block (1,3) - mean vs T
      for(j = 0; j < wj; j++){
        out.at(j, k1) = ee_r * eta.at(i, G.at(wj - 1, j));
      }
      out.at(zj, k1) += ee * swj;
      out.at(wj, k1) = ee_r;

      //Block (2,3) - logD2 vs T
      out.at(wj + d, k1) = ee_r * swj;

      //Block (3,3) - T vs T
      for(j = 0; j < t[k]; j++){
        out.at(k1, k1 + j) = -ee_r * r(z[k+j]);
      }
    }

    //Saving the nonzero elements in a vector or a matrix according if they involve the intercepts or not
    it_b1     = ib1_eta.begin();
    it_b1_end = ib1_eta.end();
    it_b3     = ib3.begin();
    it_b3_end = ib3.end();

    // Full + partial
    c1 = 0;
    for(; it_b1 !=  it_b1_end; ++it_b1){
      ik1 = b1_eta[(*it_b1)];
      for(k = 0; k < ik1.length(); k++){
        d2l.at(i, c1) = out((*it_b1), ik1[k]);
        c1 += 1;
      }
    }

    // intercepts
    c1 = 0;
    for(; it_b3 !=  it_b3_end; ++it_b3){
      ik1 = b3[(*it_b3)];
      for(k = 0; k < ik1.length(); k++){
        d2l_v[c1] += out((*it_b3), ik1[k]);
        c1 += 1;
      }
    }
  }
  return(1.0);
}

