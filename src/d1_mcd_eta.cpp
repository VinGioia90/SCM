// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(name="d1_mcd_eta")]]
double d1_mcd_eta(const arma::mat& eta, const arma::mat& y, arma::mat& d1l, arma::vec& z, arma::vec& w, arma::mat& Gm){
  using namespace arma;
  uint32_t d = y.n_cols;
  uint32_t n = y.n_rows;

  uint32_t i;
  uint32_t j;
  uint32_t k;
  double ee_s;
  double sj;
  uint32_t wj;

  double aux1 = 0.0;
  vec s(d, fill::zeros);
  rowvec r(d, fill::zeros);

  for(i = 0; i < n; i++){
    r = y.row(i) - eta(i, span(0, d - 1));

    // sum in brackets (computed one time and used below)
    s[0] = r[0];
    for(j = 1; j < d; j++){
      for(k = 0; k < j; k++){
        aux1 +=  r[k] * eta.at(i, Gm.at(j - 1, k));
      }
      s[j] = aux1 + r[j];
      aux1 = 0.0;
    }
    // Block 1 - mean; Block 2 - logD2; Block 3 - T;
    // Last element of the Block 1
    sj = s[d - 1];
    wj = w[d - 1];
    ee_s =  exp(-eta.at(i, 2 * d - 1)) * sj;
    d1l.at(i, d - 1) = ee_s;
    // Last element of the Block 2
    d1l.at(i, 2 * d - 1) = -0.5 + 0.5 * ee_s * sj;
    // One element of the Block 3
    if(d > 2){
      d1l.at(i, 3 * d - 1) = -exp(-eta.at(i, wj + d)) * s[wj] * r[z[d - 1]];
    }

    for(j = 0; j < d - 1; j++){
      sj = s[j];
      wj = w[j];
      ee_s =  exp(-eta.at(i, j + d)) * sj;
      // Elements of Block 1
      d1l.at(i, j) = ee_s;
      for(k = j + 1; k < d; k++){
        d1l.at(i, j) += exp(-eta.at(i, k + d)) * s[k] * eta.at(i, Gm.at(k - 1, j));
      }
      // Elements of Block 2
      d1l.at(i, j + d) = -0.5 + 0.5 * ee_s * sj;
      // Some elements of Block 3
      d1l.at(i, j + 2 * d) = -exp(-eta.at(i, wj + d)) * s[wj] * r[z[j]];
    }
    // Remaining elements of Block 3
    for(j = d; j < d * (d - 1)/2; j++){
      wj = w[j];
      d1l.at(i, j + 2 * d) = -exp(-eta.at(i, wj + d)) * s[wj] * r[z[j]];
    }
  }
  return(1.0);
}
