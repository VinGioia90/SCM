#include "scm.h"

// [[Rcpp::export(name="res_dev_logm")]]
double res_dev_logm(arma::mat& eta, arma::mat& y, arma::mat& resD){
  using namespace arma;
  uint32_t n = y.n_rows;
  uint32_t d = y.n_cols;

  uint32_t i;
  uint32_t j;

  rowvec r(d, fill::zeros);

  mat Sigma(d,d,fill::zeros);
  mat iSigma(d, d, arma::fill::zeros);
  mat C(d, d, arma::fill::zeros);
  rowvec etai(d, fill::zeros);

  for(i = 0; i < n; i++){
    etai = eta.row(i);
    Sigma = logM_Sigma(etai, d);
    iSigma = inv_sympd(Sigma);
    C = chol(iSigma);
    r = y.row(i) - eta(i, span(0, d - 1));
    resD.row(i) = r * C.t();
  }
  return(1.0);
}
