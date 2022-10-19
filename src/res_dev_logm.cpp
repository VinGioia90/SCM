#include "scm.h"

// [[Rcpp::export(name="res_dev_logm")]]
double res_dev_logm(arma::mat& eta, arma::mat& y, arma::mat& resD){
  using namespace arma;
  uint32_t n = y.n_rows;
  uint32_t d = y.n_cols;

  uint32_t i;
  uint32_t j;

  rowvec r(d, fill::zeros);

  mat Sigma(d,d,fill::zeros);//NumericMatrix Sigma(d, d);
  mat iSigma(d, d, arma::fill::zeros);
  mat isqrtSigma(d, d, arma::fill::zeros);
  mat isqrtSigmaj(1, d, arma::fill::zeros);
  rowvec etai(d, fill::zeros);

  for(i = 0; i < n; i++){
    etai = eta.row(i);
    Sigma = logM_Sigma(etai, d);
    iSigma = inv_sympd(Sigma);
    isqrtSigma = sqrtmat_sympd(iSigma);
    r = y.row(i) - eta(i, span(0, d - 1));
    for(j = 0; j < d; j++){
      isqrtSigmaj = isqrtSigma.row(j);
      resD(i, j) = as_scalar(isqrtSigmaj * r.t());
    }
  }
  return(1.0);
}
