#include "scm.h"

//' Residuals vector via logm parameterisation
//'
//' @param eta Linear predictor (n x (d + dx(d+1)/2) matrix).
//' @param y Outcome (n x d matrix).
//' @param y res (n x d matrix).
//' @export

// [[Rcpp::export(name="res_dev_logm")]]
double res_dev_logm(Rcpp::NumericMatrix& eta, Rcpp::NumericMatrix& y, Rcpp::NumericMatrix& res){
  using namespace Rcpp;
  uint32_t n = y.rows();
  uint32_t d = y.cols();

  uint32_t i;
  uint32_t j;

  NumericVector ri(d);
  arma::vec ri_a(d,arma::fill::zeros);
  NumericVector etai(d+d*(d+1)/2);

  NumericMatrix Sigma(d, d);
  arma::mat iSigma(d, d, arma::fill::zeros);
  arma::mat isqrtSigma(d, d, arma::fill::zeros);
  arma::mat isqrtSigmaj(1,d, arma::fill::zeros);
  for(i = 0; i < n; i++){
    etai = eta(i,_);
    Sigma = logM_Sigma(etai, d);
    iSigma = arma::inv_sympd(as<arma::mat>(Sigma));
    isqrtSigma = arma::sqrtmat_sympd(iSigma);

    for(j = d; j < 2*d; j++){
      ri(j - d) = y(i, j-d) - eta(i,j-d);
    }
    ri_a = as<arma::vec>(ri);

    for(j = d; j < 2*d; j++){
      isqrtSigmaj =isqrtSigma.row(j-d);
      res(i, j-d) = as_scalar( isqrtSigmaj*ri_a);
    }

  }
  return(1.0);
}
