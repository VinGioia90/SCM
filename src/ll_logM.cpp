// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(name="ll_logm")]]
double ll_logm(arma::mat& eta, arma::mat& y) {
  using namespace arma;

  uint32_t n = y.n_rows;
  uint32_t d = y.n_cols;
  uint32_t c1 = 0;

  uint32_t i;
  uint32_t j;
  uint32_t k;

  double nll = 0.0;

  mat Theta(d, d, fill::zeros); //(negative) log Sigma matrix
  rowvec r(d, fill::zeros); // i-th residual vector

  for(i = 0; i < n; i++){
    // i-th residual
    r = y.row(i) - eta(i, span(0, d - 1));

    //(negative) log Sigma matrix
    c1 = 0;
    for(j = 1; j < d; j++){
      for(k = 0; k < j; k++){
        Theta(j, k) = -eta(i, c1 + 2 * d);
        Theta(k, j) = Theta(j, k);
        c1 += 1;
      }
    }
    Theta.diag() = -eta(i, span(d, 2 * d - 1));
    nll += 0.5 * as_scalar(r * expmat_sym(Theta) * r.t()) - 0.5 * trace(Theta);
  }

 return -nll;
}
