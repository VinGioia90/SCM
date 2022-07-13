#include <Rcpp.h>
using namespace Rcpp;

//' Score mcd
//'
//' @param eta Linear predictor (n x (d + dx(d+1)/2) matrix).
//' @param y Outcome (n x d matrix).
//' @param z idx1
//' @param w idx2
//' @param G idx matrix
//' @export

// [[Rcpp::export(name="d1_mcd_eta")]]
double d1_mcd_eta(NumericVector& eta, NumericVector& y, NumericVector& res, IntegerVector& z, IntegerVector& w, IntegerMatrix& G){
  uint32_t d = y.length();

  uint32_t j;
  uint32_t k;

  double aux_out = 0.0;
  NumericVector s(d);

    s(0) = y(0) - eta(0);
    for(j = 1; j < d; j++){
      for(k = 0; k < j; k++){
        aux_out = aux_out + (y(k) - eta(k)) * eta(G(j - 1,k));
      }
      s(j) = aux_out + y(j) - eta(j);
      aux_out = 0.0;
    }

    res(d-1) = exp(-eta(2*d-1)) * s(d-1);
    res(2*d-1) = -0.5 + 0.5*exp(-eta(2*d-1)) * s(d-1) * s(d-1);
    if(d > 2){
      res(3 * d - 1) = -exp(-eta(w(d-1) + d)) * s(w(d-1)) * (y(z(d-1)) - eta(z(d-1)));
    }

    for(j = 0; j < d - 1; j++){
      res(j) = exp(-eta(j + d)) * s(j);
      for(k = j + 1; k < d; k++){
        res(j) = res(j) + exp(-eta(k + d)) * s(k) * eta(G(k - 1,j));
      }
      res(j + d) = -0.5 + 0.5*exp(-eta(j + d)) * s(j) * s(j);
      res(j + 2 * d) = -exp(-eta(w(j) + d)) * s(w(j)) * (y(z(j)) - eta(z(j)));
    }

    for(j = d; j < d * (d - 1)/2; j++){
      res(j + 2 * d) = -exp(-eta(w(j) + d)) * s(w(j)) * (y(z(j)) - eta(z(j)));
    }

  return(1.0);
}
