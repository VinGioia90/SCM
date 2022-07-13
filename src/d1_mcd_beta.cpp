#include "scm_der.h"

//' Score mcd in beta
//'
//' @param X model matrix
//' @param jj List of lpi index
//' @param K number of lpi
//' @param lb matrix of the derivatives (initialization)
//' @param eta lpi
//' @param y outcome matrix
//'
//' @export

// [[Rcpp::export(name="d1_mcd_beta")]]
double d1_mcd_beta(NumericMatrix& X, List& jj, uint32_t& K, NumericVector& lb, NumericMatrix& eta, NumericMatrix& y, IntegerVector& z, IntegerVector& w, IntegerMatrix& G){
  uint32_t n = y.rows();
  uint32_t d = y.cols();

  uint32_t i;
  uint32_t j;
  uint32_t k;

  NumericVector resi(K); //gradiente in etai
  NumericVector etai(K);
  NumericVector yi(d);


  for(i = 0; i < n; i++){
    etai = eta(i,_);
    yi = y(i,_);
    d1_mcd_eta(etai,yi, resi, z, w, G);

    for(j = 0; j < K; j++) {
     NumericVector idx_jj = jj[j];
     uint32_t l_jj = idx_jj.length();

     for(k = 0; k < l_jj; k++){
         lb(idx_jj(k)-1) += resi(j)*X(i,idx_jj(k)-1);
     }
    }

  }

  return(1.0);
}
