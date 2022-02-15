// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Log-likelihood matrix logM
//'
//' @param eta Linear predictor (n x (d + dx(d+1)/2) matrix).
//' @param y Outcome (n x d matrix).
//' @export

// [[Rcpp::export(name="ll_logm")]]
double ll_logm(arma::mat& eta, arma::mat& y) {
   using namespace arma;
 uint32_t n = y.n_rows;
 uint32_t d = y.n_cols;
 uint32_t count = 0;

 uint32_t i;
 uint32_t j;
 uint32_t k;

 double nll = 0.0;

 mat nAi(d,d,fill::zeros);
 rowvec ri(d,fill::zeros);

  for(i = 0; i < n; i++){

   for(j = 0; j < d; j++){
    ri(j) = y(i,j) - eta(i,j);
    nAi(j,j) = -eta(i,j + d);
   }

   count = 0;
   for(j=1; j<d; j++){
    for(k = 0; k < j; k++){
     nAi(j,k)= -eta(i,count+2*d);
     nAi(k,j) = -eta(i,count+2*d);
     count += 1;
    }
   }
   nll += 0.5 * as_scalar(ri * expmat_sym(nAi) * ri.t()) - 0.5 * trace(nAi);
 }

 return -nll;
}
