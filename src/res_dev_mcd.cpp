#include <Rcpp.h>
using namespace Rcpp;

//' Residuals vector via mcd parameterisation
//'
//' @param eta Linear predictor (n x (d + dx(d+1)/2) matrix).
//' @param y Outcome (n x d matrix).
//' @param y res (n x d matrix).
//' @export

// [[Rcpp::export(name="res_dev_mcd")]]
double res_dev_mcd(NumericMatrix& eta, NumericMatrix& y, NumericMatrix& res){
 uint32_t n = y.rows();
 uint32_t d = y.cols();
 double aux_out = 0.0;

 uint32_t i;
 uint32_t j;
 uint32_t k;
 uint32_t count;

 for(i = 0; i < n; i++){
  res(i, 0) = exp(-0.5*eta(i,d))*(y(i,0)-eta(i,0));

  count = 2*d;
  for(j = d+1; j < 2*d; j++){
   for(k = d; k < j; k++){
    aux_out = aux_out + (y(i,k-d)-eta(i,k-d))*eta(i,count);
    count = count+1;
   }
   res(i, j-d) =  exp(-0.5*eta(i,j))*(aux_out+y(i,j-d)-eta(i,j-d));
   aux_out = 0.0;
  }
 }
 return 1.0;
}
