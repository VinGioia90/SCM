#include <Rcpp.h>
// [[Rcpp::export(name="lt_inversion")]]
Rcpp::NumericMatrix lt_inversion(Rcpp::NumericMatrix& X){
 using namespace Rcpp;

 uint32_t d = X.cols();
 uint32_t j;
 uint32_t k;
 uint32_t l;
 double s;
 NumericMatrix res(d,d);

 for(j = 0; j < d; j++){
  res(j,j) = 1.0;
  if(j > 0){
   for(k = 0; k < j; k++){
    s = 0.0;
    for(l = k; l < j; l++){
     s += X(j,l) * res(l,k);
    }
    res(j,k) = -s;
   }
  } 
 }

 return(res);
}
