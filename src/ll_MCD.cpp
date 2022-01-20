#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(name = "ll_mcd")]]
double logL_mcd(NumericMatrix& eta, NumericMatrix& y){
 uint32_t n = y.rows();
 uint32_t d = y.cols();

 double out = 0.0;
 double aux_out = 0.0; 

 uint32_t i;
 uint32_t j;
 uint32_t k;
 uint32_t count;

 for(i = 0; i < n; i++){  
  out = out + 0.5*eta(i,d) + 0.5*exp(-eta(i,d))*(y(i,0)-eta(i,0))*(y(i,0)-eta(i,0));  
  
  count = 2*d;
  for(j = d+1; j < 2*d; j++){
   for(k = d; k < j; k++){
    aux_out = aux_out + (y(i,k-d)-eta(i,k-d))*eta(i,count);
    count = count+1; 
   }
   out = out + 0.5*eta(i,j) + 0.5*exp(-eta(i,j))*(aux_out+y(i,j-d)-eta(i,j-d))*(aux_out+y(i,j-d)-eta(i,j-d));
   aux_out = 0.0;
  }  
 }
 return -out;
}


