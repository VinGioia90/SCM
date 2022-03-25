#include "scm.h"

// [[Rcpp::export(name="pred_logm")]]
double pred_logm(Rcpp::NumericMatrix& eta, Rcpp::NumericMatrix& res,  uint32_t& d){
using namespace Rcpp;
 uint32_t n = eta.rows();
 uint32_t i;
 uint32_t j;
 uint32_t k;
 uint32_t count;

 NumericMatrix Sigma(d,d);
 NumericVector etai(d);

 for(i = 0; i < n; i++){
  etai = eta(i,_);

  Sigma = logM_Sigma(etai, d);
     
  count = d; 
  for(j = 0; j < d; j++){
   res(i, j) = eta(i, j);
   res(i, j + d) = Sigma(j, j); 
    if(j > 0){
     for(k = 0; k < j; k++){
      res(i, count + d) = Sigma(k, j);
      count = count + 1;
     }
    }  
   }
  }

 return (1.0);
}