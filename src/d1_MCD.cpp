#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(name = "score_mcd")]]
double d1_mcd(NumericMatrix& eta, NumericMatrix& y, NumericMatrix& res) {
 uint32_t n = y.rows();
 uint32_t d = y.cols();
 
 IntegerVector z(d * (d - 1)/2);
 IntegerVector w(d * (d - 1)/2);
 IntegerMatrix G(d - 1,d - 1);

 uint32_t count = 0;
 uint32_t count2 = 0;

 uint32_t i;
 uint32_t j;
 uint32_t k;

 for(i = 0; i < d; i++){
  if(i < d - 1){
   for(j = 0; j < (i + 1); j++){
    G(i,j) = count2 + 2 * d;
    count2 = count2 + 1;
   }
  }
  for(j = 0; j < i; j++){
   z(count) = j;
   w(count) = i;
   count = count + 1;
  }
 }
 double aux_out = 0.0;
 NumericVector s(d);

 for(i = 0; i < n; i++){
 
  s(0) = y(i,0) - eta(i,0);
  for(j = 1; j < d; j++){
   for(k = 0; k < j; k++){
    aux_out = aux_out + (y(i,k) - eta(i,k)) * eta(i,G(j - 1,k)); 
    }
   s(j) = aux_out + y(i,j) - eta(i,j); 
   aux_out = 0.0;
  }

  res(i,d-1) = exp(-eta(i,2*d-1)) * s(d-1);
  res(i,2*d-1) = -0.5 + 0.5*exp(-eta(i,2*d-1)) * s(d-1) * s(d-1);  
  if(d > 2){
   res(i,3 * d - 1) = -exp(-eta(i,w(d-1) + d)) * s(w(d-1)) * (y(i,z(d-1)) - eta(i,z(d-1)));  
  }

   for(j = 0; j < d - 1; j++){    
    res(i,j) = exp(-eta(i,j + d)) * s(j);
     for(k = j + 1; k < d; k++){
      res(i,j) = res(i,j) + exp(-eta(i,k + d)) * s(k) * eta(i,G(k - 1,j));
     }
    res(i,j + d) = -0.5 + 0.5*exp(-eta(i,j + d)) * s(j) * s(j);  
    res(i,j + 2 * d) = -exp(-eta(i,w(j) + d)) * s(w(j)) * (y(i,z(j)) - eta(i,z(j)));  
   } 

   for(j = d; j < d * (d - 1)/2; j++){
    res(i,j + 2 * d) = -exp(-eta(i,w(j) + d)) * s(w(j)) * (y(i,z(j)) - eta(i,z(j)));  
   }     
 }

 return(1.0);
}