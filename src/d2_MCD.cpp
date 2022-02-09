#include <Rcpp.h>
using namespace Rcpp;

//' Hessian mcd
//'
//' @param eta Linear predictor (n x (d + dx(d+1)/2) matrix).
//' @param y Outcome (n x d matrix).
//' @param res Memory initialization
//' @export

// [[Rcpp::export(name="d2_mcd")]]
double d2_mcd(NumericMatrix& eta, NumericMatrix& y, NumericMatrix& res){
 uint32_t n = y.rows();
 uint32_t d = y.cols();
 uint32_t l = eta.cols();

 IntegerVector z(d*(d - 1)/2);
 IntegerVector w(d*(d - 1)/2);
 IntegerMatrix G(d - 1,d - 1);

 double aux_out = 0.0;
 double aux_out2 = 0.0;
 NumericMatrix out(l,l);
 NumericVector s(d);

 uint32_t count = 0;
 uint32_t count2 = 0;

 uint32_t i;
 uint32_t j;
 uint32_t k;
 uint32_t q;

 for(i = 0; i < d; i++){
  if(i < d - 1){
   for(j = 0; j < i + 1; j++){
    G(i,j) = count2 + 2*d;
    count2 = count2 + 1;
   }
  }
  for(j = 0; j < i; j++){
   z(count) = j;
   w(count) = i;
   count = count + 1;
  }
 }

 
 IntegerVector t = w-z; 
 


 for(i = 0; i < n; i++){

  s(0) = y(i,0) - eta(i,0);
  for(j = 1; j < d; j++){
   for(k = 0; k < j; k++){
    aux_out = aux_out + (y(i,k) - eta(i,k)) * eta(i,G(j - 1,k)); 
    }
   s(j) = aux_out + y(i,j) - eta(i,j); 
   aux_out = 0.0;
  }

 for(j = 0; j < d; j++){
  // Blocco (1,1)
  out(j,j) = -exp(-eta(i,j + d));
  for(k = j + 1; k < d; k++){
  // Blocco (1,1)
   out(j,k) = -exp(-eta(i,k + d))*eta(i,G(k - 1,j));
   for(q = k + 1; q < d; q++){
    aux_out2 = aux_out2 - exp(-eta(i,q + d))*eta(i,G(q - 1,j))*eta(i,G(q - 1,k));
   }
   // Blocco (1,1)
   out(j,k) = out(j,k) + aux_out2;
   aux_out2 = 0.0;
   // Blocco (1,2)
   out(j,k + d) = -exp(-eta(i,k + d))*s(k)*eta(i,G(k - 1,j));
   aux_out = aux_out - exp(-eta(i,k + d))*eta(i,G(k - 1,j))*eta(i,G(k - 1,j));
  }
  //Blocco (1,1)
  out(j,j) = out(j,j) + aux_out;
  aux_out = 0.0;
  
  // Blocco (1,2)
  out(j,j + d) = -exp(-eta(i,j + d))*s(j);
  // Blocco (2,2) 
  out(j + d,j + d) = -0.5*exp(-eta(i,j + d))*s(j)*s(j);
 }
   
  for(k = 0; k < d*(d - 1)/2; k++){
   //Blocco (1,3)
   for(int j = 0; j < w(k); j++){
    out(j,k + 2*d) = exp(-eta(i,w(k) + d))*(y(i,z(k)) - eta(i,z(k)))*eta(i,G(w(k) - 1,j));
   } 
   out(z(k),k + 2*d) = out(z(k),k + 2*d) + exp(-eta(i,w(k) + d))*s(w(k));
   out(w(k),k + 2*d) = exp(-eta(i,w(k) + d))*(y(i,z(k)) - eta(i,z(k)));
   //Blocco (2,3)
   out(w(k) + d,k + 2*d) = exp(-eta(i,w(k) + d))*s(w(k))*(y(i,z(k)) - eta(i,z(k)));
  //Blocco (3,3)  
  for(int j = 0; j < t(k); j++){
     out(k+2*d,k+2*d+j) = -exp(-eta(i,w(k) + d))*(y(i,z(k+j)) - eta(i,z(k+j)))*(y(i,z(k)) - eta(i,z(k)));
    }
  }   

  count = 0;   
  for(j = 0; j < l; j++){
   for(k = j; k < l; k++){
    res(i, count) = out(j,k);
    count = count + 1; 
   }
  }
 }
 return(1.0);
}