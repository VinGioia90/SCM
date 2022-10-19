#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export(name="idx_zwGt")]]
double idx_zwGt(uint32_t d, IntegerVector& z, IntegerVector& w, IntegerMatrix& G, IntegerVector& t){
  uint32_t count = 0;
  uint32_t count2 = 0;
  uint32_t i;
  uint32_t j;

  for(i = 0; i < d; i++){
    if ( i < d - 1 ) {
      for(j = 0; j < i + 1; j++){
        G(i, j) = count2 + 2 * d;
        count2 += 1;
      }
    }
    for(j = 0; j < i; j++){
      z(count) = j;
      w(count) = i;
      count += 1;
    }
  }
  t = w - z;
  return(1.0);
}
