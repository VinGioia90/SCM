#ifndef _SCM_H
#define _SCM_H

#include <RcppArmadillo.h>
using namespace Rcpp;
NumericMatrix lt_inversion(NumericMatrix& X);

NumericMatrix mcd_LD(NumericVector& x, uint32_t& d);

NumericMatrix mcd_Sigma(NumericMatrix& L, NumericMatrix& D,  uint32_t& d);



#endif