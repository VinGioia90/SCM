// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(name="d1_logm_eta")]]
double  d1_logm_eta(const arma::mat& eta, const arma::mat& y, arma::mat& d1l, arma::vec& z, arma::vec& w) {
  using namespace arma;
  
  uint32_t n = y.n_rows;
  uint32_t d = y.n_cols;
  uint32_t q = eta.n_cols;
  
  uint32_t i;
  uint32_t j;
  uint32_t k;
  uint32_t c1;
  
  rowvec r(d, fill::zeros); // residual vector
  mat Theta(d, d, fill::zeros); //log Sigma matrix
  vec L(d, fill::zeros); //eigenvalues - Lambda in the paper
  mat U(d, d, fill::zeros); //eigenvector
  
  // useful vector and matrices
  vec s(d, fill::zeros);
  mat loewner(d, d, fill::zeros);
  mat F(d, d, fill::zeros);
  mat Pi(d, d, fill::zeros);
  mat Xi(d, d, fill::zeros);
  
  for(i = 0; i < n; i++){
    // i-th residuals
    r = y.row(i) - eta(i, span(0, d - 1));
    
    // log Sigma matrix
    c1 = 0;
    for(j = 1; j < d; j++){
      for(k = 0; k < j; k++){
        Theta(j, k) = eta(i, c1 + 2 * d);
        Theta(k, j) = Theta(j, k);
        c1 += 1;
      }
    }
    Theta.diag() = eta(i, span(d, 2 * d - 1));
    
    // Eigen decomposition
    eig_sym(L, U, Theta);
    
    // Loewner Matrix building
    for(j = 0; j < (d - 1); j++){
      for(k = j + 1; k < d; k++){
        loewner(j, k) = ((exp(-L(j)) - exp(-L(k)))/(L(k) - L(j)));
        loewner(k, j) = loewner(j, k);
      }
    }
    loewner.diag() =  exp(-L);
    
    // Building the needed quantities
    s = U.t() * r.t();
    F = U.each_row() % s.t();
    Pi = F * loewner;
    Xi = Pi * F.t();
    
    // Score elements 1, ..., d
    d1l(i, span(0, d - 1)) = (U * (loewner.diag() % s)).t();
    
    // Score elements d+1, ..., 2d
    d1l(i, span(d, 2 * d - 1)) = (-0.5 + 0.5 * Xi.diag()).t();
    
    // Score elements 2d+1, ..., q
    for(j = 2 * d; j < q; j++){
      d1l(i, j) = Xi(z(j - 2 * d), w(j - 2 * d));
    }
  }
  return(1.0);
}