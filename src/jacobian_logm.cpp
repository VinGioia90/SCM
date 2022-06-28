#include "scm.h"

// [[Rcpp::export(name="jacobian_logm")]]
double jacobian_logm(Rcpp::NumericMatrix& eta,  Rcpp::NumericMatrix& res, uint32_t& d,
                     uint32_t& S_row, uint32_t& S_col, uint32_t& cor_flag){
  using namespace Rcpp;
  uint32_t n = eta.rows();
  uint32_t i;
  uint32_t j;
  uint32_t k;
  uint32_t l;
  uint32_t count;

  arma::mat auxout(d,d, arma::fill::zeros);
  arma::rowvec ri(d, arma::fill::zeros);
  arma::mat Ai(d,d,arma::fill::zeros);
  arma::vec eigval(d,arma::fill::zeros);
  arma::mat eigvec(d,d,arma::fill::zeros);
  arma::mat Wbar(d,d,arma::fill::zeros);
  arma::mat Di(d,d,arma::fill::zeros);
  arma::mat loewner(d,d,arma::fill::zeros);
  arma::cube V_j (d,d,d*(d+1)/2,arma::fill::zeros);

  NumericVector etai(d);
  NumericMatrix Sigma(d,d);


  for(j = 0; j < d; j++){
    V_j(j,j,j) = 1;
  }
  j = d;
  for(k = 1; k < d; k++){
    for(l = 0; l < k; l++){
      if(k != l){
        V_j(k,l,j) = 1;
        V_j(l,k,j) = 1;
        j = j + 1;
      }
    }
  }


  for(i = 0; i < n; i++){
    etai = eta.row(i);
    Sigma = logM_Sigma(etai, d);

    for(j = 0; j < d; j++){
      Ai(j,j) = eta(i,j + d);
    }
    count = 0;
    for(j=1; j<d; j++){
      for(k = 0; k < j; k++){
        Ai(j,k)= eta(i,count+2*d);
        Ai(k,j) = eta(i,count+2*d);
        count += 1;
      }
    }
    eig_sym(eigval, eigvec, Ai);
    Di = diagmat(exp(eigval));

    for( k = 0; k < d; k++){
      loewner(k,k) = exp(eigval(k));
      if(k < d-1){
        for(l=k+1; l < d; l++ ){
          loewner(k,l)=((exp(eigval(k))-exp(eigval(l)))/
                          (eigval(k)-eigval(l)));
          loewner(l,k) = loewner(k,l);
        }
      }
    }

    for(j = 0; j < d; j++){
      if((cor_flag == 0) || (S_row == S_col)){
        Wbar = eigvec.t()*V_j.slice(j)*eigvec;
        auxout =   eigvec*(Wbar % loewner)*eigvec.t();
        res(i, j + d) = auxout(S_row, S_col);
      }
      if((S_row != S_col) && (cor_flag == 1)){
        Wbar = eigvec.t()*V_j.slice(j)*eigvec;
        auxout =   eigvec*(Wbar % loewner)*eigvec.t();
        res(i, j + d) = pow(Sigma(S_row, S_row), -0.5) * auxout(S_row,S_col) * pow(Sigma(S_col, S_col), -0.5) -
          0.5 * Sigma(S_col, S_row) * pow(Sigma(S_row,S_row), -1.5) *  pow(Sigma(S_col, S_col),-0.5) * auxout(S_row,S_row) -
          0.5 * Sigma(S_col, S_row) * pow(Sigma(S_row, S_row), -0.5) *  pow(Sigma(S_col, S_col),-1.5) * auxout(S_col,S_col) ;

      }
    }

    count = 0;
    for(j = d; j < (d * (d + 1)/2); j++){
      if((cor_flag == 0) || (S_row == S_col)){
        Wbar = eigvec.t()*V_j.slice(j)*eigvec;
        auxout =   eigvec*(Wbar % loewner)*eigvec.t();
        res(i, j + d) = auxout(S_row, S_col);
      }
      if((S_row != S_col) && (cor_flag == 1)){
        Wbar = eigvec.t()*V_j.slice(j)*eigvec;
        auxout =   eigvec*(Wbar % loewner)*eigvec.t();

        res(i, j + d) = pow(Sigma(S_row, S_row), -0.5) *  auxout(S_row,S_col) * pow(Sigma(S_col, S_col), -0.5) -
          0.5 * Sigma(S_col, S_row) * pow(Sigma(S_row,S_row), -1.5) *  pow(Sigma(S_col, S_col),-0.5) * auxout(S_row,S_row) -
          0.5 * Sigma(S_col, S_row) * pow(Sigma(S_row, S_row), -0.5) *  pow(Sigma(S_col, S_col),-1.5) * auxout(S_col,S_col)  ;
      }
      count += 1;
    }
  }
  return (1.0);
}
