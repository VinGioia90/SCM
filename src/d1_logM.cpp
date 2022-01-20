// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Score
//'
//' @param eta Linear predictor (n x (d + dx(d+1)/2) matrix).
//' @param y Outcome (n x d matrix).
//' @param res Memory initialization
//' @export
// [[Rcpp::export(name = "score_logm")]]
double  d1_logm(arma::mat& eta, arma::mat& y, arma::mat& res) {
 using namespace arma;

 uint32_t n = y.n_rows;
 uint32_t d = y.n_cols;
 uint32_t i;
 uint32_t j;
 uint32_t k;
 uint32_t l;
 uint32_t count;

 vec auxout(d,fill::zeros);

 rowvec ri(d, fill::zeros);
 mat Ai(d,d,fill::zeros);
 
 vec eigval(d,fill::zeros);
 mat eigvec(d,d,fill::zeros);
 mat Vbar(d,d,fill::zeros);
 mat Di(d,d,fill::zeros);
 mat loewner(d,d,fill::zeros);
 

 cube V_j (d,d,d*(d+1)/2,fill::zeros);
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
   
  for(j = 0; j < d; j++){
    ri(j) = y(i,j) - eta(i,j);
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
  
  Di = diagmat(exp(-eigval));
  auxout = eigvec*Di*eigvec.t()*ri.t();
  res.submat(i,0,i,d-1) =  auxout.t();
   
   for( k = 0; k < d; k++){
     loewner(k,k) = exp(-eigval(k));
     if(k < d-1){
     for(l=k+1; l < d; l++ ){
       loewner(k,l)=((exp(-eigval(k))-exp(-eigval(l)))/
                                    (eigval(l)-eigval(k)));
       loewner(l,k) = loewner(k,l);
      }
     }
    }

  for(j = d; j<(d+d*(d+1)/2) ; j++){ 
    Vbar = eigvec.t()*(-V_j.slice(j-d))*eigvec;
    res(i,j) = -0.5*as_scalar(ri*eigvec*(Vbar % loewner)*eigvec.t()*ri.t());
    if(j < 2*d ){
     res(i,j) = res(i,j)-0.5;
    }
   }
 }
 return(1.0);
}