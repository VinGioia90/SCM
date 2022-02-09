// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Hessian logm
//'
//' @param eta Linear predictor (n x (d + dx(d+1)/2) matrix).
//' @param y Outcome (n x d matrix).
//' @param res Memory initialization

double d2_logm(arma::mat& eta, arma::mat& y, arma::mat& res){
using namespace arma;
 uint32_t n = y.n_rows;
 uint32_t d = y.n_cols;
 uint32_t s = d+d*(d+1)/2;

 uint32_t i;
 uint32_t j;
 uint32_t k;
 uint32_t l;
 uint32_t r;
 uint32_t q;
 uint32_t t;
 uint32_t count;

 mat auxout(s,s,fill::zeros);

 rowvec ri(d, fill::zeros);

  mat Ai(d,d,fill::zeros);
  mat Omega(d,d,fill::zeros);
 vec eigval(d, fill::zeros);
 mat eigvec(d, d, fill::zeros);
 
 mat Vbar1(d,d, fill::zeros);
 mat Vbar2(d,d, fill::zeros);
 mat Di(d,d, fill::zeros);
 mat loewner(d,d, fill::zeros);

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
     Ai(j,k) = eta(i,count+2*d);
     Ai(k,j) = eta(i,count+2*d);
     count += 1;
    }
   } 
   eig_sym(eigval, eigvec, Ai);
   Di= diagmat(exp(-eigval));
   auxout.submat(0,0,d-1,d-1) =  -eigvec*Di*eigvec.t();
 
 for(k=0; k < d; k++){
  loewner(k,k) = exp(-eigval(k));
   if(k<d-1){
    for(l=k+1; l<d; l++){
      loewner(k,l)=((exp(-eigval(k))-exp(-eigval(l)))/
                      (eigval(l)-eigval(k)));
      loewner(l,k)=loewner(k,l);
    }
   }   
  }
  
  for(j = 0; j < d ; j++){ 
    for(k=d; k<d+d*(d+1)/2; k++){
     Vbar1 = eigvec.t()*(-V_j.slice(k-d))*eigvec;
      auxout.submat(0,k,d-1,k)= eigvec*(Vbar1 % loewner)*eigvec.t()*ri.t();
     }   
   } 

 for(j = d; j < d+d*(d+1)/2 ; j++){ 
     Vbar1 = eigvec.t()*(-V_j.slice(j-d))*eigvec;
  for(k=d; k<d+d*(d+1)/2; k++){
     Vbar2 = eigvec.t()*(-V_j.slice(k-d))*eigvec;
     for(r=0; r<d; r++){
      Omega(r,r) = Vbar1(r,r)*Vbar2(r,r)*exp(-eigval(r));
      for(q=0; q<d; q++){
       if(q != r){
        Omega(r,r)=Omega(r,r)+(Vbar1(r,q)*Vbar2(q,r)+Vbar1(q,r)*Vbar2(r,q))*(exp(-eigval(r))/(eigval(q)-eigval(r))+(exp(-eigval(q))-exp(-eigval(r)))/((eigval(q)-eigval(r))*(eigval(q)-eigval(r))));
        Omega(r,q)=(Vbar1(r,r)*Vbar2(r,q)+Vbar1(r,q)*Vbar2(r,r))*(exp(-eigval(r))/(eigval(q)-eigval(r))+(exp(-eigval(q))-exp(-eigval(r)))/((eigval(q)-eigval(r))*(eigval(q)-eigval(r))))+
                   (Vbar1(r,q)*Vbar2(q,q)+Vbar1(q,q)*Vbar2(r,q))*(exp(-eigval(q))/(eigval(r)-eigval(q))+(exp(-eigval(r))-exp(-eigval(q)))/((eigval(q)-eigval(r))*(eigval(q)-eigval(r))));
       for(t=0; t<d; t++){
        if(t != q){
         if(t != r){
          Omega(r,q)=Omega(r,q)+(Vbar1(r,t)*Vbar2(t,q)+Vbar1(t,q)*Vbar2(r,t))*((exp(-eigval(q))-exp(-eigval(r)))/((eigval(r)-eigval(q))*(eigval(t)-eigval(q)))+(exp(-eigval(t))-exp(-eigval(r)))/((eigval(r)-eigval(t))*(eigval(q)-eigval(t)))); 
         }
        }     
       }
       }
      }
     }
     auxout(j,k)=-0.5*as_scalar(ri*eigvec*Omega*eigvec.t()*ri.t());
     Omega.fill(0);
    }    
   }

   count = 0;
   for(j = 0; j<s; j++){
    for(k = j; k<s; k++){ 
     res(i,count) = auxout(j,k); 
     count = count + 1;
    }
   }
 }

 return(1.0);

}



