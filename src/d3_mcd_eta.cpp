// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


// [[Rcpp::export(name="d3_mcd_eta")]]
double d3_mcd_eta(const arma::mat& eta, const arma::mat& y, arma::mat& res,
                  arma::vec& z, arma::vec& w, arma::mat& G, arma::vec& t,
                  Rcpp::IntegerVector& h, Rcpp::IntegerVector& h2, Rcpp::IntegerVector& h3,
                  Rcpp::List& idx3_1, Rcpp::List& idx3_2,  Rcpp::List& idx3_3,
                  Rcpp::List& idx3_4, Rcpp::List& idx3_5, Rcpp::List& idx3_6, Rcpp::List& idx_neq0){
  using namespace arma;
  uint32_t n = y.n_rows;
  uint32_t d = y.n_cols;
  uint32_t q = eta.n_cols;
  uint32_t q2 = q*(q+1)/2;

  double aux_out = 0.0;

  mat out(q2,q2, fill::zeros);
  rowvec r(d, fill::zeros);
  vec s(d, fill::zeros);

  Rcpp::IntegerVector aux_idx3_1;
  Rcpp::IntegerVector aux_idx3_2;
  Rcpp::IntegerVector aux_idx3_3;
  Rcpp::IntegerVector aux_idx3_4;
  Rcpp::IntegerVector aux_idx3_5;
  Rcpp::IntegerVector aux_idx3_6;
  Rcpp::IntegerMatrix idx_k;

  uint32_t count = 0;

  uint32_t i;
  uint32_t j;
  uint32_t k;
  uint32_t l;
  uint32_t m;

  for(i = 0; i < n; i++){
    r = y.row(i) - eta(i,span(0, d-1));
    s(0) = r(0);

    for(j = 1; j < d; j++){
      for(k = 0; k < j; k++){
        aux_out = aux_out + r(k) * eta(i,G(j - 1,k));
      }
      s(j) = aux_out + r(j);
      aux_out = 0.0;
    }

    count = 0;
    for(k = 0; k < 2*d; k++){

      if(k < d){
        // Block (1,1,1) is zero everywhere

        // Block (1,1,2)
        out(k,k+d) = exp(-eta(i,k+d));
        for(j = k+1; j <d; j++){
          out(j,j+d)= exp(-eta(i,j+d))*eta(i,G(j-1,k));
        }

        if(k<(d-1)){
          // Block (1,1,2)
          for(j = 0; j < h(d-k-2); j++){
            out(z(j)+k,w(j)+k+d)=exp(-eta(i,w(j)+k+d))*eta(i,G(w(j)+k-1,k))*eta(i,G(w(j)+k-1,z(j)+k));
          }
          // Block (1,1,3)
          for(j = 0; j < d-k-1; j++){
            out(j+k+1, G(j+k,k))=-exp(-eta(i,j+k+d+1));
            out(k,G(j+k,k)) = -2*exp(-eta(i,j+k+d+1))*eta(i, G(j+k,k));
          }
          // Block (1,2,3)
          for(j = 0; j < (d*(d-1)/2+k-(h(k)-1)); j++){
            out(w(j-k+h(k)-1)+d,G(w(j-k+h(k)-1)-1,z(j+h(k)-k-1)))= -exp(-eta(i,w(j-k+h(k)-1)+d))*r(z(j+h(k)-k-1))*eta(i,G(w(j+h(k)-k-1)-1,k));
          }
          for(j = 0; j < d-k-1; j++){
            out(j+k+d+1, G(j+k, k)) -=  exp(-eta(i,j+k+d+1))*s(j+k+1);
          }
          // Block (1,3,3)
          aux_idx3_3 = idx3_3[k];
          aux_idx3_4 = idx3_4[k];
          aux_idx3_5 = idx3_5[k];
          aux_idx3_6 = idx3_6[k];
          for(j = 0; j <d-k-1; j++){
            out(G(j+k,k), G(j+k,k))=2*exp(-eta(i,k+d+j+1))*r(k);
          }
          for(j = 0; j < h3(k); j++){
            out(G(aux_idx3_3(j),aux_idx3_4(j)), G(aux_idx3_3(j),aux_idx3_5(j)))=exp(-eta(i,w(G(aux_idx3_3(j),aux_idx3_4(j))-2*d)+d))*r(aux_idx3_6(j));
          }
        }

        if(k<d-2){
          // Block (1,1,3)
          aux_idx3_1 = idx3_1[k];
          aux_idx3_2 = idx3_2[k];
          for(j = 0; j < h2(k); j++){
            out(aux_idx3_1(j), G(aux_idx3_2(j),k))=-exp(-eta(i,aux_idx3_2(j)+d+1))*eta(i,G(aux_idx3_2(j),aux_idx3_1(j)));
            out(aux_idx3_1(j), G(aux_idx3_2(j),aux_idx3_1(j)))=-exp(-eta(i,aux_idx3_2(j)+d+1))*eta(i,G(aux_idx3_2(j),k));
          }
        }

        //Block (1,2,2)
        out(k+d,k+d)=  exp(-eta(i,k+d))*s(k);
        for(j = k+d+1; j <2*d; j++){
          out(j,j) = exp(-eta(i,j))*s(j-d)*eta(i,G(j-d-1,k));
        }

        //Block (1,2,3)


        if(k > 0){
          for(j = 0; j < k; j++){
            out(k+d, G(k-1, j)) =  -exp(-eta(i,k+d))*r(j);
          }
        }

      }
      // Block (2,2,2)  - dealing only one element
      if(k==d){
        out(k,k)=0.5*exp(-eta(i,k))*r(k-d)*r(k-d);
      }
      if((k > d) && (k < 2*d)){
        // Block (2,2,2) - dealing the remaining elements
        out(k,k)=0.5*exp(-eta(i,k))*s(k-d)*s(k-d); //ok

        // Block (2,2,3)
        for(j = 0; j < k-d; j++){
          out(k,G(k-d-1,j)) = -exp(-eta(i,k))*s(k-d)*r(j);
        }

        // Block (2,3,3)
        for(j = 0; j < h(k-d-1); j++){
          out(G(k-d-1,z(j)),G(k-d-1,w(j)-1))=exp(-eta(i,k))*r(z(j))*r(w(j)-1);
        }
      }
      // Block (3,3,3) is zero everywhere

  idx_k = Rcpp::as<Rcpp::IntegerMatrix>(idx_neq0[k]);
  for(j = 0; j < idx_k.cols(); j++){
    res(i, count) = out( idx_k(0,j),idx_k(1,j));
    out( idx_k(0,j),idx_k(1,j)) = 0;
    count = count + 1;
  }

}
}

return(1.0);
}
