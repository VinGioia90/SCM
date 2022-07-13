#include "scm_der.h"

// [[Rcpp::export(name="d2_mcd_beta")]]
double d2_mcd_beta(arma::mat& X,  arma::mat& eta,  arma::mat& y,  Rcpp::List& jj,  uint32_t& K,  Rcpp::List& idx_jk, arma::mat& lbb, arma::mat& l2, arma::mat& l2_last,  arma::uvec& idx_g,  arma::vec& z, arma::vec& w, arma::mat& G, arma::vec& t){
  using namespace arma;
  uint32_t size_g = idx_g.n_elem-1;

  uint32_t i;
  uint32_t j;
  uint32_t k;
  Rcpp::IntegerVector idx_jjj;
  Rcpp::IntegerVector idx_jjk;
  Rcpp::IntegerVector idx_k;
  int ljjj;
  int ljjk;


  uint32_t count;

  for(i = 0; i < size_g; i++){
    if(i < (size_g-1)){
      d2_mcd_eta(eta.rows(idx_g(i)+1,idx_g(i+1)), y.rows(idx_g(i)+1,idx_g(i+1)), l2, idx_jk,z,w,G,t);
    } else {
       d2_mcd_eta(eta.rows(idx_g(i)+1,idx_g(i+1)), y.rows(idx_g(i)+1,idx_g(i+1)), l2_last, idx_jk,z,w,G,t);
    }
    count = 0;
     for(j = 0; j < K; j++){
       idx_jjj = jj[j];
       idx_k = idx_jk[j];
       ljjj=  idx_jjj.length()-1;
       Col<int> idx_jjja(idx_jjj.begin(),idx_jjj.length(),false);


       for(k = 0; k < idx_k.length(); k++){
           idx_jjk = jj[idx_k(k)];
           ljjk =  idx_jjk.length()-1;
           Col<int> idx_jjka(idx_jjk.begin(),idx_jjk.length(),false);
           if(i < (size_g-1)){
             lbb.submat(idx_jjja(0),idx_jjka(0),idx_jjja(ljjj),idx_jjka(ljjk)) += X.submat(idx_g(i)+1,idx_jjja(0),idx_g(i+1),idx_jjja(ljjj)).t()*(X.submat(idx_g(i)+1,idx_jjka(0),idx_g(i+1),idx_jjka(ljjk)).each_col() % l2.col(count));
           } else {
             lbb.submat(idx_jjja(0),idx_jjka(0),idx_jjja(ljjj),idx_jjka(ljjk)) += X.submat(idx_g(i)+1,idx_jjja(0),idx_g(i+1),idx_jjja(ljjj)).t()*(X.submat(idx_g(i)+1,idx_jjka(0),idx_g(i+1),idx_jjka(ljjk)).each_col() % l2_last.col(count));
           }

           if((idx_k(k) > j) & (i == (size_g-1))){
             lbb.submat(idx_jjka(0),idx_jjja(0),idx_jjka(ljjk),idx_jjja(ljjj)) = lbb.submat(idx_jjja(0),idx_jjka(0),idx_jjja(ljjj),idx_jjka(ljjk)).t();
           }
           count += 1;
         }
       }

    }

    return(1.0);
  }
