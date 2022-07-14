#include "scm_der2.h"

//' Score mcd in beta
//'
//' @param X model matrix
//' @param jj List of lpi index
//' @param K number of lpi
//' @param lb matrix of the derivatives beta (initialization)
//' @param l1 matrix of the derivatives eta (initialization)
//' @param eta lpi
//' @param y outcome matrix
//' @param z auxiliary index vector z
//' @param w auxiliary index vector w
//' @param G auxiliary index matrix G
//'
//' @export
// [[Rcpp::export(name="d1_mcd_beta")]]
double d1_mcd_beta(arma::mat& X, Rcpp::List& jj, uint32_t& K, arma::vec& lb, arma::mat& l1, arma::mat& eta,  arma::mat& y, arma::vec& z, arma::vec& w, arma::mat& G){
  using namespace arma;
  uint32_t j;
  Rcpp::IntegerVector idx_jjj;

  // eventualmente fare su gruppi  for(i = 0; i < size_g; i++){
    d1_mcd_eta(eta,y, l1, z, w, G);

    for(j = 0; j < K; j++) {
      idx_jjj = jj[j];
      uint32_t l_jjj = idx_jjj.length()-1;
      Col<int> idx_jjja(idx_jjj.begin(),idx_jjj.length(),false);
      lb.subvec(idx_jjja(0),idx_jjja(l_jjj)) = X.cols(idx_jjja(0),idx_jjja(l_jjj)).t() * l1.col(j);


  }

  return(1.0);
}
