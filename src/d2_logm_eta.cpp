// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export(name="d2_logm_eta")]]
double  d2_logm_eta(const arma::mat& eta, const arma::mat& y, arma::mat& d2l, arma::vec& d2l_v, arma::vec& z, arma::vec& w, Rcpp::List&  b1_eta, Rcpp::List&  b3, arma::vec& ib1_eta, arma::vec& ib3) {
  using namespace arma;

  uint32_t n = y.n_rows;
  uint32_t d = y.n_cols;
  uint32_t q = eta.n_cols;

  uint32_t i;
  uint32_t j;
  uint32_t k;
  uint32_t kp;
  uint32_t l;
  uint32_t m;
  uint32_t c1;
  uint32_t zl;// auxiliary access to z elements
  uint32_t wl;// auxiliary access to w elements
  uint32_t zm;// auxiliary access to z elements
  uint32_t wm;// auxiliary access to w elements

  Rcpp::IntegerVector ik1;
  vec::iterator it_b1;
  vec::iterator it_b1_end;

  vec::iterator it_b3;
  vec::iterator it_b3_end;

  rowvec r(d, fill::zeros); //i-th residual vector
  mat Theta(d, d, fill::zeros); //i-th log Sigma matrix
  vec L(d, fill::zeros); //eigenvalues - Lambda in the paper
  mat U(d, d, fill::zeros); //eigenvectors

  vec s(d, fill::zeros);
  mat out(q, q, fill::zeros);
  mat loewner(d, d, fill::zeros);
  mat F(d, d, fill::zeros);
  mat Pi(d, d, fill::zeros);
  mat tildeU(d, d * (d - 1)/2, fill::zeros); //U tilde
  mat tilde_loew(d, d, fill::zeros); // tilde loewner
  mat A(d, d, fill::zeros);
  cube star_loew(d, d, d, fill::zeros); // star loewner
  cube A1(d, d, d, fill::zeros);
  cube A2(d, d, d, fill::zeros);

  for(i = 0; i < n; i++){
    // i-th residuals vector
    r = y.row(i) - eta(i, span(0, d - 1));

    // i-th logarithm of the covariance matrix
    c1 = 0;
    for(j = 1; j < d; j++){
      for(k = 0; k < j; k++){
        Theta(j, k) = eta(i, c1 + 2 * d);
        Theta(k, j) = Theta(j, k);
        c1 += 1;
      }
    }
    Theta.diag() = eta(i, span(d, 2 * d - 1));

    // eigen decomposition
    eig_sym(L, U, Theta);

    //loewner matrix
    for(j = 0; j < (d - 1); j++){
      for(k = j + 1; k < d; k++){
        loewner(j, k) = ((exp(-L(j)) - exp(-L(k)))/(L(k) - L(j)));
        loewner(k, j) = loewner(j, k);
      }
    }
    loewner.diag() =  exp(-L);

    // Building the auxiliaries star loewner and tilde loewner
    for(j = 0; j < d; j++){
      for(k = 0; k < d; k++){
        if(k != j) {
          tilde_loew(j, k) = (loewner(j, k) - loewner(j, j))/(L(j) - L(k));
          for(kp = 0; kp < d; kp++){
            if(kp != k and kp!= j) star_loew(j, k, kp) = (loewner(kp, k) - loewner(j, k))/(L(j) - L(kp));
          }
        }
      }
    }
    // Other quantities
    s = U.t() * r.t();
    F = U.each_row() % s.t();
    Pi = F * loewner;
    A = F * tilde_loew.t();
    for(j = 0; j < d; j++){
      A1.slice(j) = F * (F.each_row() % tilde_loew.col(j).t()).t();
      A2.slice(j) = F  * (F * star_loew.slice(j)).t(); //star_loew.slice(j).t() is unuseful since star_loew is symmetric
    }
    for(l = 0; l < d * (d - 1)/2; l++){
      zl = z(l);
      wl = w(l);
      for(j = 0; j < d; j++){
        if(j == zl) tildeU(zl, l) = Pi(zl, zl) * U(wl, zl) + Pi(wl, zl) * U(zl, zl);
        if(j == wl) tildeU(wl, l) = Pi(wl, wl) * U(zl, wl) + Pi(zl, wl) * U(wl, wl);
        if(j != zl and j != wl) tildeU(j, l) = Pi(zl, j) * U(wl, j) + Pi(wl, j) * U(zl, j);
      }
    }


    for(l = 0; l < d; l++){
      //Block (1,1) + Block (2,2)
      for(m = l; m < d; m++){
        out(l, m) = -as_scalar(U.row(l) * ((U.row(m)).t() % loewner.diag()));
        for(j = 0; j < d; j++){
          out(l + d, m + d) -= U(l, j) * U(m, j) * (F(l, j) * F(m, j) * loewner(j, j)/2 + F(m, j) * A(l, j) + F(l, j) * A(m, j)+ A1(l, m, j) + A2(l, m, j));
        }
      }
      // Block(1,2)
      for(m = d; m < 2*d; m++){
        out(l, m) = -as_scalar(U.row(l) * (U.row(m - d) % Pi.row(m - d)).t());
      }
      //Block (1,3) + Block (2,3)
      for(m = 2 * d; m < q; m++){
        zm = z(m - 2 * d);
        wm = w(m - 2 * d);
        out(l, m) = -as_scalar(U.row(l) * tildeU.col(m - 2 * d));
        for(j = 0; j < d; j++){
          out(l + d, m) -= U(l, j) * (U(wm, j) * (F(zm, j) * F(l, j) * loewner(j, j)/2 + F(zm, j) * A(l, j) + A(zm, j) * F(l, j) + A1(l, zm, j) + A2(l, zm, j)) +
                                      U(zm, j) * (F(wm, j) * F(l, j) * loewner(j, j)/2 + F(wm, j) * A(l, j) + A(wm, j) * F(l, j) + A1(l, wm, j) + A2(l, wm, j)));
        }
      }
    }
    //Block (3,3)
    for(l = 2 * d; l < q; l++){
      zl = z(l - 2 * d);
      wl = w(l - 2 * d);
      for(m = l; m < q; m++){
        zm = z(m - 2 * d);
        wm = w(m - 2 * d);
        for(j = 0; j < d; j++){
          out(l, m) -= (U(wl, j) * (U(wm, j) * (F(zm, j) * F(zl, j) * loewner(j, j)/2 + F(zm, j) * A(zl, j) + A(zm, j) * F(zl, j) + A1(zl, zm, j) + A2(zm, zl, j)) +
                                    U(zm, j) * (F(wm, j) * F(zl, j) * loewner(j, j)/2 + F(wm, j) * A(zl, j) + A(wm, j) * F(zl, j) + A1(zl, wm, j) + A2(wm, zl, j)))+
                        U(zl, j) * (U(wm, j) * (F(zm, j) * F(wl, j) * loewner(j, j)/2 + F(zm, j) * A(wl, j) + A(zm, j) * F(wl, j) + A1(wl, zm, j) + A2(zm, wl, j))+
                                    U(zm, j) * (F(wm, j) * F(wl, j) * loewner(j, j)/2 + F(wm, j) * A(wl, j) + A(wm, j) * F(wl, j) + A1(wl, wm, j) + A2(wm, wl, j))));
        }
      }
    }

    //Saving the nonzero elements in a vector or a matrix according if they involve the intercepts or not
    it_b1     = ib1_eta.begin();
    it_b1_end = ib1_eta.end();
    it_b3     = ib3.begin();
    it_b3_end = ib3.end();

    // Full + partial
    c1 = 0;
    for(; it_b1 !=  it_b1_end; ++it_b1){
      ik1 = b1_eta[(*it_b1)];
      for(k = 0; k < ik1.length(); k++){
        d2l(i, c1) = out((*it_b1), ik1(k));
        c1 += 1;
      }
    }
    // intercepts
    c1 = 0;
    for(; it_b3 !=  it_b3_end; ++it_b3){
      ik1 = b3[(*it_b3)];
      for(k = 0; k < ik1.length(); k++){
        d2l_v(c1) += out((*it_b3), ik1(k));
        c1 += 1;
      }
    }
    out.fill(0); //put to zero out
  }
  return(1.0);
}
