// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// d1_beta
double d1_beta(arma::mat& X, arma::mat& eta, arma::mat& y, Rcpp::List& jj, uint32_t& K, arma::vec& lb, arma::mat& l1, arma::mat& l1_l, arma::uvec& ig, arma::vec& z, arma::vec& w, arma::mat& Gm, uint32_t& param);
RcppExport SEXP _SCM_d1_beta(SEXP XSEXP, SEXP etaSEXP, SEXP ySEXP, SEXP jjSEXP, SEXP KSEXP, SEXP lbSEXP, SEXP l1SEXP, SEXP l1_lSEXP, SEXP igSEXP, SEXP zSEXP, SEXP wSEXP, SEXP GmSEXP, SEXP paramSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type l1(l1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type l1_l(l1_lSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type ig(igSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Gm(GmSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type param(paramSEXP);
    rcpp_result_gen = Rcpp::wrap(d1_beta(X, eta, y, jj, K, lb, l1, l1_l, ig, z, w, Gm, param));
    return rcpp_result_gen;
END_RCPP
}
// d1_logm_eta
double d1_logm_eta(const arma::mat& eta, const arma::mat& y, arma::mat& d1l, arma::vec& z, arma::vec& w);
RcppExport SEXP _SCM_d1_logm_eta(SEXP etaSEXP, SEXP ySEXP, SEXP d1lSEXP, SEXP zSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type d1l(d1lSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(d1_logm_eta(eta, y, d1l, z, w));
    return rcpp_result_gen;
END_RCPP
}
// d1_mcd_eta
double d1_mcd_eta(const arma::mat& eta, const arma::mat& y, arma::mat& d1l, arma::vec& z, arma::vec& w, arma::mat& Gm);
RcppExport SEXP _SCM_d1_mcd_eta(SEXP etaSEXP, SEXP ySEXP, SEXP d1lSEXP, SEXP zSEXP, SEXP wSEXP, SEXP GmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type d1l(d1lSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Gm(GmSEXP);
    rcpp_result_gen = Rcpp::wrap(d1_mcd_eta(eta, y, d1l, z, w, Gm));
    return rcpp_result_gen;
END_RCPP
}
// d2_beta
double d2_beta(arma::mat& X, arma::mat& eta, arma::mat& y, Rcpp::List& jj, uint32_t& K, arma::mat& lbb, arma::mat& l2, arma::vec& l2_v, arma::mat& l2_l, arma::vec& l2_v_l, arma::uvec& ig, arma::vec& z, arma::vec& w, arma::mat& G, arma::vec& t, Rcpp::List& b1_eta, Rcpp::List& b1, Rcpp::List& b2, Rcpp::List& b3, arma::vec& ib1_eta, arma::vec& ib3, arma::vec& l2_el, arma::vec& l2_el2, uint32_t param);
RcppExport SEXP _SCM_d2_beta(SEXP XSEXP, SEXP etaSEXP, SEXP ySEXP, SEXP jjSEXP, SEXP KSEXP, SEXP lbbSEXP, SEXP l2SEXP, SEXP l2_vSEXP, SEXP l2_lSEXP, SEXP l2_v_lSEXP, SEXP igSEXP, SEXP zSEXP, SEXP wSEXP, SEXP GSEXP, SEXP tSEXP, SEXP b1_etaSEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP b3SEXP, SEXP ib1_etaSEXP, SEXP ib3SEXP, SEXP l2_elSEXP, SEXP l2_el2SEXP, SEXP paramSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type lbb(lbbSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type l2_v(l2_vSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type l2_l(l2_lSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type l2_v_l(l2_v_lSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type ig(igSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type t(tSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type b1_eta(b1_etaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type b3(b3SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type ib1_eta(ib1_etaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type ib3(ib3SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type l2_el(l2_elSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type l2_el2(l2_el2SEXP);
    Rcpp::traits::input_parameter< uint32_t >::type param(paramSEXP);
    rcpp_result_gen = Rcpp::wrap(d2_beta(X, eta, y, jj, K, lbb, l2, l2_v, l2_l, l2_v_l, ig, z, w, G, t, b1_eta, b1, b2, b3, ib1_eta, ib3, l2_el, l2_el2, param));
    return rcpp_result_gen;
END_RCPP
}
// d2_logm_eta
double d2_logm_eta(const arma::mat& eta, const arma::mat& y, arma::mat& d2l, arma::vec& d2l_v, arma::vec& z, arma::vec& w, Rcpp::List& b1_eta, Rcpp::List& b3, arma::vec& ib1_eta, arma::vec& ib3);
RcppExport SEXP _SCM_d2_logm_eta(SEXP etaSEXP, SEXP ySEXP, SEXP d2lSEXP, SEXP d2l_vSEXP, SEXP zSEXP, SEXP wSEXP, SEXP b1_etaSEXP, SEXP b3SEXP, SEXP ib1_etaSEXP, SEXP ib3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type d2l(d2lSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type d2l_v(d2l_vSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type b1_eta(b1_etaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type b3(b3SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type ib1_eta(ib1_etaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type ib3(ib3SEXP);
    rcpp_result_gen = Rcpp::wrap(d2_logm_eta(eta, y, d2l, d2l_v, z, w, b1_eta, b3, ib1_eta, ib3));
    return rcpp_result_gen;
END_RCPP
}
// d2_mcd_eta
double d2_mcd_eta(const arma::mat& eta, const arma::mat& y, arma::mat& d2l, arma::vec& d2l_v, arma::vec& z, arma::vec& w, arma::mat& G, arma::vec& t, Rcpp::List& b1_eta, Rcpp::List& b3, arma::vec& ib1_eta, arma::vec& ib3);
RcppExport SEXP _SCM_d2_mcd_eta(SEXP etaSEXP, SEXP ySEXP, SEXP d2lSEXP, SEXP d2l_vSEXP, SEXP zSEXP, SEXP wSEXP, SEXP GSEXP, SEXP tSEXP, SEXP b1_etaSEXP, SEXP b3SEXP, SEXP ib1_etaSEXP, SEXP ib3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type d2l(d2lSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type d2l_v(d2l_vSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type t(tSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type b1_eta(b1_etaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type b3(b3SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type ib1_eta(ib1_etaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type ib3(ib3SEXP);
    rcpp_result_gen = Rcpp::wrap(d2_mcd_eta(eta, y, d2l, d2l_v, z, w, G, t, b1_eta, b3, ib1_eta, ib3));
    return rcpp_result_gen;
END_RCPP
}
// d3_mcd_eta
double d3_mcd_eta(const arma::mat& eta, const arma::mat& y, arma::mat& res, arma::vec& z, arma::vec& w, arma::mat& G, arma::vec& t, Rcpp::IntegerVector& h, Rcpp::IntegerVector& h2, Rcpp::IntegerVector& h3, Rcpp::List& idx3_1, Rcpp::List& idx3_2, Rcpp::List& idx3_3, Rcpp::List& idx3_4, Rcpp::List& idx3_5, Rcpp::List& idx3_6, Rcpp::List& idx_neq0);
RcppExport SEXP _SCM_d3_mcd_eta(SEXP etaSEXP, SEXP ySEXP, SEXP resSEXP, SEXP zSEXP, SEXP wSEXP, SEXP GSEXP, SEXP tSEXP, SEXP hSEXP, SEXP h2SEXP, SEXP h3SEXP, SEXP idx3_1SEXP, SEXP idx3_2SEXP, SEXP idx3_3SEXP, SEXP idx3_4SEXP, SEXP idx3_5SEXP, SEXP idx3_6SEXP, SEXP idx_neq0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type res(resSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type t(tSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type h(hSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type h3(h3SEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type idx3_1(idx3_1SEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type idx3_2(idx3_2SEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type idx3_3(idx3_3SEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type idx3_4(idx3_4SEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type idx3_5(idx3_5SEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type idx3_6(idx3_6SEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type idx_neq0(idx_neq0SEXP);
    rcpp_result_gen = Rcpp::wrap(d3_mcd_eta(eta, y, res, z, w, G, t, h, h2, h3, idx3_1, idx3_2, idx3_3, idx3_4, idx3_5, idx3_6, idx_neq0));
    return rcpp_result_gen;
END_RCPP
}
// dHess_drho
double dHess_drho(arma::mat& X, arma::mat& eta, arma::mat& y, Rcpp::List& jj, uint32_t& K, arma::mat& l3, arma::mat& l3_l, arma::uvec& ig, arma::mat& d1b, arma::vec& a, arma::vec& a_l, arma::vec& d1H, arma::mat& fh, arma::vec& z, arma::vec& w, arma::mat& G, arma::vec& t, Rcpp::List& idx_l3, Rcpp::List& idx_neq0, Rcpp::List& idx_jkq);
RcppExport SEXP _SCM_dHess_drho(SEXP XSEXP, SEXP etaSEXP, SEXP ySEXP, SEXP jjSEXP, SEXP KSEXP, SEXP l3SEXP, SEXP l3_lSEXP, SEXP igSEXP, SEXP d1bSEXP, SEXP aSEXP, SEXP a_lSEXP, SEXP d1HSEXP, SEXP fhSEXP, SEXP zSEXP, SEXP wSEXP, SEXP GSEXP, SEXP tSEXP, SEXP idx_l3SEXP, SEXP idx_neq0SEXP, SEXP idx_jkqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type l3(l3SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type l3_l(l3_lSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type ig(igSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type d1b(d1bSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type a_l(a_lSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type d1H(d1HSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type fh(fhSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type t(tSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type idx_l3(idx_l3SEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type idx_neq0(idx_neq0SEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type idx_jkq(idx_jkqSEXP);
    rcpp_result_gen = Rcpp::wrap(dHess_drho(X, eta, y, jj, K, l3, l3_l, ig, d1b, a, a_l, d1H, fh, z, w, G, t, idx_l3, idx_neq0, idx_jkq));
    return rcpp_result_gen;
END_RCPP
}
// idx_zwGt
double idx_zwGt(uint32_t d, IntegerVector& z, IntegerVector& w, IntegerMatrix& G, IntegerVector& t);
RcppExport SEXP _SCM_idx_zwGt(SEXP dSEXP, SEXP zSEXP, SEXP wSEXP, SEXP GSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< uint32_t >::type d(dSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type z(zSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type G(GSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(idx_zwGt(d, z, w, G, t));
    return rcpp_result_gen;
END_RCPP
}
// jacobian_logm
double jacobian_logm(arma::mat& eta, arma::mat& res, uint32_t& d, uint32_t& S_row, uint32_t& S_col, uint32_t& cor_flag);
RcppExport SEXP _SCM_jacobian_logm(SEXP etaSEXP, SEXP resSEXP, SEXP dSEXP, SEXP S_rowSEXP, SEXP S_colSEXP, SEXP cor_flagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type res(resSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type S_row(S_rowSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type S_col(S_colSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type cor_flag(cor_flagSEXP);
    rcpp_result_gen = Rcpp::wrap(jacobian_logm(eta, res, d, S_row, S_col, cor_flag));
    return rcpp_result_gen;
END_RCPP
}
// jacobian_mcd
double jacobian_mcd(arma::mat& eta, arma::mat& res, uint32_t& d, uint32_t& S_r, uint32_t& S_c, arma::vec rc_idx_s, arma::vec& rc_idx_t, uint32_t& cor_flag);
RcppExport SEXP _SCM_jacobian_mcd(SEXP etaSEXP, SEXP resSEXP, SEXP dSEXP, SEXP S_rSEXP, SEXP S_cSEXP, SEXP rc_idx_sSEXP, SEXP rc_idx_tSEXP, SEXP cor_flagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type res(resSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type S_r(S_rSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type S_c(S_cSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rc_idx_s(rc_idx_sSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type rc_idx_t(rc_idx_tSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type cor_flag(cor_flagSEXP);
    rcpp_result_gen = Rcpp::wrap(jacobian_mcd(eta, res, d, S_r, S_c, rc_idx_s, rc_idx_t, cor_flag));
    return rcpp_result_gen;
END_RCPP
}
// ll_mcd
double ll_mcd(const arma::mat& eta, const arma::mat& y);
RcppExport SEXP _SCM_ll_mcd(SEXP etaSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(ll_mcd(eta, y));
    return rcpp_result_gen;
END_RCPP
}
// ll_logm
double ll_logm(arma::mat& eta, arma::mat& y);
RcppExport SEXP _SCM_ll_logm(SEXP etaSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(ll_logm(eta, y));
    return rcpp_result_gen;
END_RCPP
}
// logM_Sigma
arma::mat logM_Sigma(arma::rowvec& x, uint32_t& d);
RcppExport SEXP _SCM_logM_Sigma(SEXP xSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(logM_Sigma(x, d));
    return rcpp_result_gen;
END_RCPP
}
// logm_decomposition
arma::mat logm_decomposition(arma::mat& X);
RcppExport SEXP _SCM_logm_decomposition(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(logm_decomposition(X));
    return rcpp_result_gen;
END_RCPP
}
// lt_inversion
arma::mat lt_inversion(arma::mat& X);
RcppExport SEXP _SCM_lt_inversion(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(lt_inversion(X));
    return rcpp_result_gen;
END_RCPP
}
// mcd_LD
arma::mat mcd_LD(arma::rowvec& x, uint32_t& d);
RcppExport SEXP _SCM_mcd_LD(SEXP xSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(mcd_LD(x, d));
    return rcpp_result_gen;
END_RCPP
}
// mcd_Sigma
arma::mat mcd_Sigma(arma::mat& L, arma::mat& D, uint32_t& d);
RcppExport SEXP _SCM_mcd_Sigma(SEXP LSEXP, SEXP DSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(mcd_Sigma(L, D, d));
    return rcpp_result_gen;
END_RCPP
}
// mcd_decomposition
arma::mat mcd_decomposition(arma::mat& X);
RcppExport SEXP _SCM_mcd_decomposition(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(mcd_decomposition(X));
    return rcpp_result_gen;
END_RCPP
}
// precision
Rcpp::NumericMatrix precision(Rcpp::NumericMatrix& X);
RcppExport SEXP _SCM_precision(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(precision(X));
    return rcpp_result_gen;
END_RCPP
}
// pred_logm
double pred_logm(arma::mat& eta, arma::mat& pred, uint32_t& d, uint32_t& cor_flag);
RcppExport SEXP _SCM_pred_logm(SEXP etaSEXP, SEXP predSEXP, SEXP dSEXP, SEXP cor_flagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type cor_flag(cor_flagSEXP);
    rcpp_result_gen = Rcpp::wrap(pred_logm(eta, pred, d, cor_flag));
    return rcpp_result_gen;
END_RCPP
}
// pred_mcd
double pred_mcd(arma::mat& eta, arma::mat& pred, uint32_t& d, uint32_t& cor_flag);
RcppExport SEXP _SCM_pred_mcd(SEXP etaSEXP, SEXP predSEXP, SEXP dSEXP, SEXP cor_flagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type pred(predSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type cor_flag(cor_flagSEXP);
    rcpp_result_gen = Rcpp::wrap(pred_mcd(eta, pred, d, cor_flag));
    return rcpp_result_gen;
END_RCPP
}
// res_dev_logm
double res_dev_logm(arma::mat& eta, arma::mat& y, arma::mat& resD);
RcppExport SEXP _SCM_res_dev_logm(SEXP etaSEXP, SEXP ySEXP, SEXP resDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type resD(resDSEXP);
    rcpp_result_gen = Rcpp::wrap(res_dev_logm(eta, y, resD));
    return rcpp_result_gen;
END_RCPP
}
// res_dev_mcd
double res_dev_mcd(arma::mat eta, arma::mat& y, arma::mat& res);
RcppExport SEXP _SCM_res_dev_mcd(SEXP etaSEXP, SEXP ySEXP, SEXP resSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type res(resSEXP);
    rcpp_result_gen = Rcpp::wrap(res_dev_mcd(eta, y, res));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SCM_d1_beta", (DL_FUNC) &_SCM_d1_beta, 13},
    {"_SCM_d1_logm_eta", (DL_FUNC) &_SCM_d1_logm_eta, 5},
    {"_SCM_d1_mcd_eta", (DL_FUNC) &_SCM_d1_mcd_eta, 6},
    {"_SCM_d2_beta", (DL_FUNC) &_SCM_d2_beta, 24},
    {"_SCM_d2_logm_eta", (DL_FUNC) &_SCM_d2_logm_eta, 10},
    {"_SCM_d2_mcd_eta", (DL_FUNC) &_SCM_d2_mcd_eta, 12},
    {"_SCM_d3_mcd_eta", (DL_FUNC) &_SCM_d3_mcd_eta, 17},
    {"_SCM_dHess_drho", (DL_FUNC) &_SCM_dHess_drho, 20},
    {"_SCM_idx_zwGt", (DL_FUNC) &_SCM_idx_zwGt, 5},
    {"_SCM_jacobian_logm", (DL_FUNC) &_SCM_jacobian_logm, 6},
    {"_SCM_jacobian_mcd", (DL_FUNC) &_SCM_jacobian_mcd, 8},
    {"_SCM_ll_mcd", (DL_FUNC) &_SCM_ll_mcd, 2},
    {"_SCM_ll_logm", (DL_FUNC) &_SCM_ll_logm, 2},
    {"_SCM_logM_Sigma", (DL_FUNC) &_SCM_logM_Sigma, 2},
    {"_SCM_logm_decomposition", (DL_FUNC) &_SCM_logm_decomposition, 1},
    {"_SCM_lt_inversion", (DL_FUNC) &_SCM_lt_inversion, 1},
    {"_SCM_mcd_LD", (DL_FUNC) &_SCM_mcd_LD, 2},
    {"_SCM_mcd_Sigma", (DL_FUNC) &_SCM_mcd_Sigma, 3},
    {"_SCM_mcd_decomposition", (DL_FUNC) &_SCM_mcd_decomposition, 1},
    {"_SCM_precision", (DL_FUNC) &_SCM_precision, 1},
    {"_SCM_pred_logm", (DL_FUNC) &_SCM_pred_logm, 4},
    {"_SCM_pred_mcd", (DL_FUNC) &_SCM_pred_mcd, 4},
    {"_SCM_res_dev_logm", (DL_FUNC) &_SCM_res_dev_logm, 3},
    {"_SCM_res_dev_mcd", (DL_FUNC) &_SCM_res_dev_mcd, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_SCM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
