// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// d1_logm
double d1_logm(arma::mat& eta, arma::mat& y, arma::mat& res);
RcppExport SEXP _SCM_d1_logm(SEXP etaSEXP, SEXP ySEXP, SEXP resSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type res(resSEXP);
    rcpp_result_gen = Rcpp::wrap(d1_logm(eta, y, res));
    return rcpp_result_gen;
END_RCPP
}
// d1_mcd
double d1_mcd(NumericMatrix& eta, NumericMatrix& y, NumericMatrix& res);
RcppExport SEXP _SCM_d1_mcd(SEXP etaSEXP, SEXP ySEXP, SEXP resSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type res(resSEXP);
    rcpp_result_gen = Rcpp::wrap(d1_mcd(eta, y, res));
    return rcpp_result_gen;
END_RCPP
}
// d1_mcd_beta
double d1_mcd_beta(arma::mat& X, Rcpp::List& jj, uint32_t& K, arma::vec& lb, arma::mat& l1, arma::mat& eta, arma::mat& y, arma::vec& z, arma::vec& w, arma::mat& G);
RcppExport SEXP _SCM_d1_mcd_beta(SEXP XSEXP, SEXP jjSEXP, SEXP KSEXP, SEXP lbSEXP, SEXP l1SEXP, SEXP etaSEXP, SEXP ySEXP, SEXP zSEXP, SEXP wSEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lb(lbSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type l1(l1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(d1_mcd_beta(X, jj, K, lb, l1, eta, y, z, w, G));
    return rcpp_result_gen;
END_RCPP
}
// d1_mcd_eta
double d1_mcd_eta(const arma::mat& eta, const arma::mat& y, arma::mat& res, arma::vec& z, arma::vec& w, arma::mat& G);
RcppExport SEXP _SCM_d1_mcd_eta(SEXP etaSEXP, SEXP ySEXP, SEXP resSEXP, SEXP zSEXP, SEXP wSEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type res(resSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(d1_mcd_eta(eta, y, res, z, w, G));
    return rcpp_result_gen;
END_RCPP
}
// d2_logm
double d2_logm(arma::mat& eta, arma::mat& y, arma::mat& res);
RcppExport SEXP _SCM_d2_logm(SEXP etaSEXP, SEXP ySEXP, SEXP resSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type res(resSEXP);
    rcpp_result_gen = Rcpp::wrap(d2_logm(eta, y, res));
    return rcpp_result_gen;
END_RCPP
}
// d2_mcd
double d2_mcd(NumericMatrix& eta, NumericMatrix& y, NumericMatrix& res);
RcppExport SEXP _SCM_d2_mcd(SEXP etaSEXP, SEXP ySEXP, SEXP resSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type res(resSEXP);
    rcpp_result_gen = Rcpp::wrap(d2_mcd(eta, y, res));
    return rcpp_result_gen;
END_RCPP
}
// d2_mcd_beta
double d2_mcd_beta(arma::mat& X, arma::mat& eta, arma::mat& y, Rcpp::List& jj, uint32_t& K, Rcpp::List& idx_jk, arma::mat& lbb, arma::mat& l2, arma::mat& l2_last, arma::uvec& idx_g, arma::vec& z, arma::vec& w, arma::mat& G, arma::vec& t);
RcppExport SEXP _SCM_d2_mcd_beta(SEXP XSEXP, SEXP etaSEXP, SEXP ySEXP, SEXP jjSEXP, SEXP KSEXP, SEXP idx_jkSEXP, SEXP lbbSEXP, SEXP l2SEXP, SEXP l2_lastSEXP, SEXP idx_gSEXP, SEXP zSEXP, SEXP wSEXP, SEXP GSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type K(KSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type idx_jk(idx_jkSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type lbb(lbbSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type l2_last(l2_lastSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type idx_g(idx_gSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(d2_mcd_beta(X, eta, y, jj, K, idx_jk, lbb, l2, l2_last, idx_g, z, w, G, t));
    return rcpp_result_gen;
END_RCPP
}
// d2_mcd_eta
double d2_mcd_eta(const arma::mat& eta, const arma::mat& y, arma::mat& res, Rcpp::List& idx_jk, arma::vec& z, arma::vec& w, arma::mat& G, arma::vec& t);
RcppExport SEXP _SCM_d2_mcd_eta(SEXP etaSEXP, SEXP ySEXP, SEXP resSEXP, SEXP idx_jkSEXP, SEXP zSEXP, SEXP wSEXP, SEXP GSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type res(resSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type idx_jk(idx_jkSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(d2_mcd_eta(eta, y, res, idx_jk, z, w, G, t));
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
double jacobian_logm(Rcpp::NumericMatrix& eta, Rcpp::NumericMatrix& res, uint32_t& d, uint32_t& S_row, uint32_t& S_col, uint32_t& cor_flag);
RcppExport SEXP _SCM_jacobian_logm(SEXP etaSEXP, SEXP resSEXP, SEXP dSEXP, SEXP S_rowSEXP, SEXP S_colSEXP, SEXP cor_flagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type res(resSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type S_row(S_rowSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type S_col(S_colSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type cor_flag(cor_flagSEXP);
    rcpp_result_gen = Rcpp::wrap(jacobian_logm(eta, res, d, S_row, S_col, cor_flag));
    return rcpp_result_gen;
END_RCPP
}
// jacobian_mcd
double jacobian_mcd(Rcpp::NumericMatrix& eta, Rcpp::NumericMatrix& res, uint32_t& d, uint32_t& S_row, uint32_t& S_col, Rcpp::NumericVector rc_idx_s, Rcpp::NumericVector& rc_idx_t, uint32_t& cor_flag);
RcppExport SEXP _SCM_jacobian_mcd(SEXP etaSEXP, SEXP resSEXP, SEXP dSEXP, SEXP S_rowSEXP, SEXP S_colSEXP, SEXP rc_idx_sSEXP, SEXP rc_idx_tSEXP, SEXP cor_flagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type res(resSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type S_row(S_rowSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type S_col(S_colSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rc_idx_s(rc_idx_sSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type rc_idx_t(rc_idx_tSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type cor_flag(cor_flagSEXP);
    rcpp_result_gen = Rcpp::wrap(jacobian_mcd(eta, res, d, S_row, S_col, rc_idx_s, rc_idx_t, cor_flag));
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
// ll_mcd
double ll_mcd(NumericMatrix& eta, NumericMatrix& y);
RcppExport SEXP _SCM_ll_mcd(SEXP etaSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(ll_mcd(eta, y));
    return rcpp_result_gen;
END_RCPP
}
// logM_Sigma
Rcpp::NumericMatrix logM_Sigma(Rcpp::NumericVector& x, uint32_t& d);
RcppExport SEXP _SCM_logM_Sigma(SEXP xSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(logM_Sigma(x, d));
    return rcpp_result_gen;
END_RCPP
}
// logm_decomposition
Rcpp::NumericMatrix logm_decomposition(Rcpp::NumericMatrix& X);
RcppExport SEXP _SCM_logm_decomposition(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(logm_decomposition(X));
    return rcpp_result_gen;
END_RCPP
}
// lt_inversion
Rcpp::NumericMatrix lt_inversion(Rcpp::NumericMatrix& X);
RcppExport SEXP _SCM_lt_inversion(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(lt_inversion(X));
    return rcpp_result_gen;
END_RCPP
}
// mcd_LD
Rcpp::NumericMatrix mcd_LD(Rcpp::NumericVector& x, uint32_t& d);
RcppExport SEXP _SCM_mcd_LD(SEXP xSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(mcd_LD(x, d));
    return rcpp_result_gen;
END_RCPP
}
// mcd_Sigma
Rcpp::NumericMatrix mcd_Sigma(Rcpp::NumericMatrix& L, Rcpp::NumericMatrix& D, uint32_t& d);
RcppExport SEXP _SCM_mcd_Sigma(SEXP LSEXP, SEXP DSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type L(LSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type D(DSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(mcd_Sigma(L, D, d));
    return rcpp_result_gen;
END_RCPP
}
// mcd_decomposition
Rcpp::NumericMatrix mcd_decomposition(Rcpp::NumericMatrix& X);
RcppExport SEXP _SCM_mcd_decomposition(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type X(XSEXP);
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
double pred_logm(Rcpp::NumericMatrix& eta, Rcpp::NumericMatrix& res, uint32_t& d);
RcppExport SEXP _SCM_pred_logm(SEXP etaSEXP, SEXP resSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type res(resSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(pred_logm(eta, res, d));
    return rcpp_result_gen;
END_RCPP
}
// pred_mcd
double pred_mcd(Rcpp::NumericMatrix& eta, Rcpp::NumericMatrix& res, uint32_t& d, uint32_t& cor_flag);
RcppExport SEXP _SCM_pred_mcd(SEXP etaSEXP, SEXP resSEXP, SEXP dSEXP, SEXP cor_flagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type res(resSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type d(dSEXP);
    Rcpp::traits::input_parameter< uint32_t& >::type cor_flag(cor_flagSEXP);
    rcpp_result_gen = Rcpp::wrap(pred_mcd(eta, res, d, cor_flag));
    return rcpp_result_gen;
END_RCPP
}
// res_dev_logm
double res_dev_logm(Rcpp::NumericMatrix& eta, Rcpp::NumericMatrix& y, Rcpp::NumericMatrix& res);
RcppExport SEXP _SCM_res_dev_logm(SEXP etaSEXP, SEXP ySEXP, SEXP resSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type res(resSEXP);
    rcpp_result_gen = Rcpp::wrap(res_dev_logm(eta, y, res));
    return rcpp_result_gen;
END_RCPP
}
// res_dev_mcd
double res_dev_mcd(NumericMatrix& eta, NumericMatrix& y, NumericMatrix& res);
RcppExport SEXP _SCM_res_dev_mcd(SEXP etaSEXP, SEXP ySEXP, SEXP resSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type res(resSEXP);
    rcpp_result_gen = Rcpp::wrap(res_dev_mcd(eta, y, res));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SCM_d1_logm", (DL_FUNC) &_SCM_d1_logm, 3},
    {"_SCM_d1_mcd", (DL_FUNC) &_SCM_d1_mcd, 3},
    {"_SCM_d1_mcd_beta", (DL_FUNC) &_SCM_d1_mcd_beta, 10},
    {"_SCM_d1_mcd_eta", (DL_FUNC) &_SCM_d1_mcd_eta, 6},
    {"_SCM_d2_logm", (DL_FUNC) &_SCM_d2_logm, 3},
    {"_SCM_d2_mcd", (DL_FUNC) &_SCM_d2_mcd, 3},
    {"_SCM_d2_mcd_beta", (DL_FUNC) &_SCM_d2_mcd_beta, 14},
    {"_SCM_d2_mcd_eta", (DL_FUNC) &_SCM_d2_mcd_eta, 8},
    {"_SCM_idx_zwGt", (DL_FUNC) &_SCM_idx_zwGt, 5},
    {"_SCM_jacobian_logm", (DL_FUNC) &_SCM_jacobian_logm, 6},
    {"_SCM_jacobian_mcd", (DL_FUNC) &_SCM_jacobian_mcd, 8},
    {"_SCM_ll_logm", (DL_FUNC) &_SCM_ll_logm, 2},
    {"_SCM_ll_mcd", (DL_FUNC) &_SCM_ll_mcd, 2},
    {"_SCM_logM_Sigma", (DL_FUNC) &_SCM_logM_Sigma, 2},
    {"_SCM_logm_decomposition", (DL_FUNC) &_SCM_logm_decomposition, 1},
    {"_SCM_lt_inversion", (DL_FUNC) &_SCM_lt_inversion, 1},
    {"_SCM_mcd_LD", (DL_FUNC) &_SCM_mcd_LD, 2},
    {"_SCM_mcd_Sigma", (DL_FUNC) &_SCM_mcd_Sigma, 3},
    {"_SCM_mcd_decomposition", (DL_FUNC) &_SCM_mcd_decomposition, 1},
    {"_SCM_precision", (DL_FUNC) &_SCM_precision, 1},
    {"_SCM_pred_logm", (DL_FUNC) &_SCM_pred_logm, 3},
    {"_SCM_pred_mcd", (DL_FUNC) &_SCM_pred_mcd, 4},
    {"_SCM_res_dev_logm", (DL_FUNC) &_SCM_res_dev_logm, 3},
    {"_SCM_res_dev_mcd", (DL_FUNC) &_SCM_res_dev_mcd, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_SCM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
