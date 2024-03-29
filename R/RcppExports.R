# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

d1_beta <- function(X, eta, y, jj, K, lb, l1, l1_l, ig, z, w, Gm, param) {
    .Call(`_SCM_d1_beta`, X, eta, y, jj, K, lb, l1, l1_l, ig, z, w, Gm, param)
}

d1_logm_eta <- function(eta, y, d1l, z, w) {
    .Call(`_SCM_d1_logm_eta`, eta, y, d1l, z, w)
}

d1_mcd_eta <- function(eta, y, d1l, z, w, Gm) {
    .Call(`_SCM_d1_mcd_eta`, eta, y, d1l, z, w, Gm)
}

d2_beta <- function(X, eta, y, jj, K, lbb, l2, l2_v, l2_l, l2_v_l, ig, z, w, G, t, b1_eta, b1, b2, b3, ib1_eta, ib3, l2_el, l2_el2, param) {
    .Call(`_SCM_d2_beta`, X, eta, y, jj, K, lbb, l2, l2_v, l2_l, l2_v_l, ig, z, w, G, t, b1_eta, b1, b2, b3, ib1_eta, ib3, l2_el, l2_el2, param)
}

d2_logm_eta <- function(eta, y, d2l, d2l_v, z, w, b1_eta, b3, ib1_eta, ib3) {
    .Call(`_SCM_d2_logm_eta`, eta, y, d2l, d2l_v, z, w, b1_eta, b3, ib1_eta, ib3)
}

d2_mcd_eta <- function(eta, y, d2l, d2l_v, z, w, G, t, b1_eta, b3, ib1_eta, ib3) {
    .Call(`_SCM_d2_mcd_eta`, eta, y, d2l, d2l_v, z, w, G, t, b1_eta, b3, ib1_eta, ib3)
}

d3_mcd_eta <- function(eta, y, res, z, w, G, t, h, h2, h3, idx3_1, idx3_2, idx3_3, idx3_4, idx3_5, idx3_6, idx_neq0) {
    .Call(`_SCM_d3_mcd_eta`, eta, y, res, z, w, G, t, h, h2, h3, idx3_1, idx3_2, idx3_3, idx3_4, idx3_5, idx3_6, idx_neq0)
}

dHess_drho <- function(X, eta, y, jj, K, l3, l3_l, ig, d1b, a, a_l, d1H, fh, z, w, G, t, idx_l3, idx_neq0, idx_jkq) {
    .Call(`_SCM_dHess_drho`, X, eta, y, jj, K, l3, l3_l, ig, d1b, a, a_l, d1H, fh, z, w, G, t, idx_l3, idx_neq0, idx_jkq)
}

idx_zwGt <- function(d, z, w, G, t) {
    .Call(`_SCM_idx_zwGt`, d, z, w, G, t)
}

jacobian_logm <- function(eta, res, d, S_row, S_col, cor_flag) {
    .Call(`_SCM_jacobian_logm`, eta, res, d, S_row, S_col, cor_flag)
}

jacobian_mcd <- function(eta, res, d, S_r, S_c, rc_idx_s, rc_idx_t, cor_flag) {
    .Call(`_SCM_jacobian_mcd`, eta, res, d, S_r, S_c, rc_idx_s, rc_idx_t, cor_flag)
}

ll_mcd <- function(eta, y) {
    .Call(`_SCM_ll_mcd`, eta, y)
}

ll_logm <- function(eta, y) {
    .Call(`_SCM_ll_logm`, eta, y)
}

logM_Sigma <- function(x, d) {
    .Call(`_SCM_logM_Sigma`, x, d)
}

logm_decomposition <- function(X) {
    .Call(`_SCM_logm_decomposition`, X)
}

lt_inversion <- function(X) {
    .Call(`_SCM_lt_inversion`, X)
}

mcd_LD <- function(x, d) {
    .Call(`_SCM_mcd_LD`, x, d)
}

mcd_Sigma <- function(L, D, d) {
    .Call(`_SCM_mcd_Sigma`, L, D, d)
}

mcd_decomposition <- function(X) {
    .Call(`_SCM_mcd_decomposition`, X)
}

precision <- function(X) {
    .Call(`_SCM_precision`, X)
}

pred_logm <- function(eta, pred, d, cor_flag) {
    .Call(`_SCM_pred_logm`, eta, pred, d, cor_flag)
}

pred_mcd <- function(eta, pred, d, cor_flag) {
    .Call(`_SCM_pred_mcd`, eta, pred, d, cor_flag)
}

res_dev_logm <- function(eta, y, resD) {
    .Call(`_SCM_res_dev_logm`, eta, y, resD)
}

res_dev_mcd <- function(eta, y, res) {
    .Call(`_SCM_res_dev_mcd`, eta, y, res)
}

