#' A pool of internal functions used by the family
#'
#' @description Such function masks some functions used internally (please do not modify this function)
#'
#' @return a list of building functions
#' @export
internal <- function(){
  return(list(# Common functions (between mcd and logm):
              "gamlss.gH_scm" = gamlss.gH_scm,
              "mat2vec" = mat2vec,
              "Sigma_mat" = Sigma_mat,
              "idx_zwGt" = idx_zwGt,
              "idxHess_no0" = idxHess_no0,
              "formula_init" = formula_init,
              "aux_idx" = aux_idx,
              "d1_beta" = d1_beta,
              "d2_beta"  = d2_beta,
              "labTh" = labTh,
              # mcd functions:
              "ll_mcd" =  ll_mcd,
              "res_dev_mcd" = res_dev_mcd,
              "pred_mcd" = pred_mcd,
              "jacobian_mcd" = jacobian_mcd,
              "mcd_LD" = mcd_LD,
              "dHess_drho" = dHess_drho,
              "d3_mcd_eta" = d3_mcd_eta,
              "aux_idx_l3" = aux_idx_l3,
              "il3" = il3,
              "il3_no0_mcd" = il3_no0_mcd ,
              # logm functions:
              "ll_logm" =  ll_logm,
              "res_dev_logm" = res_dev_logm,
              "pred_logm" = pred_logm,
              "jacobian_logm" = jacobian_logm,
              "logM_Sigma" = logM_Sigma))
}




