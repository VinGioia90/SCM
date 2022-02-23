#' Building formula for the covariance model
#' @description This function return the formula for the covariance model, namely the formula for the elements of the modified cholesky decomposition
#'
#' @param data_boost The gradient boosting object
#' @param res The summary of the gradient boosting
#' @param stop_elem The number of elements considered for modelling purposes
#'
#' @return The formula (a list) for the covariance model
#' @export
#' @importFrom stats as.formula
#'
#' @examples
#'
formula_mcd <- function(data_boost, res, stop_elem=20){
  Theta_formula <- formula_mcd_int(data_boost$d)
  name_pred  = unlist(data_boost$effects)
  data_boost_unique <- unique(data_boost$idx)
  data_boost_unique_val <- rep(NA, stop_elem)
  data_boost_unique_form <- rep(NA, stop_elem)

  for(i in 1:stop_elem){
    data_boost_unique_val[i]<- res[i, "ele_rows"]
    data_boost_unique_form[i] <-data_boost$effects[[res[i,"label"]]]
  }
  idx_uni_value<-unique(data_boost_unique_val)

  for(i in 1:length(idx_uni_value)){
    Theta_formula[[idx_uni_value[i]]] <- as.formula(paste("~1+",paste(data_boost_unique_form[which(data_boost_unique_val %in% idx_uni_value[i])], collapse="+")))
  }
  return(Theta_formula)
}
