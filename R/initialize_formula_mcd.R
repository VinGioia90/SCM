#' covariance model initialization formula
#' @description This function return the covariance model initialization formula
#'
#' @param d The dimension of the outcome
#'
#' @return The initializing formula (a list) for the covariance model
#' @export
#'
#' @examples
#'
formula_mcd_int <- function(d){
  Theta_formula_int <-  list()
  for(ii in 1:(d*(d + 1)/2)){
    Theta_formula_int[[ii]] <- as.formula("~ 1")
  }
  return(Theta_formula_int)
}
