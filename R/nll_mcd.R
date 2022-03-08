#' Negative log-likelihood MCD
#'
#' @description The function returns the negative log-likelihood according to the MCD parameterization
#' @param eta The linear predictor, a matrix of dimension n x (d+d x (d+1)/2))
#' @param y The outcome, that is a matrix n x d 
#'
#' @return The function returns the negative log-likelihood according to the MCD parameterization
#' @export
#'
#' @examples
nll_mcd <- function(eta, y){
  n <- nrow(y)
  d <- ncol(y)
  res <- -ll_mcd(eta, y) + 0.5 * n * d * log(2 * pi)
  return(res)
}



