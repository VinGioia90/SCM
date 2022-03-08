#' Negative log-likelihood logM
#'
#' @description The function returns the negative log-likelihood according to the logM parameterization
#' @param eta The linear predictor, a matrix of dimension n x (d+d x (d+1)/2))
#' @param y The outcome, that is a matrix n x d 
#'
#' @return The function returns the negative log-likelihood according to the logM parameterization
#' @export
#'
#' @examples
nll_logm <- function(eta, y){
  n <- nrow(y)
  d <- ncol(y)
  res <- -ll_logm(eta, y) + 0.5 * n * d * log(2 * pi)
  return(res)
}



