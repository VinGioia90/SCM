#' Matrix logarithm
#'
#' @description The logm function allows to obtain the unconstrained elements of the matrix logarithm transform
#' @param X A covariance matrix
#'
#' @return The logarithm of the covariance matrix
#' @export
#'
#' @examples
logm <- function(X){
  res <- logm_decomposition(X)
  return(res)
}
