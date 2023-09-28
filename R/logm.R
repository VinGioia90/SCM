#' Matrix logarithm
#'
#' @description The logm function allows to obtain the unconstrained elements of the matrix logarithm transform
#' @param X A covariance matrix
#'
#' @return A matrix with the unconstrained elements resulting from the matrix logarithm
#' @export
logm <- function(X){
  res <- logm_decomposition(X)
  return(res)
}
