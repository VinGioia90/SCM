#' Matrix logarithm
#'
#' @description The logm function allows to obtain the unconstrained elements of the matrix logarithm transform
#' @param X A covariance matrix
#'
#' @return A matrix the unconstrained elements of the matrix logarithm transform...
#' @export
#'
#' @examples
logm <- function(X){
  d <- ncol(X)
  res <- matrix(0, d, d)
  ed <- eigen(X)
  res <- ed$vectors%*%diag(log(ed$values))%*%t(ed$vectors)
  return(res)
}
