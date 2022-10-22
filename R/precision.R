#' Precision matrix (inverse of the covariance matrix)
#'
#' @description The inverse function allows to obtain the precision matrix
#' @param X A covariance matrix
#'
#' @return The precision matrix
#' @export
inverse <- function(X){
 res <- precision(X)
 return(res)
}

