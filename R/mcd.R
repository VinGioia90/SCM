#' Modified Cholesky decomposition 2
#'
#' @description The mcd function allows to obtain the unconstrained elements of the modified Cholesky decomposition of the precision matrix
#' @param X A covariance matrix
#'
#' @return A matrix with the unconstrained elements of the modified Cholesky decomposition of the precision matrix
#' @export
#'
#' @examples

mcd <- function(X){
 res <- mcd_decomposition(X)
 return(res)
}

