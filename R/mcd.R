#' Modified Cholesky decomposition
#'
#' @description The mcd function allows to obtain the unconstrained elements of the modified Cholesky decomposition of the precision matrix
#' @param X A covariance matrix
#'
#' @return A matrix with the unconstrained elements of the modified Cholesky decomposition of the precision matrix
#' @export
mcd <- function(X){
 #such function returns a matrix containing the elements of log(D^2) on diagonal and the elements of T out the diagonal
 res <- mcd_decomposition(X) #call to cpp function
 return(res)
}

