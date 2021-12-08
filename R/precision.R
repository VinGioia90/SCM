#' Precision matrix
#'
#' @description The inverse function allows to obtain the precision matrix
#' @param X A covariance matrix
#'
#' @return The precision matrix
#' @export
#'
#' @examples
inverse <- function(X){
  return(solve(X))
}
