#' Lower triangular matrix inversion
#'
#' @description ...
#' @param X A lower triangular matrix
#'
#' @return the inverse of a lower triangular matrix
#' @export
#'
#' @examples

lowertri_inv <- function(X){ #X is a lower triangular matrix
   d <- ncol(X)
   res <- matrix(0, d, d)
    for(j in 1:d){
    res[j, j] <- 1
    if(j > 1){
      for(k in 1:(j - 1)){
        s <- 0
        for(l in k:(j - 1)){
          s <- s + X[j, l] * res[l, k]
        }
        res[j, k] <- -s
      }
    }
  }
  return(res)
}
