#' Modified Cholesky decomposition
#'
#' @description The mcd function allows to obtain the unconstrained elements of the modified Cholesky decomposition of the precision matrix
#' @param X A covariance matrix
#'
#' @return A matrix with the unconstrained elements of the modified Cholesky decomposition of the precision matrix
#' @export
#'
#' @examples
mcd <- function(X){
  d <- ncol(X)
  res <- matrix(0, d, d)
  chol_dec <- chol(X)
  D <- diag(chol_dec)
  L <- matrix(0, d, d)
  diag(L) <- 1
  for(j in 2:d){
    for(k in 1:(j - 1)){
      L[j, k] <- chol_dec[k, j]/D[k]
    }
  }

  for(j in 1:d){
    res[j, j] <- 1
    if(j > 1){
      for(k in 1:(j - 1)){
        s <- 0
        for(l in k:(j - 1)){
          s <- s + L[j, l] * res[l, k]
        }
        res[j, k] <- -s
      }
    }
  }
  res <- res + t(res)
  diag(res) <- log(D^2)
  return(res)
}
