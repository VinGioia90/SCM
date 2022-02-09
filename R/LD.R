#' LD decomposition from mcd
#'
#' @description ...
#' @param x The additive predictor vector
#' @param d The dimension of the outcome
#'
#' @return The LD matrix via mcd
#' @export
#'
#' @examples


LD <- function(x,d){
  T_mat <- matrix(0, d, d)
  L <- matrix(0, d, d)
  D <- matrix(0, d, d)
  count <- d + 1
  for(j in 1:d){
    D[j,j] <- exp(0.5 * x[j + d])
    T_mat[j,j] <- 1
    if(j > 1){
      for(k in 1:(j - 1)){
        T_mat[j, k] <- x[count + d]
        count <- count + 1
      }
    }
  }

  L <- lowertri_inv(T_mat)
  return(list(L = L, D = D))
}
