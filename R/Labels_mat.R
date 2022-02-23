#' Labels of mcd decomposition
#' @description This function return the labels of the unconstrained mcd elements
#'
#' @param d The dimension of the outcome
#'
#' @return A matrix with the labels of the uncostrained mcd elements
#' @export
#'
#' @examples
#'
l_mat <- function(d){
  mat <- matrix(NA, d, d)
  for(i in 1:d){
    for(j in 1:i){
      if(i == j){
        mat[i, i] <- paste0("logD2[", i, ",", i, "]")
      } else {
        mat[i, j] <- mat[j, i] <- paste0("T[", i, ",", j, "]")
      }
    }
  }
  return(mat = mat)
}
