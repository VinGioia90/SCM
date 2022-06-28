#' Prediction Element Index
#'
#' @description Select the prediction element
#' @param d The dimension of the outcome
#'
#' @return The index selected
#' @export
#'
#' @examples
# idx <- sel_elem(4)
# idx(3,1)
sel_elem <- function(d){
  pred_elem <- function(j,k){

    CovCor_mat <- function(d){
      G <- matrix(0, d, d)
      diag(G) <- 1:d
      count <- d + 1
      for(i in 2:d){
        for(j in 1:(i - 1)){
          G[i,j] <- G[j,i] <- count
          count <- count + 1
        }
      }
      return(G+d)
    }

    return(CovCor_mat(d)[j,k])
  }
  return(pred_elem)
}
