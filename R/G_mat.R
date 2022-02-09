#' Matrice G
#'
#' @description The G matrix
#' @param d The dimension of the outcome
#'
#' @return The G matrix
#' @export
#'
#' @examples

G_mat <- function(d){
  G <- matrix(0, d, d)
  diag(G) <- 1:d
  count <- d + 1
  for(i in 2:d){
   for(j in 1:(i - 1)){
    G[i,j] <- count
    count <- count + 1
   }
  }
  return(G)
}
