#' Index Zero components
#'
#' @description The G matrix
#' @param no_eta The dimension of the outcome
#' @param z auxiliary index vector z
#' @param w auxiliary index vector w
#' @param Gmat auxiliary index matrix G (not used here)
#' @param t auxiliary index vector t (not used here)
#'
#' @return list of zero elements for each row
#' @export
#'
#' @examples
index0s <- function(no_eta, z, w, Gmat, t){
  idx_jk_no0 <- list()
  d <- -3/2 + sqrt(9/4 + 2 * no_eta)

  z1 <- z + 1
  w1 <- w + 1

  for(j in 1 : no_eta){
    idx_jk_no0[[j]] <- rep(0, no_eta)
    for(k in  j : no_eta){
      if(j <= d & k<= d){ idx_jk_no0[[j]][k]<-k}
      if(j <= d & k> d & k<=2*d){if(k-j-d>=0) {idx_jk_no0[[j]][k]<-k}}
      if(j <= d & k> 2*d){if(w1[k-2*d] >= j | z1[k-2*d] == j) {idx_jk_no0[[j]][k]<-k} }
      if(j>d & j <=2*d & k>d & k<=2*d & j==k){idx_jk_no0[[j]][k]<-k }
      if(j>d & j <=2*d & k>2*d){if(w1[k-2*d]==j-d){idx_jk_no0[[j]][k]<-k} }
      if(j>2*d & k>2*d){if(w1[j-2*d]==w1[k-2*d]){idx_jk_no0[[j]][k]<-k} }
    }
  }
  return(lapply(1 : no_eta, function( x ) which( idx_jk_no0[[x]] != 0) - 1 ))
}
