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
