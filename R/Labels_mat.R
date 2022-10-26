l_mat <- function(d){
  # Maybe not useful for the family: we used it in the gradient boosting initial version
  # However, it produces the labels for the T and logD2 matrices
  mat <- matrix(NA, d, d)
  diag(mat) <- paste0("logD2[", 1 : d, ",", 1 : d, "]")
  for (i in 2 : d) {
    for ( j in 1 : (i - 1) ) {
      mat[i, j] <- mat[j, i] <- paste0("T[", i, ",", j, "]")
    }
  }
  return(mat = mat)
}

