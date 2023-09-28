labTh <- function(d, k){
  # This function  creates the labels used in the modified print summary function
  Thmat <- matrix(NA, d, d)
  d_len <- nchar(d)
  if(d_len == 1){
    for ( i in 1 : d ) {
      Thmat[i, i] <- paste0("Th_", i, i)
      if ( i > 1 ) {
        for(j in 1:(i-1)){
          Thmat[j,i] <- paste0("Th_", j, i)
        }
      }
    }
  } else {
    for ( i in 1 : d ) {
      i_len <- nchar(i)
      if(i_len < d_len){
        Thmat[i, i] <- paste0("Th_", rep("0",d_len-i_len), i, rep("0",d_len-i_len),i)
        if ( i > 1 ) {
          for(j in 1:(i-1)){
            j_len <- nchar(j)
            if(j_len < d_len){
              Thmat[j,i] <- paste0("Th_", rep("0",d_len-j_len), j, rep("0",d_len-i_len),i)
            } else {
              Thmat[j,i] <- paste0("Th_",j, rep("0",d_len-j_len),i)
            }
          }
        }
      } else {
        Thmat[i, i] <- paste0("Th_", i, i)
        if ( i > 1 ) {
          for(j in 1:(i-1)){
            j_len <- nchar(j)
            if(j_len < d_len){
              Thmat[j,i] <- paste0("Th_", rep("0",d_len-j_len), j,i)
            } else {
              Thmat[j,i] <- paste0("Th_",j,i)
            }
          }
        }
      }

    }
  }

  Thvec <- c(diag(Thmat), Thmat[upper.tri(Thmat, diag = FALSE)])
  return(Thvec[k - d])
}


