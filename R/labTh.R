labTh <- function(d, k){
  d<-4
  Thmat <- matrix(NA, d, d)
  for(i in 1:d){
    Thmat[i,i] <- paste0("Th_",i,i)
    if(i >1){
      for(j in 1:(i-1)){
        Thmat[j,i] <- paste0("Th_",j,i)
      }
    }
  }
  Thvec <- c(diag(Thmat),Thmat[upper.tri(Thmat, diag=FALSE)])
  return(Thvec[k-d])
}


