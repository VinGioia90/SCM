sel_elem <- function(d){
  pred_elem <- function(j,k){

    CovCor_mat <- function(d){
      Cm <- matrix(0, d, d)
      diag(Cm) <- 1:d
      count <- d + 1
      for(i in 2:d){
        for(j in 1:(i - 1)){
          Cm[i,j] <- Cm[j,i] <- count
          count <- count + 1
        }
      }
      return(Cm+d)
    }

    return(CovCor_mat(d)[j,k])
  }
  return(pred_elem)
}
