mat2vec <- function(d){
  C <- matrix(0, d, d)
  diag(C) <- 1 : d
  count <- d + 1
  for(i in 2 : d){
   for(j in 1 : (i - 1)){
    C[i, j] <- count
    count <- count + 1
   }
  }
  return(C)
}
