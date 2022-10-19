aux_idx_l3 <- function(d,z,w,G){
  no_eta <- d+d*(d+1)/2
  idx_jk <- list()
  #d <- -3/2 + sqrt(9/4 + 2 * no_eta)

  z <- z+1
  w <- w+1


  for(k in 1:(2*d)){
    count <- 1
    idx_jk[[k]] <- matrix(0, 2, no_eta*(no_eta+1)/2)
    if(k <= d){
      for(l in  k:d){
        for(m in  (d+1):(2*d)){
          if((l >=k) & ((m-d) >= l)){
            idx_jk[[k]][,count] <- c(l,m)
            count <- count + 1
          }
        }
        for(m in  (2*d+1):no_eta){
          if(((w[m-2*d] >=l) & (l > k) & (z[m-2*d]==k))|
             ((w[m-2*d] > l) & (l > k) & (z[m-2*d]==l))|
             (( k == l) & (l == z[m-2*d]))){
            idx_jk[[k]][,count] <- c(l,m)
            count <- count + 1
          }
        }
      }
      for(l in  (d+1):(2*d)){
        for(m in  l:(2*d)){
          if((l == m) & ((l-d) >= k)){
            idx_jk[[k]][,count] <- c(l,m)
            count <- count + 1
          }
        }
        for(m in (2*d+1):no_eta){
          if((((l-d) >= k ) | (k==z[m-2*d]) ) & ((l-d) == w[m-2*d])){
            idx_jk[[k]][,count] <- c(l,m)
            count <- count + 1
          }
        }
      }
      for(l in  (2*d+1):no_eta){
        for(m in l:no_eta){
          if(((k==z[l-2*d]) | (k==z[m-2*d]) ) & (w[l-2*d] == w[m-2*d])){
            idx_jk[[k]][,count] <- c(l,m)
            count <- count + 1
          }
        }
      }
    } else {
      for(l in  k:(2*d)){
        for(m in  l:(2*d)){
          if((l == k) & ( m == l)){
            idx_jk[[k]][,count] <- c(l,m)
            count <- count + 1
          }
        }
        for(m in  (2*d+1):no_eta){
          if(((l-d) == (k-d)) & ( w[m-2*d] == (l-d))){
            idx_jk[[k]][,count] <- c(l,m)
            count <- count + 1
          }
        }
      }
      for(l in  (2*d+1):no_eta){
        for(m in  l:no_eta){
          if((w[l-2*d] == (k-d)) & ( w[m-2*d] == (k-d))){
            idx_jk[[k]][,count] <- c(l,m)
            count <- count + 1
          }
        }
      }
    }
  }
  return(lapply(1:(2*d),function(x) matrix(idx_jk[[x]][1:2,which(idx_jk[[x]][1,]!=0)]-1,nrow=2)))
}


