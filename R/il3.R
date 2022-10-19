il3 <- function(d){
  h <- choose(2:d,2)
  h2 <- choose((d-1):2,2)
  h3 <- rep(0,d-2)
  count <- 0
  for(k in count:(d-2)){
    h3[count+1] <- sum(count:(d-2))
    count <- count + 1
  }

  idx3_1 <- list()
  idx3_2 <- list()

  for(i in 1:length(h2)){
    idx3_1[[i]]<-rep(0,h2[i])
    idx3_2[[i]]<-rep(0,h2[i])
    count <- 1
    for(j in i: (d-2)){
      for(k in j: (d-2)){
        idx3_1[[i]][count] <- j
        idx3_2[[i]][count] <- k
        count <- count+1
      }
    }
  }

  idx3_3 <- list()
  idx3_4 <- list()
  idx3_5 <- list()
  idx3_6 <- list()

  for(i in 1:length(h3)){
    idx3_3[[i]]<-rep(0,h3[i])
    idx3_4[[i]]<-rep(0,h3[i])
    idx3_5[[i]]<-rep(0,h3[i])
    idx3_6[[i]]<-rep(0,h3[i])
    count <- 1
    for(j in 1:(d-1)){
      if(j > i){
        for(k in j:(d-1)){
          idx3_3[[i]][count] <- k-1
          idx3_4[[i]][count] <- i-1
          idx3_5[[i]][count] <- j-1
          idx3_6[[i]][count] <- j-1
          count <- count + 1
        }
      }
      if(j < i){
        for(k in i:(d-1)){
          idx3_3[[i]][count] <- k-1
          idx3_4[[i]][count] <- j-1
          idx3_5[[i]][count] <- i-1
          idx3_6[[i]][count] <- j-1
          count <- count + 1
        }
      }
    }
  }
 return(list(h = h, h2 = h2, h3=h3, idx3_1=idx3_1, idx3_2 = idx3_2, idx3_3 = idx3_3, idx3_4 = idx3_4, idx3_5 = idx3_5, idx3_6 = idx3_6 ))
}
