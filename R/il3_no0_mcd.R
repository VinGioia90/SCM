il3_no0_mcd <-  function(d, z, w){
  no_eta <- d+d*(d+1)/2
  z1 <- z+1
  w1 <- w+1

  idxj<-list()
  idxk<-list()
  for( i in 1:no_eta){
    idxj[[i]] <- rep(i, no_eta+1-i)
    idxk[[i]] <- i:no_eta
  }

  idxj<-unlist(idxj)
  idxk<-unlist(idxk)

  ixx<-list()

  idxq<-list()
  for(j in 1:(no_eta*(no_eta+1)/2)){
    ixx[[j]]<-idxq[[j]]<-rep(0, no_eta)
  }
  c1<-1

  csum<- cumsum(c(0,(no_eta-1):1))
  count <-1
  for(j in 1:no_eta){
    for(k in j:no_eta){
      for(m in 1:no_eta){
        if(m >= (d+1) & m <= 2*d & j <= d  & k <= d) {if( k>= j & (m-d) >= k){ixx[[c1]][m]<- count;  ixx[[m+csum[j]]][k] <- count; ixx[[m+csum[k]]][j]<- count ;  count <- count + 1;  idxq[[c1]][[m]]<- m; idxq[[m+csum[j]]][k] <- k; idxq[[m+csum[k]]][j]<- j}}
        if(m >= 2*d + 1 & j  < d & k <= d){if((((w1[m-2*d] >=k) & (k > j  ) & (z1[m-2*d]==j  ))|
                                               ((w1[m-2*d] > k) & (k > j  ) & (z1[m-2*d]==k))|
                                               (( j   == k) & (k == z1[m-2*d])))){ ixx[[c1]][m]<- count;  ixx[[m+csum[j]]][k] <- count;  ixx[[m+csum[k]]][j]<- count ; count <- count + 1;   idxq[[c1]][[m]]<- m; idxq[[m+csum[j]]][k] <- k; idxq[[m+csum[k]]][j]<- j}}
        if(m >= k & m <= 2*d & j  <= d  & k > d & k <= 2*d) { if(m==k & (k-d) >=j ){ ixx[[c1]][m]<- count; ixx[[m+csum[j]]][k] <- count; ixx[[m+csum[k]]][j]<- count ; count <- count + 1; idxq[[c1]][[m]]<- m; idxq[[m+csum[j]]][k] <- k; idxq[[m+csum[k]]][j]<- j}}
        if(m >= 2*d + 1 & j  <= d  & k > d & k<= 2*d) { if(( (k-d) >=j  | j ==z1[m-2*d] ) & (k-d) == w1[m-2*d]){ ixx[[c1]][m]<- count; ixx[[m+csum[j]]][k] <- count; ixx[[m+csum[k]]][j] <- count; count <- count + 1;  idxq[[c1]][[m]]<- m; idxq[[m+csum[j]]][k] <- k; idxq[[m+csum[k]]][j]<- j}}
        if(k > 2*d & m >= k){if(((j ==z1[k -2*d]) | (j ==z1[m-2*d]) ) & (w1[k -2*d] == w1[m-2*d])){ixx[[c1]][m]<- count;ixx[[m+csum[j]]][k] <- count; ixx[[m+csum[k]]][j] <- count; count <- count + 1;  idxq[[c1]][[m]]<- m; idxq[[m+csum[j]]][k] <- k; idxq[[m+csum[k]]][j]<- j}}
        if(j > d & j <=2*d & k>=j & k <= 2*d & m>=k & m<= 2*d){if(k==j & m==k){ixx[[c1]][m]<- count; ixx[[m+csum[j]]][k] <- count;  count <- count + 1;  idxq[[c1]][[m]]<- m; idxq[[m+csum[j]]][k] <- m}}
        if(j > d & j <=2*d & k>=j & k <= 2*d & m> 2*d){if((k-d)==(j-d)&  w1[m-2*d] == (k-d)){ixx[[c1]][m]<- count; ixx[[m+csum[j]]][k] <- count;   count <- count + 1;  idxq[[c1]][[m]]<- m; idxq[[m+csum[j]]][k] <- k}}
        if(j > d & j <=2*d & k> 2*d & m>= k){if((w1[k-2*d] == (j-d)) & ( w1[m-2*d] == (j-d))){ixx[[c1]][m]<- count; ixx[[m+csum[j]]][k] <- count; ixx[[m+csum[k]]][j] <- count; count <- count + 1;  idxq[[c1]][[m]]<- m; idxq[[m+csum[j]]][k] <- k; idxq[[m+csum[k]]][j]<- j}}
      }
      c1 <- c1+1
    }
  }
  aaa<-lapply(1 : (no_eta*(no_eta+1)/2), function(x) ixx[[x]][which(ixx[[x]] != 0)])
  bbb<-lapply(1 : (no_eta*(no_eta+1)/2), function(x) idxq[[x]][which(idxq[[x]] != 0)])
  il3_no0<-which(unlist(lapply(1 : (no_eta*(no_eta+1)/2), function(x) identical(aaa[[x]], numeric(0)))) == 0)
  return(list(idxl3 = lapply(il3_no0, function(x) aaa[[x]]),
              idxq = lapply(il3_no0, function(x) bbb[[x]]),
              idxj=idxj[il3_no0],
              idxk=idxk[il3_no0],
              csj = as.vector(cumsum(c(0,table(idxj[il3_no0]))))))
}


