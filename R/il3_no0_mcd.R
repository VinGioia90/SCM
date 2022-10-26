il3_no0_mcd <-  function(d, z, w){
  # Such function allows to access to the indices needed to compute the derivative of the Hessian w.r.t. the smoothing parameters
  # (taking into account the sparsity). It returns:
  # -) idxl3: column indices of the third derivatives matrix useful to build the V matrix
  # -) idxq: column indices of the derivatives deta/drho useful to build the V matrix
  # -) idxj: indices related to the passage to the linear predictor to the Xj elements
  # -) idxk: indices related to the passage to the linear predictor to the Xk elements

  # After Matteo's improvement, consider or not the intercept and/or block cases.

  no_eta <- d + d * (d + 1)/2
  z1 <- z + 1
  w1 <- w + 1

  idxj <- idxk <- list()
  for ( i in 1 : no_eta ) {
    idxj[[i]] <- rep(i, no_eta + 1 - i)
    idxk[[i]] <- i : no_eta
  }

  idxj <- unlist(idxj)
  idxk <- unlist(idxk)

  ixx <- idxq <- list()
  K <- no_eta * (no_eta + 1)/2
  for ( j in 1 : K )  ixx[[j]] <- idxq[[j]] <- rep(0, no_eta)

  c1 <- 1
  csum <- cumsum(c(0, (no_eta - 1) : 1))
  count <- 1
  for ( j in 1 : no_eta ) {
    for ( k in j : no_eta ) {
      for ( m in 1 : no_eta ) { #Take a look to the formulas: the if statements are related to the indicator function
        if ( m >= (d + 1) & m <= 2 * d & j <= d  & k <= d) { #Block (1, 1, 2)
          if ( k>= j & (m-d) >= k ) {
            ixx[[c1]][m] <- ixx[[m+csum[j]]][k] <- ixx[[m+csum[k]]][j]<- count
            count <- count + 1
            idxq[[c1]][[m]] <- m
            idxq[[m + csum[j]]][k] <- k
            idxq[[m + csum[k]]][j] <- j
          }
        }
        if ( m >= (2 * d + 1) & j  < d & k <= d ) { #Block (1, 1, 3)
          if ( (((w1[m - 2 * d] >= k) & k > j & (z1[m - 2 * d] == j)) |
               ((w1[m - 2 * d] > k) & (k > j) & (z1[m - 2 * d] == k)) |
               ((j == k) & (k == z1[m - 2 * d])))){
            ixx[[c1]][m] <- ixx[[m+csum[j]]][k] <- ixx[[m+csum[k]]][j] <- count
            count <- count + 1
            idxq[[c1]][[m]] <- m
            idxq[[m + csum[j]]][k] <- k
            idxq[[m + csum[k]]][j] <- j
          }
        }
        if ( m >= k & m <= 2*d & j  <= d  & k > d & k <= 2*d) { if(m==k & (k-d) >=j ){ ixx[[c1]][m]<- count; ixx[[m+csum[j]]][k] <- count; ixx[[m+csum[k]]][j]<- count ; count <- count + 1; idxq[[c1]][[m]]<- m; idxq[[m+csum[j]]][k] <- k; idxq[[m+csum[k]]][j]<- j}}
        if ( m >= 2*d + 1 & j  <= d  & k > d & k<= 2*d) { if(( (k-d) >=j  | j ==z1[m-2*d] ) & (k-d) == w1[m-2*d]){ ixx[[c1]][m]<- count; ixx[[m+csum[j]]][k] <- count; ixx[[m+csum[k]]][j] <- count; count <- count + 1;  idxq[[c1]][[m]]<- m; idxq[[m+csum[j]]][k] <- k; idxq[[m+csum[k]]][j]<- j}}
        if ( k > 2*d & m >= k){if(((j ==z1[k -2*d]) | (j ==z1[m-2*d]) ) & (w1[k -2*d] == w1[m-2*d])){ixx[[c1]][m]<- count;ixx[[m+csum[j]]][k] <- count; ixx[[m+csum[k]]][j] <- count; count <- count + 1;  idxq[[c1]][[m]]<- m; idxq[[m+csum[j]]][k] <- k; idxq[[m+csum[k]]][j]<- j}}
        if ( j > d & j <=2*d & k>=j & k <= 2*d & m>=k & m<= 2*d){if(k==j & m==k){ixx[[c1]][m]<- count; ixx[[m+csum[j]]][k] <- count;  count <- count + 1;  idxq[[c1]][[m]]<- m; idxq[[m+csum[j]]][k] <- m}}
        if ( j > d & j <=2*d & k>=j & k <= 2*d & m> 2*d){if((k-d)==(j-d)&  w1[m-2*d] == (k-d)){ixx[[c1]][m]<- count; ixx[[m+csum[j]]][k] <- count;   count <- count + 1;  idxq[[c1]][[m]]<- m; idxq[[m+csum[j]]][k] <- k}}
        if ( j > d & j <=2*d & k> 2*d & m>= k){if((w1[k-2*d] == (j-d)) & ( w1[m-2*d] == (j-d))){ixx[[c1]][m]<- count; ixx[[m+csum[j]]][k] <- count; ixx[[m+csum[k]]][j] <- count; count <- count + 1;  idxq[[c1]][[m]]<- m; idxq[[m+csum[j]]][k] <- k; idxq[[m+csum[k]]][j]<- j}}
      }
      c1 <- c1+1
    }
  }

  il3_no01 <- lapply(1 : K, function(x) ixx[[x]][which(ixx[[x]] != 0)])

  il3_no02 <- lapply(1 : K, function(x) idxq[[x]][which(idxq[[x]] != 0)])

  il3_no0 <- which(unlist(lapply(1 : K, function(x) identical(il3_no01[[x]], numeric(0)))) == 0)

  return(list(idxl3 = lapply(il3_no0, function(x) il3_no01[[x]]),
              idxq = lapply(il3_no0, function(x) il3_no02[[x]]),
              idxj = idxj[il3_no0],
              idxk = idxk[il3_no0]
              ))
}


