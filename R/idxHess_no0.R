idxHess_no0 <- function(no_eta, z, w, param) {
  idx_jk_no0 <- list()

  # Dealing with sparsity of the mcd
  if( param == 1 ) {
    d <- -3/2 + sqrt(9/4 + 2 * no_eta)
    z1 <- z + 1
    w1 <- w + 1

    for(j in 1 : no_eta){
      idx_jk_no0[[j]] <- rep(0, no_eta)
      for(k in  j : no_eta){
        if ( j <= d & k <= d ) idx_jk_no0[[j]][k] <- k
        if ( j <= d & k > d & k <= 2 * d & k - j - d >= 0 ) idx_jk_no0[[j]][k] <- k
        if ( j <= d & k > 2 * d){ if ( w1[k - 2 * d] >= j | z1[k - 2 * d] == j ) idx_jk_no0[[j]][k] <- k }
        if ( j > d & j <= 2 * d & k > d & k <= 2 * d & j == k ) idx_jk_no0[[j]][k] <- k
        if ( j > d & j <= 2 * d & k > 2 * d){ if ( w1[k - 2 * d] == j - d ) idx_jk_no0[[j]][k] <- k }
        if ( j > 2 * d & k > 2 * d ){ if ( w1[j - 2 * d] == w1[k - 2 * d] ) idx_jk_no0[[j]][k] <- k }
      }
    }
  }

  # if I consider no sparsity in the mcd parametrisation
  #if(param == 1){
  #  for(j in 1 : no_eta){
  #    idx_jk_no0[[j]] <- rep(0, no_eta)
  #    for(k in  j : no_eta){
  #      idx_jk_no0[[j]][k] <- k
  #    }
  #  }
  #}

  if(param == 2){
    for(j in 1 : no_eta){
      idx_jk_no0[[j]] <- rep(0, no_eta)
      for(k in  j : no_eta){
        idx_jk_no0[[j]][k] <- k
      }
    }
  }
  return(lapply(1 : no_eta, function(x) which(idx_jk_no0[[x]] != 0) - 1))
}
