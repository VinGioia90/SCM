aux_idx <- function(jj, idx_jk, no_eta){
  ll_jjx <- unlist(lapply(1 : length(jj), function(x) length(jj[[x]])))
  idx_jjx_l1 <- which(ll_jjx == 1) - 1
  idx_jjx_lg1 <- which(ll_jjx != 1) - 1

  b1_eta <- b1 <- b2 <- b3 <- list()

  for(j in 1 : length(idx_jk)){
    if ( (j - 1) %in% idx_jjx_lg1 ) {
      b1_eta[[j]] <- idx_jk[[j]]
      b1[[j]] <- idx_jk[[j]][idx_jk[[j]] %in% idx_jjx_lg1]
      b2[[j]] <- idx_jk[[j]][idx_jk[[j]] %in% idx_jjx_l1]
      b3[[j]] <- numeric(0)
    }  else {
      b1_eta[[j]] <- idx_jk[[j]][!(idx_jk[[j]] %in% idx_jjx_l1)]
      b1[[j]] <- numeric(0)
      b2[[j]] <- idx_jk[[j]][idx_jk[[j]] %in% idx_jjx_lg1]
      b3[[j]] <- idx_jk[[j]][idx_jk[[j]] %in% idx_jjx_l1]
    }
  }

  idx_b1_eta <- which(unlist(lapply(1 : length(idx_jk), function(x) identical(b1_eta[[x]], numeric(0)))) == 0) - 1
  idx_b3 <- idx_jjx_l1

  llls <- length(unlist(b3))
  idx_l2_full <- rep(0, length(unlist(b1)))
  idx_l2_partial <- rep(0, length(unlist(b2)))

  count_f1 <- count_f2 <- count_p1 <- count_p2 <- 1

  for(i in 1 : no_eta){
    for(j in idx_jk[[i]]){
      if ( !( j %in% idx_jjx_l1 ) & ((i - 1) %in% idx_jjx_lg1) ) {
        idx_l2_full[count_f2] <- count_f1 - 1
        count_f1 <- count_f1 + 1
        count_f2 <- count_f2 + 1
      }
      if ( ((j %in% idx_jjx_l1) & ((i - 1) %in% idx_jjx_lg1)) | (((i - 1) %in% idx_jjx_l1) & (j %in% idx_jjx_lg1)) ) {
        count_f1 <- count_f1 + 1
        idx_l2_partial[count_p2] <- count_p1 - 1
        count_p1 <- count_p1 + 1
        count_p2 <- count_p2 + 1
      }
      if ( (j %in% idx_jjx_lg1) & ((i - 1) %in% idx_jjx_lg1) ) count_p1 <- count_p1 + 1
    }
  }
  return(list(b1_eta = b1_eta, b1 = b1, b2 = b2, b3 = b3, idx_b1_eta = idx_b1_eta, idx_b3 = idx_b3,
              l2_el = idx_l2_full, l2_el2 = idx_l2_partial, llls = llls))
}
