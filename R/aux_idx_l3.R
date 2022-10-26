aux_idx_l3 <- function(d,z,w,G){
  # Such function allows to identify which elements are different from zero.
  # The idea is that, given k in l_klm the results are stored in a q x q matrix; with these indices
  # we are able in the final for loop of the l3_mcd_eta function  to save the results of the elements different from
  # zero in each of such matrices (k=1,..., q)
  no_eta <- d + d * (d + 1)/2
  idx_jk <- list()

  z <- z + 1
  w <- w + 1

  # Take a look to the formulas: we are saving the positions (l,m) of the k-th matrix where the indicator functions are equal to 1
  # Despite there is a more efficient way to extract such indices (by considering | and & operators), the computational cost of
  # computing them is low w.r.t. to rest of the procedure (such indices are calculated only one time when they are called to deriv > 1
  # and they are saved in memory); however, it should be possible to improve it
  for ( k in 1 : (2 * d) ) {
    count <- 1
    idx_jk[[k]] <- matrix(0, 2, no_eta * (no_eta + 1)/2)
    if ( k <= d ) {
      for ( l in  k : d ) {
        # Block (1,1,2)
        for ( m in  (d + 1) : (2 * d) ) {
          if ( (l >= k) & ((m - d) >= l) ) {
            idx_jk[[k]][, count] <- c(l, m)
            count <- count + 1
          }
        }
        # Block (1,1,3)
        for ( m in  (2 * d + 1) : no_eta ) {
          mp <- m - 2 * d
          if ( ((w[mp] >= l) & (l > k) & (z[mp] == k)) |
               ((w[mp] > l) & (l > k) & (z[mp] == l))|
               (( k == l) & (l == z[mp])) ) {
            idx_jk[[k]][, count] <- c(l, m)
            count <- count + 1
          }
        }
      }

      for ( l in  (d + 1) : (2 * d) ) {
        # Block (1,2,2)
        for ( m in  l : (2 * d) ) {
          if ( (l == m) & ((l - d) >= k) ) {
            idx_jk[[k]][, count] <- c(l, m)
            count <- count + 1
          }
        }
        # Block (1,2,3)
        for (m in (2 * d + 1) : no_eta ) {
          mp <- m - 2 * d
          if ( (((l - d) >= k) | (k == z[mp])) & ((l - d) == w[mp]) ) {
            idx_jk[[k]][, count] <- c(l, m)
            count <- count + 1
          }
        }
      }
      for (l in  (2 * d + 1) : no_eta ) {
        lp <- l - 2 * d
        # Block (1,3,3)
        for ( m in l : no_eta ) {
          mp <- m - 2 * d
          if ( ((k == z[lp]) | (k == z[mp]) ) & (w[lp] == w[mp]) ) {
            idx_jk[[k]][, count] <- c(l, m)
            count <- count + 1
          }
        }
      }
    } else {
      for ( l in  k : (2 * d) ) {
        # Block (2,2,2)
        for ( m in  l : (2 * d) ) {
          if ( (l == k) & (m == l) ) {
            idx_jk[[k]][, count] <- c(l, m)
            count <- count + 1
          }
        }
        # Block (2,2,3)
        for ( m in  (2 * d + 1) : no_eta ) {
          if ( ((l - d) == (k - d)) & (w[m - 2 * d] == (l - d)) ){ # l-d == k-d should be simply l == k
            idx_jk[[k]][, count] <- c(l, m)
            count <- count + 1
          }
        }
      }
      for ( l in  (2 * d + 1) : no_eta ) {
        # Block (2,3,3)
        for ( m in  l : no_eta ) {
          if ( (w[l - 2 * d] == (k - d)) & (w[m - 2 * d] == (k - d)) ) {
            idx_jk[[k]][, count] <- c(l, m)
            count <- count + 1
          }
        }
      }
    }
  }
  return(lapply(1 : (2 * d), function(x) matrix(idx_jk[[x]][1 : 2, which(idx_jk[[x]][1,] != 0)] - 1, nrow = 2)))
}


