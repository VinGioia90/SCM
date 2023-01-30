gamlss.gH_scm <- function (X, jj, eta, y, w, z, t, Gm, # l3 = 0,
                           l1, l1_l, l2, l2_v, l2_l, l2_v_l, idx_b, idx_aux, param = 1, #by default param=1 for the mcd parameterisation
                           l3, l3_l, idx_l3, idx_neq0, idx_jkq,
                           #d1eta, d1eta_l, V, V_l, d1H, #Such elements were initialised in the cpp version
                           l4 = 0, i4 = 0, d1b = 0, d2b = 0, deriv = 0, fh = NULL, D = NULL){
  K <- length(jj) # no_eta
  sparse <- discrete <- FALSE
  p <- ncol(X)
  n <- nrow(X)

  trHid2H <- d2H <- NULL
  if( deriv < 1 ) d1H <- NULL


  lpi <- lapply(jj, function(x) x - 1) # trick to allow the use of the c++functions
  #####################
  # Score w.r.t. beta #
  #####################
  lb <- rep(0, p)
  internal()$d1_beta(X, eta, y, lpi, K, lb, l1, l1_l, idx_b, z, w, Gm, param) #first derivatives w.r.t beta

  #######################
  # Hessian w.r.t. beta #
  #######################
  lbb <- if ( sparse )
    Matrix(0, p, p)
  else matrix(0, p, p) # dovremmo inizializzarla?????

  internal()$d2_beta(X, eta, y, lpi, K, lbb, l2, l2_v, l2_l, l2_v_l, idx_b, z, w, Gm, t, #second derivatives w.r.t beta
                     idx_aux$b1_eta, idx_aux$b1, idx_aux$b2, idx_aux$b3, idx_aux$idx_b1_eta, idx_aux$idx_b3,
                     idx_aux$l2_el, idx_aux$l2_el2 , param)
  #######################################
  # Derivative of the Hessian w.r.t eta #
  #######################################
  if ( deriv > 0 & param == 1) {
    if (discrete) stop("er... no discrete methods for higher derivatives")
      ## the first derivative of the Hessian, using d1b
      ## the first derivates of the coefficients wrt the sps
      m <- ncol(d1b) ## number of smoothing parameters
      ## stack the derivatives of the various linear predictors on top
      ## of each other...
      d1eta <- matrix(0,n*K,m)
      ind <- 1:n
      for (i in 1:K) {
        d1eta[ind,] <- if (discrete) Xbd(X$Xd,d1b,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[i]]) else
          X[,jj[[i]],drop=FALSE]%*%d1b[jj[[i]],]
        ind <- ind + n
      }

      #new R version: it works in mgcv version 42

          internal()$d3_mcd_eta(eta, y, l3, z, w,  Gm,  t,  idx_l3$h, idx_l3$h2, idx_l3$h3, idx_l3$idx3_1, idx_l3$idx3_2, idx_l3$idx3_3, idx_l3$idx3_4, idx_l3$idx3_5, idx_l3$idx3_6, idx_neq0);

         #Inefficient
         #d1H <- rep(0,m)
         #    for (l in 1:m) { ## sp loop
         #     #v <- rep(0,n)
         #      for(t in 1:length(idx_jkq$idxl3)){
         #        v <- rep(0,n)
         #
         #        for (q in 1:length(idx_jkq$idxl3[[t]])) {
         #          v <- v + l3[,idx_jkq$idxl3[[t]][q]] * d1eta[((idx_jkq$idxq[[t]][q]-1)*n + 1):(idx_jkq$idxq[[t]][q]*n),l]
         #        }
         #     a <- rowSums((X[,jj[[idx_jkq$idxj[t]]],drop=FALSE] %*% fh[jj[[idx_jkq$idxj[t]]],jj[[idx_jkq$idxk[t]]]]) * X[,jj[[idx_jkq$idxk[t]]],drop=FALSE])
         #     mult <- if (idx_jkq$idxk[t]==idx_jkq$idxj[t]) 1 else 2
         #     d1H[l] <- d1H[l] + mult * (sum(a*v)) ## accumulate tr(Hp^{-1}dH/drho_l)
         #    }
         #}
         #Efficient
          d1H <- rep(0,m)
           for(t in 1:length(idx_jkq$idxl3)){
               a <- rowSums((X[,jj[[idx_jkq$idxj[t]]],drop=FALSE] %*% fh[jj[[idx_jkq$idxj[t]]],jj[[idx_jkq$idxk[t]]]]) * X[,jj[[idx_jkq$idxk[t]]],drop=FALSE])
               mult <- if (idx_jkq$idxk[t]==idx_jkq$idxj[t]) 1 else 2
               for (l in 1:m) { ## sp loop
                 v <- rep(0,n)
                 for (q in 1:length(idx_jkq$idxl3[[t]])) {
                   v <- v + l3[,idx_jkq$idxl3[[t]][q]] * d1eta[((idx_jkq$idxq[[t]][q]-1)*n + 1):(idx_jkq$idxq[[t]][q]*n),l]
                 }
                 d1H[l] <- d1H[l] + mult * (sum(a*v)) ## accumulate tr(Hp^{-1}dH/drho_l)
               }
             }


      #old R version: it works with mgcv version 41

      # d3_mcd_eta(eta, y, l3, z, w,  Gm,  t, idx_l3$h, idx_l3$h2, idx_l3$h3, idx_l3$idx3_1, idx_l3$idx3_2, idx_l3$idx3_3, idx_l3$idx3_4, idx_l3$idx3_5, idx_l3$idx3_6, idx_neq0)
        # d1H <- list()
        # for (l in 1:m) {
        #   d1H[[l]] <- matrix(0,p,p)
        #
        #   for(j in 1:length(idx_jkq$idxl3)){
        #     v <- rep(0,n)
        #
        #     for (q in 1:length(idx_jkq$idxl3[[j]])) {
        #       v <- v + l3[,idx_jkq$idxl3[[j]][q]] * d1eta[((idx_jkq$idxq[[j]][q]-1)*n + 1):(idx_jkq$idxq[[j]][q]*n),l]
        #     }
        #     A <- crossprod(X[,jj[[idx_jkq$idxj[j]]],drop=FALSE],v*X[,jj[[idx_jkq$idxk[j]]],drop=FALSE])
        #     d1H[[l]][jj[[idx_jkq$idxj[j]]],jj[[idx_jkq$idxk[j]]]] <- d1H[[l]][jj[[idx_jkq$idxj[j]]],jj[[idx_jkq$idxk[j]]]] + A
        #     if (idx_jkq$idxk[j]>idx_jkq$idxj[j]) d1H[[l]][jj[[idx_jkq$idxk[j]]],jj[[idx_jkq$idxj[j]]]] <- d1H[[l]][jj[[idx_jkq$idxk[j]]],jj[[idx_jkq$idxj[j]]]] + t(A)
        #   }
        # }
  }
  if ( deriv > 0 & param == 2) {
    warning("bfgs is not available for the logm parametrisation: the fitting is done by using efs optimizer")
    d1H <- NULL
  }
  list(lb = lb, lbb = lbb, d1H = d1H, d2H = d2H, trHid2H = trHid2H)
}
