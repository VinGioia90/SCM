#' Multivariate Gaussian Additive model (family) with covariance modelling
#'
#' @description  Such family contains all the needed functions to fit and post-process the model
#' @param d dimension of the outcome
#' @param nb number of observations' blocks
#' @param param type of parameterisation: MCD or logM
#'
#' @return family needed quantities
#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom mvnfast rmvn
mvn_scm <- function(d = 2, nb = 1, param = NULL){ # manage internally the blocks division !!!

  if ( d < 2 ) stop("mvn_scm requires to or more dimensional data")
  if ( nb < 1 ) stop("the number of observations' blocks must be greater than 0")

  no_eta <- d + d * (d + 1)/2 # dimension of the linear predictor

  # Set "mcd parametrisation" to 1 and logm to 2
  if ( is.null(param) ) param <- "mcd"
  param2 <- param # copy to print the name of the parametrisation used
  if ( param == "mcd" )     param <- 1
  else if ( param == "logm" ) param <- 2
  else stop( "Wrong parametrisation chosen" )

  # Set the available derivatives according to the parametrisation chosen
  # It allows to choose the type of smoothing parameter optimization
  if( param == 1 ) a.der <- 1
  if( param == 2 ) a.der <- 0

  stats <- list()
  for ( j in 1 : no_eta ) stats[[j]] <- make.link("identity")
  validmu <- function(mu) all(is.finite(mu))

  # z, w, G, t: auxiliary indices used in the derivative's computation (all of them used in the mcd and only z, w in the logm)
  z <- w <- t <- rep(0, (d * (d - 1)/2))
  Gm <- matrix(0, d - 1, d - 1)
  mode(Gm) <- mode(z) <- mode(w) <- mode(t) <- "integer"
  internal()$idx_zwGt(d, z, w, Gm, t)

  # indices of elements different from zero in the hessian matrix (logm is not sparse)
  idx_jk <- internal()$idxHess_no0(no_eta, z, w, param)

  #number of (non-redundant) elements different from zero in the hessian matrix
  if ( param == 1 ) { # See the computational paper
    nHel <- d * (d^2 + 15 * d + 2)/6
    if ( d > 2 ) nHel <- nHel + d * (d - 1) * (d - 2)/3
  }
  if ( param == 2 ) nHel <- no_eta * (no_eta + 1)/2 #logm: no sparsity

  # Quantities defined in the environment:
  # l1 and l1_l: matrix of the 1st derivatives
  getL1 <- function() get(".l1")   # only the first nb - 1 blocks of observations
  putL1 <- function(.l1) assign(".l1", .l1, envir = environment(sys.function()))
  getL1_l <- function() get(".l1_l")   # Last block of observations
  putL1_l <- function(.l1_l) assign(".l1_l", .l1_l, envir = environment(sys.function()))

  # l2 and l2_l: matrices of the 2nd derivatives excluding the intercepts' blocks
  getL2 <- function() get(".l2") # only the first nb - 1 blocks of observations
  putL2 <- function(.l2) assign(".l2", .l2, envir = environment(sys.function()))
  getL2_l <- function() get(".l2_l")   # Last block of observations
  putL2_l <- function(.l2_l) assign(".l2_l", .l2_l, envir = environment(sys.function()))

  # l2_v and l2_v_l: vectors of the (cumulated) second derivatives considering only the  intercepts' blocks
  getL2_v <- function() get(".l2_v")  # First nb - 1 blocks of observations
  putL2_v <- function(.l2_v) assign(".l2_v", .l2_v, envir = environment(sys.function()))
  getL2_v_l <- function() get(".l2_v_l")   # Last block of observations
  putL2_v_l <- function(.l2_v_l) assign(".l2_v_l", .l2_v_l, envir = environment(sys.function()))

  # l3 and l3_l: is the matrix of the third derivatives w.r.t. eta
  getL3 <- function() get(".l3")  # First nb - 1 blocks of observations
  putL3 <- function(.l3) assign(".l3", .l3, envir = environment(sys.function()))
  getL3_l <- function() get(".l3_l")  # Last block of observations
  putL3_l <- function(.l3_l) assign(".l3_l", .l3_l, envir = environment(sys.function()))

  # d1H is the matrix of the derivatives of the Hessian w.r.t smoothing parameters
  getd1H <- function() get(".d1H")
  putd1H <- function(.d1H) assign(".d1H", .d1H, envir = environment(sys.function()))

  # d1eta and d1_eta_l is the matrix of the derivatives of eta w.r.t  smoothing parameters
  getd1eta <- function() get(".d1eta") # First nb - 1 blocks of observations
  putd1eta <- function(.d1eta) assign(".d1eta", .d1eta, envir = environment(sys.function()))
  getd1eta_l <- function() get(".d1eta_l")  # Last block of observations
  putd1eta_l <- function(.d1eta_l) assign(".d1eta_l", .d1eta_l, envir = environment(sys.function()))

  # V matrix: see Wood et al. 2016
  getV <- function() get(".V") # First nb - 1 blocks of observations
  putV <- function(.V) assign(".V", .V, envir = environment(sys.function()))
  getV_l <- function() get(".V_l")  # Last block of observations
  putV_l <- function(.V_l) assign(".V_l", .V_l, envir = environment(sys.function()))

  # Indices for observations' blocks
  getidx_b <- function() get(".idx_b")
  putidx_b <- function(.idx_b) assign(".idx_b", .idx_b, envir = environment(sys.function()))

  # List of auxiliary indices for hessian blocks  (intercept/partial/full taking into account the sparsity if it exists)
  getidx_aux <- function() get(".idx_aux")
  putidx_aux <- function(.idx_aux) assign(".idx_aux", .idx_aux, envir = environment(sys.function()))

  # indices for the third derivatives building
  getidxl3 <- function() get(".idxl3")
  putidxl3 <- function(.idxl3) assign(".idxl3", .idxl3, envir = environment(sys.function()))

  getidxl3_no0 <- function() get(".idxl3_n00")
  putidxl3_no0 <- function(.idxl3_no0) assign(".idxl3_no0", .idxl3_no0, envir = environment(sys.function()))

  getidxl3_jkq <- function() get(".idxl3_jkq")
  putidxl3_jkq <- function(.idxl3_jkq) assign(".idxl3_jkq", .idxl3_jkq, envir = environment(sys.function()))


  # Correlation flag: allowing to obtain the correlations (TRUE) or the covariances (FALSE)
  assign(".cflag", TRUE, envir = environment())
  getcflag <- function() get(".cflag")
  putcflag <- function(.cflag) assign(".cflag", .cflag, envir = environment(sys.function()))

  # Dimension of the outcome
  assign(".d", d, envir = environment())
  getd <- function() get(".d")
  putd <- function(.d) assign(".d", .d, envir = environment(sys.function()))

  # Dimension of the linear predictor vector
  assign(".no_eta", no_eta, envir = environment())
  getno_eta <- function() get(".no_eta")
  putno_eta <- function(.no_eta) assign(".no_eta", .no_eta, envir = environment(sys.function()))

  # Type of parametrisation
  assign(".param", param, envir = environment())
  getparam <- function() get(".param")
  putparam <- function(.param) assign(".param", .param, envir = environment(sys.function()))


  # Initialize: we initialize the mean vector parameters and only the intercepts related to the decomposition or transformation involved
  initialize <- expression({
    my_init_fun <- function(y, nobs, E, x, family, offset){
      d <- family$getd()
      no_eta <- family$getno_eta()
      param <- family$getparam()

      use.unscaled <- if ( !is.null(attr(E,"use.unscaled")) ) TRUE else FALSE
      jj <- attr(x, "lpi")
      resid <-  matrix(0, nrow(y), d)

      if ( is.list(x) ) { # discrete case: incomplete !!!! (Initialization of the intercepts is missing in the discrete case)
        start <- rep(0, max(unlist(jj)))
        for ( k in 1 : d ) {
          yt1 <- y[, k]
          e1 <- E[, jj[[k]], drop = FALSE] ## square root of total penalty
          R <- suppressWarnings(chol(XWXd(x$Xd, w = rep(1, length(yt1)),
                                          k = x$kd, ks = x$ks, ts = x$ts,
                                          dt = x$dt, v = x$v, qc = x$qc, nthreads = 1,
                                          drop = x$drop, lt = x$lpid[[k]]) +
                                       crossprod(E[, jj[[k]]]), pivot = TRUE))
          Xty <- XWyd(x$Xd, rep(1, nrow(y)), yt1, x$kd,
                      x$ks, x$ts, x$dt, x$v, x$qc, x$drop, lt = x$lpid[[k]])
          piv <- attr(R, "pivot")
          rrank <- attr(R, "rank")
          startji <- rep(0, ncol(R))
          if ( rrank < ncol(R) ) {
            R <- R[1 : rrank, 1 : rrank]
            piv <- piv[1 : rrank]
          }

          startMu <- numeric(ncol(R))
          startMu[piv] <- backsolve(R, forwardsolve(t(R), Xty[piv]))
          startMu[!is.finite(startji)] <- 0
          start[jj[[k]]] <- startMu
          eta1 <- Xbd(x$Xd, start, k = x$kd, ks = x$ks,
                      ts = x$ts, dt = x$dt, v = x$v, qc = x$qc, drop = x$drop,
                      lt = x$lpid[[k]])
          resid[, k] <- y[, k] - eta1
        }

        Svcov <- cov(resid) #empirical covariance matrix of the residual vector
        if ( param == 1 ) dec <- mcd(Svcov) # mcd
        if ( param == 2 ) dec <- logm(Svcov) # logm
        Theta_el <- c(diag(dec), dec[upper.tri(dec, diag = FALSE)])
      } else { #regular case
        start <- rep(0, ncol(x))
        for ( k in 1 : d ) { # mean vector initialization
          yt1 <- y[, k]
          x1 <- x[, jj[[k]], drop = FALSE]
          e1 <- E[, jj[[k]], drop = FALSE] ## square root of total penalty

          if ( use.unscaled ) {
            x1 <- rbind(x1, e1)
            startMu <- qr.coef(qr(x1), c(yt1, rep(0, nrow(E))))
            startMu[!is.finite(startMu)] <- 0
          } else {
            startMu <- pen.reg(x1, e1, yt1)
          }
          if ( !is.matrix(x[, jj[[k]]]) ) { # To take into account the case when  x[, jj[[k]]] is a vector (--> intercepts)
            resid[, k] <- y[, k] - x[, jj[[k]]] * startMu
          } else {
            resid[, k] <- y[,k] - x[, jj[[k]]] %*% startMu
          }
          start[jj[[k]]] <- startMu
        }

        Svcov <- cov(resid) #empirical covariance matrix of the residual vector
        if ( param == 1 ) dec <- mcd(Svcov) #mcd
        if ( param == 2 ) dec <- logm(Svcov) #logm
        Theta_el <- c(diag(dec), dec[upper.tri(dec, diag = FALSE)])

        for ( k in (d + 1) : no_eta ) { #Initialization of the intercepts of the decomposition/transformation
          x1 <- x[, jj[[k]], drop = FALSE]
          startji <- qr.coef(qr(x1), c(rep(Theta_el[k - d], nrow(x1))))
          startji[!is.finite(startji)] <- 0
          start[jj[[k]]] <- startji
        }
      }
      return(start)
    }
    if ( is.null(start) ) {
      start <- my_init_fun(y = y, nobs = nobs, E = E, x = x, family = family, offset = offset)
    }
  }) ## initialize

  ## Residuals
  residuals <- function(object, type = c("response", "deviance")) { #by defualt deviance residuals
    type <- match.arg(type)

    if ( type == "deviance" ) {
      n <- dim(object$fitted.values)[1]
      res <- matrix(0, n, d)
      if ( param == 1 ) internal()$res_dev_mcd(object$fitted.values, object$y, res) #Deviance residuals for mcd
      if ( param == 2 ) internal()$res_dev_logm(object$fitted.values, object$y, res) #Deviance residuals for logm
    } else {
      res <- object$y - object$fitted.values[, 1 : d]
    }
    res
  } ## residuals

  # ll function
  ll <- function(y, X, coef, wt, family,  offset = NULL,
                 deriv = 0, d1b = NULL, d2b = NULL, Hp = NULL,
                 rank = 0, fh = NULL, D = NULL) {
    discrete <- is.list(X)
    jj <- attr(X, "lpi") ## extract linear predictor index
    p <- ncol(X)
    n <- nrow(y)
    eta <- matrix(0, n, no_eta) #linear predictor matrix

    # number of observations in each block (nb - 1 blocks with the same observations; the last could differ)
    nlast <- n %% nb
    nset <- n %/% nb

    # Memory initialization of several quantities:
    l1 <- try(getL1(), TRUE) #First derivatives matrix (first nb - 1 blocks of observations)
    if ( "try-error" %in% class(l1) ) {
      l1 <- matrix(0, nset, no_eta)
      putL1(l1)
    }

    idx_aux <- try(getidx_aux(), TRUE) #list of auxiliary indices for hessian blocks (intercept/partial/full taking into account the sparsity if it exists)
    if ( "try-error" %in% class(idx_aux) ) {
      idx_aux <- internal()$aux_idx(jj, idx_jk, no_eta)
      putidx_aux(idx_aux)
    }

    if(nb > 1){
      l2 <- try(getL2(), TRUE)  #Second derivatives matrix (first nb - 1 blocks of observations and Hessian blocks not involving the intercepts)
      if ( "try-error" %in% class(l2) ) {
        l2 <- matrix(0, nset, nHel - idx_aux$llls)
        putL2(l2)
      }
    } else { # Such trick allows to pass the matrix l2 in the case nb = 1, avoiding cpp issues
      l2 <- try(getL2(), TRUE)  #Second derivatives matrix (first nb - 1 blocks of observations and Hessian blocks not involving the intercepts)
      if ( "try-error" %in% class(l2) ) {
        l2 <- matrix(0, 1, 1)
        putL2(l2)
      }
    }

    l2_v <- try(getL2_v(), TRUE) #Second derivatives vector (first nb - 1 blocks of observations and Hessian blocks involving the intercepts)
    if ( "try-error" %in% class(l2_v) ) {
      l2_v <- rep(0, idx_aux$llls)
      putL2_v(l2_v)
    }

    l2_v_l <- try(getL2_v_l(), TRUE)  #Second derivatives vector (last blocks of observations and Hessian blocks not involving the intercepts)
    if ( "try-error" %in% class(l2_v_l) ) {
      l2_v_l <- rep(0, idx_aux$llls)
      putL2_v_l(l2_v_l)
    }

    # Create the last block quantities
    if ( nlast == 0 ) {
      nobs_b <- nset
      idx_b_seq <- cumsum(rep(nset, nb - 1)) - 1
    } else {
      nobs_b <- nlast
     idx_b_seq <- cumsum(rep(nset, nb)) - 1
    }

    l1_l <- try(getL1_l(), TRUE)  #First derivatives matrix (Last block of observations)
    if ( "try-error" %in% class(l1_l) ) {
      l1_l <- matrix(0, nobs_b,  no_eta)
      putL1_l(l1_l)
    }

    l2_l <- try(getL2_l(), TRUE) #Second derivatives matrix (Last block of observations)
    if("try-error" %in% class(l2_l)){
      l2_l <- matrix(0, nobs_b, nHel - idx_aux$llls)
      putL2_l(l2_l)
    }

    idx_b <- try(getidx_b(), TRUE)  # Indices of observations' blocks
    if("try-error" %in% class(idx_b)){
      idx_b <- c(-1, idx_b_seq, n - 1)
      putidx_b(idx_b)
    }

    d1H <- NULL
    if ( deriv > 1 ) { # only for the mcd parametrisation
      idxl3 <- try(getidxl3(), TRUE) # indices for the third derivatives building
      if("try-error" %in% class(idxl3)){
        idxl3 <- internal()$il3(d)
        putidxl3(idxl3)
      }

      idxl3_no0 <- try(getidxl3_no0(), TRUE) # indices for the third derivatives building
      if("try-error" %in% class(idxl3_no0)){
        idxl3_no0 <- internal()$aux_idx_l3(d,z,w,Gm)
        putidxl3_no0(idxl3_no0)
      }

      idxl3_jkq <- try(getidxl3_jkq(), TRUE) # indices for the third derivatives building
      if("try-error" %in% class(idxl3_jkq)){
        idxl3_jkq <- internal()$il3_no0_mcd(d, z, w)
        putidxl3_jkq(idxl3_jkq)
      }

      if(nb > 1){
        l3 <- try(getL3(), TRUE)  # Third derivatives matrix: First nb - 1 blocks of observations (the intercepts block case is not considered at the moment)
        if("try-error" %in% class(l3)){
          l3 <- matrix(0, nset, d * (4 * d^2 + 3 * d + 2)/3)
          putL3(l3)
        }
      } else { # Such trick allows to pass the matrix l3 in the case nb = 1, avoiding cpp issues
        l3 <- try(getL3(), TRUE)  # Third derivatives matrix: First nb - 1 blocks of observations
        if("try-error" %in% class(l3)){
          l3 <- matrix(0, 1, 1)
          putL3(l3)
        }
      }

      d1H <- try(getd1H(), TRUE) # derivative of hessian w.r.t. smoothing parameters
      if("try-error" %in% class(d1H)){
        d1H <- list()
        m <- ncol(d1b)
        for ( l in 1 : m ) d1H[[l]] <- matrix(0, p, p)
        putd1H(d1H)
      }

      d1eta <- try(getd1eta(), TRUE) #d1eta is the matrix (First nb - 1 observations' blocks) of the derivatives of eta w.r.t  smoothing parameters
      if("try-error" %in% class(d1eta)){
        d1eta <-  matrix(0, nset, no_eta)
        putd1eta(d1eta)
      }

      V <- try(getV(), TRUE) # V matrix: see Wood et. al (2016)
      if("try-error" %in% class(V)){
        V <-  rep(0, nset)
        putV(V)
      }

      if ( nlast == 0 ) nobs_b <- nset
      else  nobs_b <- nlast

      l3_l <- try(getL3_l(), TRUE)  #Third derivatives matrix: Last observations' blocks
      if("try-error" %in% class(l3_l)){
        l3_l <- matrix(0, nobs_b, d * (4 * d^2 + 3 * d + 2)/3)
        putL3_l(l3_l)
      }

      d1eta_l <- try(getd1eta_l(), TRUE) #d1eta_l is the matrix (Last observations' blocks) of the derivatives of eta w.r.t  smoothing parameters
      if("try-error" %in% class(d1eta_l)){
        d1eta_l <-  matrix(0, nobs_b,no_eta)
        putd1eta_l(d1eta_l)
      }

      V_l <- try(getV_l(), TRUE) # V matrix: see Wood et. al (2016)
      if("try-error" %in% class(V_l)){
        V_l <-  rep(0, nobs_b)
        putV_l(V_l)
      }
    }

    # Building the linear predictor vector
    for ( k in 1 : no_eta )  eta[, k] <- if ( discrete ) Xbd(X$Xd, coef, k = X$kd, ks= X$ks,ts = X$ts, dt =X$dt,v = X$v, qc = X$qc, drop = X$drop, lt = X$lpid[[k]])
                                         else X[, jj[[k]], drop = FALSE] %*% coef[jj[[k]]]


    ## log-likelihood
    if ( param == 1 )  l <- internal()$ll_mcd(eta, y)  #mcd
    if ( param == 2 )  l <- internal()$ll_logm(eta, y) #logm
    l <- l - 0.5 * n * d * log(2 * pi)

    # To return first, second derivatives w.r.t. beta and the derivatives of the Hessian w.r.t. to the smoothing parameters (if bfgs)
    if ( deriv ) {
      ret <- internal()$gamlss.gH_scm(X, jj, eta, y, w, z, t, Gm,
                                      l1, l1_l, l2, l2_v, l2_l, l2_v_l,
                                      idx_b, idx_aux, param = param,
                                      l3, l3_l, idxl3, idxl3_no0, idxl3_jkq,
                                      d1eta, d1eta_l, V, V_l, d1H,
                                      d1b = d1b, deriv = deriv - 1, fh = fh, D = D)
    } else ret <- list()
    ret$l <- l
    ret
  } ## end ll


  predict <- function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
                      beta = NULL, off = NULL, Vb = NULL) {
    ## optional function to give predicted values - idea is that
    ## predict.gam(...,type="response") will use this, and that
    ## either eta will be provided, or {X, beta, off, Vb}. family$data
    ## contains any family specific extra information.
    ## if se = FALSE returns one item list containing matrix otherwise
    ## list of two matrices "fit" and "se.fit"...

    if ( is.null(eta) ) {
      discrete <- is.list(X)
      lpi <- attr(X, "lpi")
      if (is.null(lpi)) lpi <- list(1 : ncol(X))

      K <- length(lpi) ## number of linear predictors --> no_eta
      nobs <- if ( discrete ) nrow(X$kd) else nrow(X)
      eta <- matrix(0, nobs, K)

      if ( se ) {
        ve <- matrix(0, nobs, K) ## variance of eta
        ce <- matrix(0, nobs, K * (K - 1)/2) ## covariance of eta_i eta_j
      }

      for ( i in 1 : K ) {
        if ( discrete ) eta[, i] <- Xbd(X$Xd, beta, k = X$kd, ks = X$ks, ts= X$ts, dt = X$dt, v = X$v, qc = X$qc, drop = X$drop, lt = X$lpid[[i]])
        else  {
          Xi <- X[, lpi[[i]], drop = FALSE]
          eta[, i] <- Xi %*% beta[lpi[[i]]] ## ith linear predictor
        }
        if ( !is.null(off[[i]]) ) eta[, i] <- eta[, i] + off[[i]]
        if ( se ) { ## variance and covariances for kth l.p.
          ve[, i] <- if ( discrete ) diagXVXd(X$Xd, Vb, k = X$kd, ks = X$ks, ts = X$ts, dt = X$dt, v = X$v,qc = X$qc, drop = X$drop, nthreads = 1, lt = X$lpid[[i]], rt = X$lpid[[i]])
                     else drop( pmax(0, rowSums((Xi %*% Vb[lpi[[i]], lpi[[i]]])* Xi)))
          ii <- 0
          if ( i < K ) for ( j in (i + 1) : K ) {
            ii <- ii + 1
            ce[, ii] <- if ( discrete ) diagXVXd(X$Xd, Vb, k = X$kd, ks = X$ks, ts = X$ts, dt = X$dt, v = X$v, qc = X$qc, drop = X$drop, nthreads = 1,
                                                 lt = X$lpid[[i]], rt = X$lpid[[j]]) else drop( pmax(0, rowSums((Xi %*% Vb[lpi[[i]],lpi[[j]]]) * X[,lpi[[j]]])))
          } # end for loop over ii
        } # end if se
      } # end for loop over i
    } else { # end if is.null eta
      se <- FALSE
    }

    out <- matrix(0, nrow(X), no_eta) # matrix of the predicted values
    cor_flag <- as.integer(getcflag()) #correlation by default
    if ( param == 1) internal()$pred_mcd(eta, out, d, cor_flag) #predicted values for the mcd
    if ( param == 2) internal()$pred_logm(eta, out, d, cor_flag) #predicted values for the logm

    if (se) { ## need to loop to find se of the mean vector, variances and covariances/correlations
      vp <- matrix(0, nrow(X), no_eta)
      for ( j in 1 : K ) {
        ## get dout_j/deta_k where dout_is a member of the mean vector or variances or correlations/covariances
        dout <- family$jacobian(eta, j)   #call to the jacobian function for the j-th component
        ## variances...
        vp[, j] <- rowSums(dout^2 * ve)
        ii <- 0
        for ( i in 1 : (K - 1) )
          for (k in (i + 1) : K) {
            ii <- ii + 1
            vp[, j] <- vp[, j] + 2 * dout[, i] * dout[, k] * ce[, ii]
          } # end for loop k
        vp[, j] <- sqrt( pmax(0, vp[, j]) ) ## transform to se
      } # end for loop j
      return(list(fit = out, se.fit = vp))
    } # end if se
    list(fit = out)
  } ## end mvn_mcd predict

  jacobian <- function(eta, jj){

    res <- matrix(0, nrow(eta), no_eta) # Jacobian matrix
    Cm <- internal()$mat2vec(d) #Matrix C (see the UK paper)

    cor_flag <-  as.integer(getcflag())

    if ( jj <= d ) {
      res[, jj] <- 1 # Mean vector case
    } else { # Covariance/Correlation matrix case
      jj <- jj - d
      idx_jj <- which(Cm == jj, arr.ind = TRUE) # identify the row and the column of the C matrix associated the linear predictor index
      S_r <- as.numeric(idx_jj[1, 1]) - 1 # extract the row
      S_c <- as.numeric(idx_jj[1, 2]) - 1 # extract the column

      if ( param == 1 ) { #mcd
        rc_idx_s <- rc_idx_t <- rep(NA, d * (d - 1)/2)
        count <- 1
        for ( j in (d + 1) : (d * (d + 1)/2) ) {  # identify the rows and the columns of the C matrix associated the linear predictor index
          rc_idx_s[count] <- as.numeric(which(Cm == j, arr.ind = TRUE)[1, 1]) - 1
          rc_idx_t[count] <- as.numeric(which(Cm == j, arr.ind = TRUE)[1, 2]) - 1
          count <- count + 1
        }
        internal()$jacobian_mcd(eta, res, d, S_r, S_c, rc_idx_s, rc_idx_t, cor_flag)
      }

      if ( param == 2 ) internal()$jacobian_logm(eta, res, d, S_r, S_c, cor_flag) #logm
    }
    return(res)
  } ## end jacobian

  rd <- function(mu, wt, scale) { # we need to consider the case of the logm
    ## simulate data given fitted linear predictor matrix in mu
    no_eta <- ncol(mu)
    d <- -3/2 + sqrt(9/4 + 2 * no_eta)
    out <- matrix(NA, nrow(mu), d)
    if ( param == 1 ) {
      for ( i in 1 : nrow(mu) ) {
        LD <- internal()$mcd_LD(mu[i,], d)
        C <- t(t(LD)/diag(LD))
        diag(C) <- diag(LD)
        u <- rmvn(1, rep(0, d), diag(rep(1, d)))
        out[i,] <- t(mu[i, 1 : d] + C %*% t(u))
      }
    }
    if ( param == 2 ) {
      for ( i in 1 : nrow(mu) ) {
        Sigma <- internal()$logM_Sigma(mu[i,], d)
        C <- t(chol(Sigma))
        u <- rmvn(1, rep(0, d), diag(rep(1, d)))
        out[i,] <- t(mu[i, 1 : d] + C %*% t(u))
      }
    }
    return(out)
  } ## rd

  structure(list(family = paste("Multivariate normal - Covariance modelling via", param2),
                 ll = ll, nlp = no_eta,
                 initialize = initialize,
                 getd = getd, putd = putd,
                 getno_eta = getno_eta, putno_eta = putno_eta,
                 getcflag = getcflag, put_cflag = putcflag,
                 getL1 = getL1, putL1 = putL1,
                 getL1_l = getL1_l, putL1_l = putL1_l,
                 getL2 = getL2, putL2 = putL2,
                 getL2_v = getL2_v, putL2_v = putL2_v,
                 getL2_l = getL2_l, putL2_l = putL2_l,
                 getL2_v_l = getL2_v_l, putL2_v_l = putL2_v_l,
                 getL3 = getL3, putL3 = putL3,
                 getL3_l = getL3_l, putL3_l = putL3_l,
                 getd1H = getd1H, putd1H = putd1H,
                 getd1eta = getd1eta, putd1eta = putd1eta,
                 getd1eta_l = getd1eta_l, putd1eta_l = putd1eta_l,
                 getV = getV, putV = putV,
                 getV_l = getV_l, putV_l = putV_l,
                 getidx_b = getidx_b, putidx_b = putidx_b,
                 getidx_aux =  getidx_aux, putidx_aux =  putidx_aux,
                 getidxl3 = getidxl3, putidxl3 = putidxl3,
                 getidxl3_no0 = getidxl3_no0, putidxl3_no0 = putidxl3_no0,
                 getidxl3_jkq = getidxl3_jkq, putidxl3_jkq = putidxl3_jkq,
                 getparam =  getparam, putparam =  putparam,
                 #postproc=postproc, ##to do???
                 rd = rd,
                 residuals = residuals,
                 predict = predict,
                 jacobian = jacobian,
                 linfo = stats, ## link information list
                 validmu = validmu,
                 d2link = 1, d3link = 1, d4link = 1, ## signals to fix.family.link that all done
                 ls = 1, ## signals that ls not needed here
                 available.derivs = a.der, ## can use full Newton here for the mcd parametrisation
                 discrete.ok = TRUE
  ),class = c("general.family", "extended.family", "family"))
} ## end mvn_scm
