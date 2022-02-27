#' Family definition
#'
#' @description  ...
#' @param d dimension of the outcome
#'
#' @return family
#' @export
#'
#' @importFrom Rcpp evalCpp
#'
#' @examples

mvn_mcd <- function(d = 2){
  #Multivariate model with covariance modelling via modified Cholesky decomposition
  if(d < 2) stop("mvn_mcd requires to or more dimensional data")
  no_eta <- d + d * (d + 1)/2
  stats <- list()
  for (j in 1:no_eta) stats[[j]] <- make.link("identity")

  validmu <- function(mu) all(is.finite(mu))

  assign(".d", d, envir = environment())
  getd <- function() get(".d")
  putd <- function(.x) assign(".d", .x, envir = environment(sys.function()))

  assign(".no_eta", no_eta, envir = environment())
  getno_eta <- function() get(".no_eta")
  putno_eta <- function(.x) assign(".no_eta", .x, envir = environment(sys.function()))

  initialize <- expression({
    my_init_fun <- function(y, nobs, E, x, family, offset){
      d <- family$getd() # d <- ncol(y)
      no_eta <- family$getno_eta()
      use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
      jj <- attr(x,"lpi")
      resid <-  matrix(0, nrow(y), d)
      if(is.list(x)){ # discrete case
        start <- rep(0,max(unlist(jj)))
        for(k in 1:d){
          yt1 <- y[,k]
          e1 <- E[ , jj[[k]], drop=FALSE] ## square root of total penalty
          R <- suppressWarnings(chol(XWXd(x$Xd, w = rep(1, length(yt1)),
                                          k = x$kd, ks = x$ks, ts = x$ts,
                                          dt = x$dt, v = x$v, qc = x$qc, nthreads = 1,
                                          drop = x$drop, lt = x$lpid[[k]]) +
                                     crossprod(E[ , jj[[k]]]), pivot = TRUE))
          Xty <- XWyd(x$Xd, rep(1, nrow(y)), yt1, x$kd,
                      x$ks, x$ts, x$dt, x$v, x$qc, x$drop, lt = x$lpid[[k]])
          piv <- attr(R, "pivot")
          rrank <- attr(R, "rank")
          startji <- rep(0, ncol(R))
          if (rrank < ncol(R)) {
            R <- R[1:rrank, 1:rrank]
            piv <- piv[1:rrank]
          }
          startMu <- numeric(ncol(R))
          startMu[piv] <- backsolve(R, forwardsolve(t(R), Xty[piv]))
          startMu[!is.finite(startji)] <- 0
          start[jj[[k]]] <- startMu
          eta1 <- Xbd(x$Xd, start, k = x$kd, ks = x$ks,
                      ts = x$ts, dt = x$dt, v = x$v, qc = x$qc, drop = x$drop,
                       lt = x$lpid[[k]])
          resid[,k] <-  y[,k] -   eta1
        }
        Svcov <- cov(resid)
        MCD_dec <- mcd(Svcov)
        Theta_el<-c(diag(MCD_dec), MCD_dec[upper.tri(MCD_dec, diag=FALSE)])

        #for(k in (d+1):no_eta){
          #x1 <-  x[ , jj[[k]],drop=FALSE]
          #startji <- qr.coef(qr(x1), c(rep(Theta_el[k-d],nrow(x1))))
          #startji[!is.finite(startji)] <- 0
          #start[jj[[k]]] <- startji
        #}
      } else { #regular case
        start <- rep(0,ncol(x))
        for(k in 1:d){
          yt1 <- y[,k]
          x1 <- x[ , jj[[k]], drop=FALSE]
          e1 <- E[ , jj[[k]], drop=FALSE] ## square root of total penalty

          if(use.unscaled){
            x1 <- rbind(x1, e1)
            startMu <- qr.coef(qr(x1), c(yt1,rep(0,nrow(E))))
            startMu[ !is.finite(startMu) ] <- 0
          } else {
            startMu <- pen.reg(x1, e1, yt1)
          }

          resid[,k] <-  y[,k] - x[ , jj[[k]]]%*%startMu
          start[jj[[k]]] <- startMu
        }

        Svcov <- cov(resid)
        MCD_dec <- mcd(Svcov)
        Theta_el<-c(diag(MCD_dec), MCD_dec[upper.tri(MCD_dec, diag=FALSE)])

        for(k in (d+1):no_eta){
          x1 <-  x[ , jj[[k]],drop=FALSE]
          startji <- qr.coef(qr(x1), c(rep(Theta_el[k-d],nrow(x1))))
          startji[!is.finite(startji)] <- 0
          start[jj[[k]]] <- startji
        }
      }
      return(start)
    }
    if(is.null(start)){
      start <- my_init_fun(y = y, nobs = nobs, E = E, x = x, family = family, offset = offset)
    }
  }) ## initialize

  ##To do residuals and postproc
  residuals <- function(object,type=c("response","deviance")) {
    type <- match.arg(type)
    res <- object$y[,1:d] - object$fitted.values[,1:d]
    res
  } ## residuals

  trind.generator <- function (K = 2, deriv = 0) {
    m.start <- 1
    i2 <- array(0, dim = c(K, K))
    m <- m.start
    for (k in 1:K) for (l in k:K) {
      i2[k, l] <- i2[l, k] <- m
      m <- m + 1
    }

    m <- m.start
    i2r <- rep(0, max(i2))
    for (k in 1:K) for (l in k:K) {
      i2r[m] <- l + (k - 1) * K
      m <- m + 1
    }
    if(deriv == 0){
      i3 <- i3r <- i4 <- i4r <- NULL
    } else {
      i3 <- array(0, dim = c(K, K, K))
      m <- m.start
      for (j in 1:K) for (k in j:K) for (l in k:K) {
        i3[j, k, l] <- i3[j, l, k] <- i3[k, l, j] <- i3[k, j,
                                                        l] <- i3[l, j, k] <- i3[l, k, j] <- m
        m <- m + 1
      }
      m <- m.start
      i3r <- rep(0, max(i3))
      for (j in 1:K) for (k in j:K) for (l in k:K) {
        i3r[m] <- l + (k - 1) * K + (j - 1) * K^2
        m <- m + 1
      }
    i4 <- array(0, dim = c(K, K, K, K))
    m <- m.start
    for (i in 1:K) for (j in i:K) for (k in j:K) for (l in k:K) {
      i4[i, j, k, l] <- i4[i, j, l, k] <- i4[i, k, l, j] <- i4[i,
                                                               k, j, l] <- i4[i, l, j, k] <- i4[i, l, k, j] <- i4[j,
                                                                                                                  i, k, l] <- i4[j, i, l, k] <- i4[j, k, l, i] <- i4[j,
                                                                                                                                                                     k, i, l] <- i4[j, l, i, k] <- i4[j, l, k, i] <- i4[k,
                                                                                                                                                                                                                        j, i, l] <- i4[k, j, l, i] <- i4[k, i, l, j] <- i4[k,
                                                                                                                                                                                                                                                                           i, j, l] <- i4[k, l, j, i] <- i4[k, l, i, j] <- i4[l,
                                                                                                                                                                                                                                                                                                                              j, k, i] <- i4[l, j, i, k] <- i4[l, k, i, j] <- i4[l,
                                                                                                                                                                                                                                                                                                                                                                                 k, j, i] <- i4[l, i, j, k] <- i4[l, i, k, j] <- m
      m <- m + 1
    }
    m <- m.start
    i4r <- rep(0, max(i4))
    for (i in 1:K) for (j in i:K) for (k in j:K) for (l in k:K) {
      i4r[m] <- l + (k - 1) * K + (j - 1) * K^2 + (i - 1) * K^3
      m <- m + 1
     }
    }
    list(i2 = i2, i3 = i3, i4 = i4, i2r = i2r, i3r = i3r, i4r = i4r)
  }


  ll <- function(y, X, coef, wt, family,  offset = NULL,
                 deriv = 0, d1b = NULL, d2b = NULL, Hp = NULL,
                 rank = 0, fh = NULL, D = NULL) {
    discrete <- is.list(X)
    jj <- attr(X,"lpi") ## extract linear predictor index

    n <- nrow(y)
    eta <- matrix(0, n, no_eta)
    for(k in 1 : no_eta){
      eta[,k] <- if (discrete) Xbd(X$Xd,coef,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[k]])
      else X[,jj[[k]],drop=FALSE]%*%coef[jj[[k]]]
    }

    l1 <- matrix(0, n, no_eta) #initialization
    ## log-likelihood: eta is a matrix n*w and y is a matrix n*d

    l <- ll_mcd(eta, y[,1:d]) - 0.5 * n * d * log(2 * pi)


    if (deriv>0) {
      ## the first derivative: eta is a matrix n*w, y is a matrix n*d,
      ##                         l1 is a matrix n*w
      d1_mcd(eta, y[,1:d], l1)

      ## the second derivatives
      l2 <- matrix(0, n, no_eta * (no_eta + 1)/2) #initialization
      d2_mcd(eta,y[,1:d], l2)

     }

     l3 <- 0 ## defaults

     if (deriv) {
       i2 <- family$tri$i2
       i3 <- family$tri$i3

       ## get the gradient and Hessian...
       ret <- gamlss.gH(X, jj, l1, l2, i2, l3 = l3, i3 = i3,
                        d1b = d1b, deriv = deriv - 1, fh = fh, D = D)
     } else ret <- list()
    ret$l <- l
    ret
   } ## end ll mvn_mcd

  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                      beta=NULL,off=NULL,Vb=NULL) {
    ## optional function to give predicted values - idea is that
    ## predict.gam(...,type="response") will use this, and that
    ## either eta will be provided, or {X, beta, off, Vb}. family$data
    ## contains any family specific extra information.
    ## if se = FALSE returns one item list containing matrix otherwise
    ## list of two matrices "fit" and "se.fit"...


    if (is.null(eta)) {
      discrete <- is.list(X)
      lpi <- attr(X,"lpi")
      if (is.null(lpi)) {
        lpi <- list(1:ncol(X))
      }
      K <- length(lpi) ## number of linear predictors
      nobs <- if (discrete) nrow(X$kd) else nrow(X)
      eta <- matrix(0,nobs,K)
      for (i in 1:K) {
        if (discrete) {
          eta[,i] <- Xbd(X$Xd,beta,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[i]])
        } else {
          Xi <- X[,lpi[[i]],drop=FALSE]
          eta[,i] <- Xi%*%beta[lpi[[i]]] ## ith linear predictor
        }
        if (!is.null(off[[i]])) eta[,i] <- eta[,i] + off[[i]]
      }
    }

    out <- matrix(0, nrow(eta), ncol(eta))
    pred_mcd(eta, out,d)
    list(fit = out)
  } ## mvncm predict

  jacobian <- function(eta, jj){
    #The following two lines could be unuseful
    #no_eta <- ncol(eta) # delete in SCM
    #d <- -3/2 + sqrt(9/4 + 2 * no_eta) # delete in SCM
    res <- matrix(0, nrow(eta), no_eta)
    G <- mat2vec(d) #put in cpp???

    if(jj <= d){
      res[,jj] <- 1
    } else {
      jj <- jj - d
      idx_jj <- which(G == jj, arr.ind = TRUE)
      S_row <- as.numeric(idx_jj[1, 1])
      S_col <- as.numeric(idx_jj[1, 2])

      rc_idx_s<-rep(NA, d*(d-1)/2)
      rc_idx_s<-rep(NA, d*(d-1)/2)

      count <- 1
      for(j in (d + 1):(d * (d + 1)/2)){
        rc_idx_s[count] <- which(G == j, arr.ind=TRUE)[1, 1]
        rc_idx_t[count] <- which(G == j, arr.ind=TRUE)[1, 2]
        count <- count + 1
      }

      jacobian_mcd(eta, res, d, S_row - 1, S_col - 1,
                   as.numeric(rc_idx_s) - 1, as.numeric(rc_idx_t) - 1)
    }
    return(res)
  }

  structure(list(family = "Multivariate normal (MCD)", ll = ll, nlp = no_eta,
                 ##link=paste(link), ## unuseful?
                 tri = trind.generator(no_eta, deriv = 0), ## symmetric indices for accessing derivative arrays
                 initialize = initialize,
                 getd=getd, putd=putd,
                 getno_eta=getno_eta, putno_eta=putno_eta,
                 #postproc=postproc, ##to do
                 residuals=residuals,
                 predict = predict,
                 jacobian = jacobian,
                 linfo = stats, ## link information list
                 validmu = validmu,
                 #rd=rd,
                 d2link = 1, d3link = 1, d4link = 1, ## signals to fix.family.link that all done
                 ls = 1, ## signals that ls not needed here
                 available.derivs = 0, ## can use full Newton here
                 discrete.ok = TRUE
                 ),class = c("general.family", "extended.family", "family"))
} ## end mvn_cmd

