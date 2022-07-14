#' Family definition (internal use)
#'
#' @description  ...
#' @param d dimension of the outcome
#' @param nb number of blocks
#'
#' @return family
#' @export
#'
#' @importFrom Rcpp evalCpp
#'
#' @examples

mvn_mcd_ef <- function(d = 2, nb=1){
  #Multivariate model with covariance modelling via modified Cholesky decomposition
  if(d < 2) stop("mvn_mcd requires to or more dimensional data")
  no_eta <- d + d *(d + 1)/2

   z<-rep(0,(d*(d - 1)/2))
   w<-rep(0,(d*(d - 1)/2))
   t<-rep(0,(d*(d - 1)/2))
   Gmat<-matrix(0,d - 1,d - 1)
   mode(Gmat) <- mode(z) <- mode(w) <- mode(t) <- "integer"
   idx_zwGt(d,z,w,Gmat,t)
   idx_jk <- index0s(no_eta, z,w,Gmat,t)


   nel_Heta<- as.integer(length(unlist(idx_jk)))

   getL1 <- function() get(".l1")
   putL1 <- function(.l1) assign(".l1", .l1, envir = environment(sys.function()))

   getL2 <- function() get(".l2")
   putL2 <- function(.l2) assign(".l2", .l2, envir = environment(sys.function()))

   getL2_last <- function() get(".l2_last")
   putL2_last <- function(.l2_last) assign(".l2_last", .l2_last, envir = environment(sys.function()))

   getidx_block <- function() get(".idx_block")
   putidx_block <- function(.idx_block) assign(".idx_block", .idx_block, envir = environment(sys.function()))

  stats <- list()
  for (j in 1:no_eta) stats[[j]] <- make.link("identity")

  validmu <- function(mu) all(is.finite(mu))

  assign(".cflag", TRUE, envir = environment())
  getcflag <- function() get(".cflag")
  putcflag <- function(.cflag) assign(".cflag", .cflag, envir = environment(sys.function()))


  assign(".d", d, envir = environment())
  getd <- function() get(".d")
  putd <- function(.d) assign(".d", .d, envir = environment(sys.function()))


  assign(".no_eta", no_eta, envir = environment())
  getno_eta <- function() get(".no_eta")
  putno_eta <- function(.no_eta) assign(".no_eta", .no_eta, envir = environment(sys.function()))

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
          if(!is.matrix(x[ , jj[[k]]])){
            resid[,k] <-  y[,k] - x[ , jj[[k]]]*startMu
          } else{
            resid[,k] <-  y[,k] - x[ , jj[[k]]]%*%startMu
          }
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
  residuals <- function(object, type=c("response","deviance")) {
    type <- match.arg(type)
    if(type == "deviance"){
      n <- dim(object$fitted.values)[1]
      res <- matrix(0, n, d)
      res_dev_mcd(object$fitted.values, object$y, res)
    } else {
      res <- object$y - object$fitted.values[,1:d]
    }
    res
  } ## residuals


  ll <- function(y, X, coef, wt, family,  offset = NULL,
                 deriv = 0, d1b = NULL, d2b = NULL, Hp = NULL,
                 rank = 0, fh = NULL, D = NULL) {
    discrete <- is.list(X)
    jj <- attr(X,"lpi") ## extract linear predictor index

    n <- nrow(y)
    eta <- matrix(0, n, no_eta)

    l1 <- try(getL1(), TRUE)
    if("try-error" %in% class(l1)){
      l1 <-   matrix(0, n, no_eta) #initialization
      putL1(l1)
    }



    nlast<- n%%nb
    nset<- n%/%nb

    l2 <- try(getL2(), TRUE)
    if("try-error" %in% class(l2)){
      l2 <-   matrix(0, nset, nel_Heta) #initialization
      putL2(l2)
    }

    if(nlast == 0){
      l2_last <- try(getL2_last(), TRUE)
      if("try-error" %in% class(l2_last)){
        l2_last <- matrix(0,nset,nel_Heta)
        putL2_last(l2_last)
      }

      idx_block <- try(getidx_block(), TRUE)
      if("try-error" %in% class(idx_block)){
        idx_block <- c(-1,cumsum(rep(nset, nb-1))-1, n-1)
        putidx_block(idx_block)
      }
    } else {
      l2_last <- try(getL2_last(), TRUE)
      if("try-error" %in% class(l2_last)){
        l2_last <- matrix(0, nlast,nel_Heta)
        putL2_last(l2_last)
      }

      idx_block <- try(getidx_block(), TRUE)
      if("try-error" %in% class(idx_block)){
        idx_block <- c(-1,cumsum(rep(nset, nb))-1, n-1)
        putidx_block(idx_block)
      }

    }


    for(k in 1 : no_eta){
      eta[,k] <- if (discrete) Xbd(X$Xd,coef,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[k]])
      else X[,jj[[k]],drop=FALSE]%*%coef[jj[[k]]]
    }

    ## log-likelihood: eta is a matrix n*q and y is a matrix n*d
    l <- ll_mcd(eta, y) - 0.5 * n * d * log(2 * pi)

    l3 <- 0 ## defaults

    if (deriv) {
      i2 <- family$tri$i2
      i3 <- family$tri$i3

      ## get the gradient and Hessian...
     #print("ciao")
     #print(microbenchmark( ret = gamlss.gH_ef(X, jj, eta, y, i2, n_block=nb,   w, z, t, Gmat, idx_jk,   l3 = l3, i3 = i3, l2, l2_last, idx_block,
     #                                     d1b = d1b, deriv = deriv - 1, fh = fh, D = D), times=10L)) # si può rimuovere il passaggio di i2 (non più utilizzzato))
     #print("saluti")
      ret <- gamlss.gH_ef(X, jj, eta, y, i2, n_block=nb,   w, z, t, Gmat, idx_jk,   l3 = l3, i3 = i3, l1, l2, l2_last, idx_block,
                        d1b = d1b, deriv = deriv - 1, fh = fh, D = D) # si può rimuovere il passaggio di i2 (non più utilizzzato)

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
    cor_flag <- getcflag()
    pred_mcd(eta, out,d, as.integer(cor_flag))
    list(fit = out)
  } ## mvncm predict

  jacobian <- function(eta, jj){
    #The following two lines could be unuseful
    #no_eta <- ncol(eta) # delete in SCM
    #d <- -3/2 + sqrt(9/4 + 2 * no_eta) # delete in SCM
    res <- matrix(0, nrow(eta), no_eta)
    G <- mat2vec(d) #put in cpp???

    cor_flag <- getcflag()

    if(jj <= d){
      res[,jj] <- 1
    } else {
      jj <- jj - d
      idx_jj <- which(G == jj, arr.ind = TRUE)
      S_row <- as.numeric(idx_jj[1, 1])
      S_col <- as.numeric(idx_jj[1, 2])

      rc_idx_s <- rep(NA, d * (d - 1)/2)
      rc_idx_t <- rep(NA, d * (d - 1)/2)

      count <- 1
      for(j in (d + 1):(d * (d + 1)/2)){
        rc_idx_s[count] <- which(G == j, arr.ind=TRUE)[1, 1]
        rc_idx_t[count] <- which(G == j, arr.ind=TRUE)[1, 2]
        count <- count + 1
      }

      jacobian_mcd(eta, res, d, S_row - 1, S_col - 1,
                   as.numeric(rc_idx_s) - 1, as.numeric(rc_idx_t) - 1, as.integer(cor_flag))
    }
    return(res)
  }

  structure(list(family = "Multivariate normal (MCD)", ll = ll, nlp = no_eta,
                 ##link=paste(link), ## unuseful?
                 tri = trind.generator(no_eta,ifunc=FALSE),
                 initialize = initialize,
                 getd = getd, putd = putd,
                 getno_eta = getno_eta, putno_eta = putno_eta,
                 getcflag = getcflag, put_cflag = putcflag,
                 #get_z=get_z, put_z=put_z,get_w=get_w, put_w=put_w,
                 #get_t=get_t, put_t=put_t, get_G=get_G, put_G=put_G,
                 #getidx_jk = getidx_jk, putidx_jk = putidx_jk,
                 getL2 = getL2, putL2 = putL2,
                 getL1 = getL1, putL1 = putL1,
                 getL2_last = getL2_last, putL2_last = putL2_last,
                 getidx_block = getidx_block, putidx_block = putidx_block,
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
} ## end mvn_cmd_ef
