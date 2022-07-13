#' Family definition (internal use)
#'
#' @description  ...
#' @param d dimension of the outcome
#' @param order_AD order of the ante dependence model
#'
#' @return family
#' @export
#'
#' @importFrom Rcpp evalCpp
#'
#' @examples

mvn_mcd_AD <- function(d = 2, order_AD=1){
  #Multivariate model with covariance modelling via modified Cholesky decomposition
  if(d < 2) stop("mvn_mcd requires to or more dimensional data")

    no_eta <- d + d *(d + 1)/2

  if(order_AD >= 0 & order_AD < d) {
    #no_eta_noc <- d + d *(d + 1)/2 - (d-order_AD)*(d-order_AD-1)/2
    constraint_AD<-function(d,order_AD){
      # only if orderAD < d-1
      res<-matrix(0,(d-order_AD)*(d-order_AD-1)/2,2)
      count<-1
      for(j in (2+order_AD):d){
        for(k in 1:(j-order_AD-1)){
          res[count,1] <- j
          res[count,2] <- k
          count <- count + 1
        }
      }
      return(res)
    }
    CAD <- constraint_AD(d, order_AD)

    constrlpi <- rep(0, dim(CAD)[1])
    for(k in 1:length(constrlpi)){
      constrlpi[k]<-   (mat2vec(d)+d)[CAD[k,1], CAD[k,2]]
    }

    unconstr_lpi <- (1:no_eta)[!((1:no_eta) %in% constrlpi)]
    assign(".constrlpi", constrlpi, envir = environment())
    getconstrlpi <- function() get(".constrlpi")
    putconstrlpi <- function(.constrlpi) assign(".constrlpi", .constrlpi, envir = environment(sys.function()))
    assign(".unconstr_lpi", unconstr_lpi, envir = environment())
    getunconstr_lpi <- function() get(".unconstr_lpi")
    putunconstr_lpi <- function(.unconstr_lpi) assign(".unconstr_lpi", .unconstr_lpi, envir = environment(sys.function()))

  } else {
    print("The order of the ante-dependence model lies in [0,d]")
  }

  stats <- list()
  for (j in 1:(d + d *(d + 1)/2 - (d-order_AD)*(d-order_AD-1)/2)) stats[[j]] <- make.link("identity")

  validmu <- function(mu) all(is.finite(mu))

  assign(".cflag", TRUE, envir = environment())
  getcflag <- function() get(".cflag")
  putcflag <- function(.cflag) assign(".cflag", .cflag, envir = environment(sys.function()))


  assign(".d", d, envir = environment())
  getd <- function() get(".d")
  putd <- function(.d) assign(".d", .d, envir = environment(sys.function()))


  assign(".order_AD", order_AD, envir = environment())
  getorder_AD <- function() get(".order_AD")
  putorder_AD <- function(.order_AD) assign(".order_AD", .order_AD, envir = environment(sys.function()))

  getL2 <- function() get(".l2")
  putL2 <- function(.l2) assign(".l2", .l2, envir = environment(sys.function()))

  assign(".no_eta", no_eta, envir = environment())
  getno_eta <- function() get(".no_eta")
  putno_eta <- function(.no_eta) assign(".no_eta", .no_eta, envir = environment(sys.function()))

  initialize <- expression({
    my_init_fun <- function(y, nobs, E, x, family, offset){
      d <- family$getd() # d <- ncol(y)
      order_AD <- family$getorder_AD() # d <- ncol(y)
      no_eta <- family$getno_eta()
      unconstr_lpi <- family$getunconstr_lpi()
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

        count<-1
        for(k in (d+1):no_eta){
        #for(k in (d+1):no_eta){
          if(k %in% unconstr_lpi){
          x1 <-  x[ , jj[[count+d]],drop=FALSE]
          startji <- qr.coef(qr(x1), c(rep(Theta_el[k-d],nrow(x1))))
          startji[!is.finite(startji)] <- 0
          start[jj[[count+d]]] <- startji
          count<-count+1
          }
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
    #if(type == "deviance"){
    #  n <- dim(object$fitted.values)[1]
    #  res <- matrix(0, n, d)
    #  res_dev_mcd(object$fitted.values, object$y, res)
    #} else {
      res <- object$y - object$fitted.values[,1:d]
    #}
    res
  } ## residuals


  ll <- function(y, X, coef, wt, family,  offset = NULL,
                 deriv = 0, d1b = NULL, d2b = NULL, Hp = NULL,
                 rank = 0, fh = NULL, D = NULL) {
    discrete <- is.list(X)
    jj <- attr(X,"lpi") ## extract linear predictor index

    n <- nrow(y)
    eta <- matrix(0, n, no_eta)

    l2 <- try(getL2(), TRUE)
    if("try-error" %in% class(l2)){
      l2 <-  matrix(0, n, no_eta * (no_eta + 1)/2) #initialization
      putL2(l2)
    }

    count<-1
    for(k in 1 : no_eta){
      if(k %in% unconstr_lpi){
        eta[,k] <- if (discrete) Xbd(X$Xd,coef,k=X$kd,ks=X$ks,ts=X$ts,dt=X$dt,v=X$v,qc=X$qc,drop=X$drop,lt=X$lpid[[count]])
        else X[,jj[[count]],drop=FALSE]%*%coef[jj[[count]]]
        count<- count + 1
      } else {
        eta[,k] <- rep(0, n)
      }
    }

    l1 <- matrix(0, n, no_eta) #initialization
    ## log-likelihood: eta is a matrix n*w and y is a matrix n*d

    #l <- ll_mcd(eta, y[,1:d]) - 0.5 * n * d * log(2 * pi)
    l <- ll_mcd(eta, y) - 0.5 * n * d * log(2 * pi)


    if (deriv>0) {
      ## the first derivative: eta is a matrix n*w, y is a matrix n*d,
      ##                         l1 is a matrix n*w
      #d1_mcd(eta, y[,1:d], l1)
      d1_mcd(eta, y, l1)

      ## the second derivatives
      #l2 <- matrix(0, n, no_eta * (no_eta + 1)/2) #initialization
      #d2_mcd(eta,y[,1:d], l2)
      d2_mcd(eta,y, l2)
      constrlpi <- family$getconstrlpi()
      #idx_constr_l2 <- rep(0,sum(constrlpi))
      lclpi<-length(constrlpi)
      nT<-d*(d-1)/2
      idx_constr_l2 <- rep(0,2*d*lclpi+ (nT*(nT+1)/2-(nT-lclpi)*(nT-lclpi+1)/2))
      length(idx_constr_l2)
      count <- 1
      count2 <- 1

      for(k in 1:length(constrlpi)){
        count <- 1
        for(i in 1:no_eta){
          for( j in i:no_eta){
          if(i != constrlpi[k] & j == constrlpi[k] & i != j){
              idx_constr_l2[count2] <- count
              count2 <- count2+1
            }
            if(i==constrlpi[k]){
              idx_constr_l2[count2] <- count
              count2 <- count2+1
            }
            count <- count + 1
          }
        }
      }
    }


    l3 <- 0 ## defaults
    if (deriv) {
      i2 <- family$tri$i2
      i3 <- family$tri$i3

      ## get the gradient and Hessian...
      ret <- gamlss.gH(X, jj, l1[,-constrlpi], l2[,-sort(unique(sort(idx_constr_l2)))], i2, l3 = l3, i3 = i3,
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


    jj <- attr(X,"lpi") ## extract linear predictor index
    #if (is.null(eta)) {
    #  discrete <- is.list(X)
    #  lpi <- attr(X,"lpi")
    #  if (is.null(lpi)) {
    #    lpi <- list(1:ncol(X))
    #  }
      #K <- length(lpi) ## number of linear predictors
    nobs <- nrow(X)
      eta <- matrix(0,nobs,no_eta)

      count<-1
      for(k in 1 : no_eta){
        if(k %in% unconstr_lpi){
          eta[,k] <- X[,jj[[count]],drop=FALSE]%*%beta[jj[[count]]]
          count<- count + 1
        } else {
          eta[,k] <- rep(0, nobs)
        }
      }



    out <- matrix(0, nrow(eta), ncol(eta))
    cor_flag <- getcflag()
    pred_mcd(eta, out,d, as.integer(cor_flag))
    list(fit = out)
  } ## mvncm predict
  #
  # jacobian <- function(eta, jj){
  #   #The following two lines could be unuseful
  #   #no_eta <- ncol(eta) # delete in SCM
  #   #d <- -3/2 + sqrt(9/4 + 2 * no_eta) # delete in SCM
  #   res <- matrix(0, nrow(eta), no_eta)
  #   G <- mat2vec(d) #put in cpp???
  #
  #   cor_flag <- getcflag()
  #
  #   if(jj <= d){
  #     res[,jj] <- 1
  #   } else {
  #     jj <- jj - d
  #     idx_jj <- which(G == jj, arr.ind = TRUE)
  #     S_row <- as.numeric(idx_jj[1, 1])
  #     S_col <- as.numeric(idx_jj[1, 2])
  #
  #     rc_idx_s <- rep(NA, d * (d - 1)/2)
  #     rc_idx_t <- rep(NA, d * (d - 1)/2)
  #
  #     count <- 1
  #     for(j in (d + 1):(d * (d + 1)/2)){
  #       rc_idx_s[count] <- which(G == j, arr.ind=TRUE)[1, 1]
  #       rc_idx_t[count] <- which(G == j, arr.ind=TRUE)[1, 2]
  #       count <- count + 1
  #     }
  #
  #     jacobian_mcd(eta, res, d, S_row - 1, S_col - 1,
  #                  as.numeric(rc_idx_s) - 1, as.numeric(rc_idx_t) - 1, as.integer(cor_flag))
  #   }
  #   return(res)
  # }

  structure(list(family = "Multivariate normal (MCD - AD(r))", ll = ll, nlp = d + d *(d + 1)/2 - (d-order_AD)*(d-order_AD-1)/2,
                 ##link=paste(link), ## unuseful?
                 tri = trind.generator(d + d *(d + 1)/2 - (d-order_AD)*(d-order_AD-1)/2,ifunc=TRUE),
                 initialize = initialize,
                 getd = getd, putd = putd,
                 getorder_AD = getorder_AD, putorder_AD = putorder_AD,
                 getno_eta = getno_eta, putno_eta = putno_eta,
                 getconstrlpi = getconstrlpi, putconstrlpi = putconstrlpi,
                 getunconstr_lpi = getunconstr_lpi, putunconstr_lpi = putunconstr_lpi,
                 getcflag = getcflag, put_cflag = putcflag,
                 getL2 = getL2, putL2 = putL2,
                 #postproc=postproc, ##to do
                 residuals=residuals,
                 predict = predict,
                 #jacobian = jacobian,
                 linfo = stats, ## link information list
                 validmu = validmu,
                 #rd=rd,
                 d2link = 1, d3link = 1, d4link = 1, ## signals to fix.family.link that all done
                 ls = 1, ## signals that ls not needed here
                 available.derivs = 0, ## can use full Newton here
                 discrete.ok = TRUE
  ),class = c("general.family", "extended.family", "family"))
} ## end mvn_cmd
