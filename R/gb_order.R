#' Perform the gradient boosting
#' @description This function allows to obtain an ordering of importance of the effects on the unconstrained elements of the modified cholesky decomposition
#'
#' @param y outcome (residual)
#' @param X covariates
#' @param effects list of effects included in the boosting procedure
#' @param constraints list of constranints on the effects
#' @param nstep numbers of step of the boosting procedure
#' @param verbose to print the iteration summary of the procedure
#' @param rate learning rate
#' @param ncores number of cores for parallel computation
#'
#' @return A list of summaries
#' @export
#' @importFrom parallel mclapply
#' @importFrom mgcv gam
#' @importFrom mgcv predict.gam
#' @importFrom stats as.formula
#' @importFrom stats cov
#'
#' @examples
#'


boost_order <- function(y, X, effects, constraints,
                              nstep = 1, verbose = FALSE,
                              rate = 1, ncores = 1){
  d <- ncol(y)
  n <- nrow(y)
  neta <- d * (d + 1)/2
  neff <- length(effects)

  # Intercept initialization
  Svcov <- cov(y)

  MCD_dec <- mcd(Svcov)
  Theta_el<-c(diag(MCD_dec), MCD_dec[upper.tri(MCD_dec, diag=FALSE)])

  score_init <- matrix(0, n, neta + d)
  eta <- matrix(0, n, neta + d)

  for(j in 1:(d*(d + 1)/2)){
    eta[,j+d]  <- rep(Theta_el[j], n)
  }

  l0 <- ll_mcd(eta, y)

  dll <- list()
  idx <- list()

  for(ii in 1:nstep){
    d1_mcd(eta, y, score_init)

    dll[[ii]] <- mclapply(effects, function(eff){
      delta <- rep(0, neta)
      index <- constraints[[eff]]
      if( is.null(index) ) index <- 1:neta
      for(jj in index){
        form_gam <- as.formula(paste0("score_init[,jj + d] ~ ", eff))
        fit <- try( gam(form_gam, data=X) )
        if( inherits(fit, "try-error") ){
          delta[jj] <- -Inf
        } else {
          # Improve the following steps: avoid creating a copy
          eta1 <- eta
          eta1[,jj+d] <- eta1[ , jj+d] + rate * predict.gam(fit)
          delta[jj] <- ll_mcd(eta1, y) - l0
        }
      }
      return(delta)
    }, mc.cores = ncores)

    dll[[ii]] <- do.call("cbind", dll[[ii]])

    l0 <- l0 + max(dll[[ii]])

    idx[[ii]] <- which(dll[[ii]] == max(dll[[ii]]), arr.ind = TRUE)

    yyy <- score_init[, idx[[ii]][1, 1] + d]
    form_gam <- as.formula(paste0("yyy ~ ", effects[idx[[ii]][1, 2]]))

    eta[, d + idx[[ii]][1,1]] <- eta[, d + idx[[ii]][1,1]] + rate * predict.gam(gam(form_gam, data=X))

    if(verbose == TRUE){
      print(paste("iter" = ii, "delta" = max(dll[[ii]]),
                  "id" = length(unique(sapply(idx, function(.x) paste0(.x[1,1], ".", .x[1,2]))))))
      }
    if(max(dll[[ii]]) < 0) break
  }

  return(list(dll = dll, idx = idx, d = d, nstep = nstep, effects = effects))
}




