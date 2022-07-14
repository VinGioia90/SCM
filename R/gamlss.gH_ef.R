#' Modified gamlss.gH function
#' @description This function allows to obtain an ordering of importance of the effects on the unconstrained elements of the modified cholesky decomposition
#'
#' @param X model matrix covariates
#' @param jj lpi index
#' @param eta lpi
#' @param y outcome
#' @param i2 idx der2
#' @param n_block number of blocks partitioning the observations
#' @param z auxiliary index vector z
#' @param w auxiliary index vector w
#' @param Gmat auxiliary index matrix G (not used here)
#' @param t auxiliary index vector t (not used here)
#' @param idx_jk list of zero indices
#' @param l3 3 - derivatives w.r.t. eta
#' @param i3 idx der3
#' @param l2 initialization l2
#' @param l2_last initialization l2_last
#' @param idx_block indici block
#' @param l4 4 - derivatives w.r.t. eta
#' @param i4 idx der4
#' @param d1b d1b??
#' @param d2b d2b??
#' @param deriv numeber of derivatives
#' @param fh fh??
#' @param D D??
#'
#' @return A list of summaries
#' @export
#'
#' @importFrom Rcpp evalCpp
#'
#' @examples
#'
gamlss.gH_ef <- function (X, jj, eta, y, i2, n_block = 1, w, z, t, Gmat, idx_jk, l3 = 0, i3 = 0, l1,  l2, l2_last, idx_block,
                          l4 = 0, i4 = 0, d1b = 0, d2b = 0, deriv = 0, fh = NULL, D = NULL){
  K <- length(jj) # no_eta
  sparse <- discrete <- FALSE
  p <- ncol(X)
  n <- nrow(X)

  trHid2H <- d1H <- d2H <- NULL
  ifunc <- !is.array(i2)
  #####################
  # Score w.r.t. beta #
  #####################
  lb <- rep(0, p)
  # Chiamata funzione cpp per ottenere lb (derivate parziali in beta)
  #print("ciao")
  #print(microbenchmark({
   d1_mcd_beta(X,lapply(jj,function(x) x -1),K,lb,l1, eta,y, z, w, Gmat)
  #}, times=2L))

  #print("saluti")
  #######################
  # Hessian w.r.t. beta #
  #######################
  lbb <- if (sparse)
    Matrix(0, p, p)
  else matrix(0, p, p)

  # Spostare fuori da gamlss.GH_ef
  #nel_Heta<- as.integer(length(unlist(idx_jk)))
  #nlast<- n%%n_block
  #nset<- n%/%n_block
  #if(nlast == 0){
  #  l2 <- matrix(0, nset, nel_Heta) #si può tirare fuori dall'if
  #  l2_last <- matrix(0,nset,nel_Heta)
  #  idx_block <- c(-1,cumsum(rep(nset, n_block-1))-1, n-1)
  #} else {
  #  l2 <- matrix(0, nset, nel_Heta) #si può tirare fuori dall'if
  #  l2_last <- matrix(0, nlast,nel_Heta)
  #  idx_block <- c(-1,cumsum(rep(nset, n_block))-1, n-1)
  #}
  #print("ciao")
  #print(microbenchmark(d2_mcd_beta(X,eta,y,lapply(jj,function(x) x -1),K,idx_jk,lbb,l2, l2_last, idx_block, z, w, Gmat, t), times=10L))
  d2_mcd_beta(X,eta,y,lapply(jj,function(x) x -1),K,idx_jk,lbb,l2, l2_last, idx_block, z, w, Gmat, t)
  #print("saluti")
  ###############################
  # Parte da lasciare invariata #
  ###############################
  if (deriv > 0) {
    stop("error")
  }

  list(lb = lb, lbb = lbb, d1H = d1H, d2H = d2H, trHid2H = trHid2H)
}
