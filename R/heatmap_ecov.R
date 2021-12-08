#' Heatmap Empirical Covariance, Precision and Decomposition Matrix.
#' @description plotECPD allows to plot the heatmap of the covariance,  precision and the decomposition matrix. The decomposition matrix  is obtained by using the modified cholesky decomposition or the matrix logarithm approach.
#'
#' @param ecov empirical covariance matrix
#' @param type the type of matrix to be plotted. It is allowed to plot the covariance matrix by type = "cov", the precision matrix by type = "prec" and the elements of the decomposition by type="dec". By default type = "cov".
#' @param repar the parameterization considered if type = "dec". The modified Cholesky decomposition is obtained with repar = "mcd" and the matrix logarithm with repar = "logm".
#' @param ord an optional vector that allows to specify the order of the variables.
#' @param lab an optional vector that allows to specify the label of the variables.
#' @param neimat an optional matrix to be used to marks neighbourhood variables.
#' @param di logical. If FALSE (default) the diagonal elements are not showed.
#' @param lev an optional integer that allows to specify the number of color levels used in legend strip.
#'
#' @return The empirical covariance, precision or decomposition matrix.
#' @export
#'
#' @importFrom fields image.plot
#' @importFrom grDevices hcl.colors
#' @importFrom graphics axis points title
#'
#' @examples
#' A <- matrix(c(1, 0.5, 0.25, 0.5, 2, 0.2, 0.25, 0.2, 3),3,3)
#' plotECPD(A, type = "dec", repar="mcd")
plotECPD <- function(ecov,
                     type = c("cov", "prec", "dec"),
                     repar = c("mcd", "logm"),
                     ord = NULL,
                     lab = NULL,
                     neimat  = NULL,
                     di = FALSE,
                     lev = NULL){
  d <- ncol(ecov)
  type <- match.arg(type)
  repar <- match.arg(repar)
  if(is.null(ord))  ord <- 1:d
  if(is.null(lab)) lab <- 1:d
  if(is.null(lev)) lev <- d * (d + 1)/2

  lab <- lab[ord]
  ecov <- ecov[ord, ord]
  res <- matrix(NA, d, d)

  if(type == "cov"){
    res[upper.tri(res, diag = di)] <- ecov[upper.tri(ecov, diag = di)]
    fields::image.plot(res, axes = FALSE,
                       col = hcl.colors(lev, "YlOrRd", rev = TRUE))
    axis(2, at = seq(0, 1, length = d), labels = lab)
    axis(3, at = seq(0, 1, length = d), labels = lab)
    title(sub = "Empirical Covariance Matrix", cex.sub = 1.5, line = 1)
  } else if (type == "prec") {
    pre <- solve(ecov)
    res[upper.tri(res, diag = di)] <- pre[upper.tri(pre, diag = di)]
    fields::image.plot(res, axes = FALSE,
                       col = hcl.colors(lev, "YlOrRd", rev = TRUE))
    axis(2, at = seq(0, 1, length = d), labels = lab)
    axis(3, at = seq(0, 1, length = d), labels = lab)
    title(sub = "Empirical Precision Matrix", cex.sub = 1.5, line = 1)
  } else {
    if(repar == "mcd"){
      chol_dec <- chol(ecov)
      D <- diag(chol_dec)
      L <- matrix(0, d, d)
      diag(L) <- 1
      for(j in 2:d){
        for(k in 1:(j - 1)){
          L[j, k] <- chol_dec[k, j]/D[k]
        }
      }

      res<-matrix(0, d, d)

      for(j in 1:d){
        res[j, j] = 1
        if(j > 1){
          for(k in 1:(j - 1)){
            s <- 0
            for(l in k:(j - 1)){
              s = s + L[j, l] * res[l, k]
            }
            res[j, k] = -s
          }
        }
      }
      res <- t(res)
      if(di == TRUE){
        diag(res) <- log(D^2)
      }
      res[lower.tri(res, diag = !di)] <- NA
      fields::image.plot(res, axes = FALSE,
                         col = hcl.colors(lev, "YlOrRd", rev = TRUE))
      axis(2, at = seq(0, 1, length = d), labels = lab)
      axis(3, at = seq(0, 1, length = d), labels = lab)
      if(di == TRUE){
        title(sub = "Theta matrix (logD2 diagonal, T outdiagonal)", cex.sub = 1.5, line = 1)
      } else {
        title(sub = "Theta matrix (T outdiagonal)", cex.sub = 1.5, line = 1)
      }
    } else {
      eig_dec <- eigen(ecov)
      res <- eig_dec$vectors%*%diag(log(eig_dec$values))%*%t(eig_dec$vectors)
      res[lower.tri(res, diag = !di)] <- NA
      fields::image.plot(res, axes = FALSE,
                         col = hcl.colors(lev, "YlOrRd", rev = TRUE))
      axis(2, at = seq(0, 1, length = d), labels = lab)
      axis(3, at = seq(0, 1, length = d), labels = lab)
      if(di == TRUE){
        title(sub = "Theta matrix (logSigma elements)", cex.sub = 1.5, line = 1)
      } else {
        title(sub = "Theta matrix (logSigma elements without diagonal)", cex.sub = 1.5, line=1)
      }
    }
  }
  if(!is.null(neimat)){
    neimat <- neimat[ord, ord]
    neighb <- seq(0, 1, length = d)
    for(i in 2:d){
      for(j in 1:(i - 1)){
        if(neimat[i, j] == TRUE) points(neighb[j], neighb[i], pch=8)
      }
    }
  }
  return(invisible(res))
}
