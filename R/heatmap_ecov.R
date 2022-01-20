#' Heatmap of the Empirical Covariance transformation.
#' @description plotECPD allows to plot the heatmap of the covariance matrix transformation.
#'
#' @param ecov empirical covariance matrix
#' @param trans an optional function that allows to obtain a transformation of the empirical covariance matrix
#' @param lab an optional vector that allows to specify the label of the variables.
#' @param di logical. If FALSE (default) the diagonal elements are not showed.
#' @param ... optional arguments to pass to hcl.colors of the image.plot function
#'
#' @return The empirical covariance, precision or decomposition matrix.
#' @export
#'
#' @importFrom fields image.plot
#' @importFrom grDevices hcl.colors
#' @importFrom graphics axis points title
#'
#' @examples
#' A <- matrix(c(1, 0.5, 0.25, 0.5, 2, 0.1, 0.25, 0.1, 3),3,3)
#' plotECPD(A)
#' plotECPD(A, trans = inverse)
#' plotECPD(A, trans = mcd)
#' plotECPD(A, trans = logm)
plotECPD <- function(ecov,
                      trans = NULL,
                      lab = NULL,
                      di = FALSE,
                      ...){
  d <- ncol(ecov)
  if(is.null(lab)) lab <- 1:d
  if(is.null(trans)){
    res <- ecov
  } else {
    res <- trans(ecov)
  }
  arg <- list(...)
  if(is.null(arg$col) ){
    arg$col <- hcl.colors(ifelse(di == TRUE, d * (d + 1)/2, d * (d - 1)/2),
                          "YlOrRd", rev = TRUE)
  }

  res[lower.tri(res, diag = !di)] <- NA
  do.call("image.plot", c(list(res, "axes" = FALSE), arg))
  axis(2, at = seq(0, 1, length = d), labels = lab)
  axis(3, at = seq(0, 1, length = d), labels = lab)

  return(invisible(res))
}

