#' List of predicted Sigma matrix
#'
#' @description Such function returns a list of the predicted covariance matrix
#' @param pred_Sigma the matrix n x d(d+1)/2 of the estimated covariance matrix elements (from predict function with type="response" and excluding the mean vector)
#'
#' @return a list of n elements, each corresponding to the estimated d x d covariance matrix
#' @export
#'
Sigma_mat <- function(pred_Sigma){
  no_eta <- dim(as.matrix(pred_Sigma))[2]
  nobs <- dim(as.matrix(pred_Sigma))[1]
  #d <- -3/2 + sqrt(9/4 + 2 * no_eta)
  d <- (-1 + sqrt(1 + 8 * n_el))/2
  out <- list()
  for ( i in 1 : nobs){
    out[[i]] <-  matrix(0, d, d)
    out[[i]][upper.tri(out[[i]], diag = FALSE)] <- pred_Sigma[i, (d + 1) : no_eta]
    out[[i]] <- out[[i]] + t(out[[i]])
    diag(out[[i]]) <- pred_Sigma[i, 1 : d]
  }
  return (out)
}

