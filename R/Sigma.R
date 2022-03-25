#' Return the covariance matrix in matrix
#' @description The function returns the covariance matrix
#' @param pred_Sigma A vector with elements needed to fill the covariance matrix
#'
#' @return The function returns the covariance matrix in matrix form
#' @export
#'
#' @examples
Sigma_mat <- function(pred_Sigma){
  n_el <- length(pred_Sigma)
  d <- (-1 + sqrt(1 + 8 * n_el))/2
  res <- matrix(0, d, d)
  res[upper.tri(res, diag = FALSE)] <-  pred_Sigma[(d+1):n_el]
  res <- res + t(res)
  diag(res) <- pred_Sigma[1:d]
  return(res)
}

