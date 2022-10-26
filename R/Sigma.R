Sigma_mat <- function(pred_Sigma){
  #!!!! rendere visibile???? Si costruire una funzione che prenda in input la matrice degli lpi e restituisce una lista
  # Such function takes the i-th predicted linear predictor element of the covariance modelling
  #(so excluding the mean vector and when predict is used in combination with type="response")
  # and return the covariance/correlation matrix (according to the flag)
  n_el <- length(pred_Sigma)
  d <- (-1 + sqrt(1 + 8 * n_el))/2
  res <- matrix(0, d, d)
  res[upper.tri(res, diag = FALSE)] <-  pred_Sigma[(d + 1) : n_el]
  res <- res + t(res)
  diag(res) <- pred_Sigma[1 : d]
  return(res)
}

