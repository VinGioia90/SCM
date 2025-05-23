#' Fit a multivariate GAM model with covariance modelling
#'
#' @description This is a function that fits a GAM model via covariance modelling using [mgcv::gam]
#'
#' @param formula formula provided by the user
#' @param family set the family to mvn_scm (whose parameters are d - dimension of the outcome; nb - number of observations' blocks; param - the type of parametrisation: cuurently only mcd (default) and logm are available)
#' @param optimizer set to efs or bfgs (avialable for the mcd)   to optimize the smoothing parameters
#' @param data  same arguments as in [mgcv::gam]
#' @param aGam list of further arguments to be passed to [mgcv::gam]
#' @return An object of class "scm"
#' @name gam_scm
#' @export
gam_scm <- function(formula, family = mvn_scm(d = 2, nb = 1, param = NULL), optimizer = NULL, data = list(), aGam = list()){
  
  if( is.null(aGam$G) ){
    d <- family$getd()
    q <- d + d * (d + 1)/2 #number of lpi
    
    param <- family$getparam()
    if (is.null(optimizer) || param == "logm") opt <- "efs" else opt <- optimizer
    
    # Check on the model formula provided by the users: check if there are variables in the dataset called  Th_<something>
    idxTh <- which(grepl( "Th_", colnames(data), fixed = TRUE))
    if ( !(identical(idxTh, integer(0))) ) {
      Th_v <- rep(0, length(idxTh))
      for(j in 1 : length(idxTh)) Th_v[j] <- gregexpr("Th_", colnames(data[idxTh[j]]))[[1]][1]
      if(1 %in% Th_v) stop("Some variables in the data frame are specified as Th_<something>: Please change the label of such variables")
    }
    
    # Processing the model formula to the formula  used by the family mvn_scm
    foo <- get_foo(formula, d = d)
    
    call_list <- list("formula" = foo$foo_eval, "family" = family,
                      "data" = data, "optimizer" = opt)
  } else {
    
    param <- aGam$G$family$getparam()
    if (is.null(optimizer) || param == "logm") opt <- "efs" else opt <- optimizer
    call_list <- list(optimizer = opt)
    
  }
  
  # Fit the model
  obj <- do.call("gam", c(call_list, aGam))
  
  if( is.null(aGam$G) ){
    # foo_print and foo_summary are the formulas by the summary function
    obj$foo_print <- foo$foo_print
    obj$foo_summary <- foo$foo_summary
  } else {
    obj$foo_print <- aGam$G$foo_print
    obj$foo_summary <- aGam$G$foo_summary
  }
  
  class(obj) <- c("scm", class(obj))
  return(obj)
}
