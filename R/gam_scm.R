#' Fit a multivariate GAM model with covariance modelling
#'
#' @description These is a function that fits a GAM model via covariance modelling using [mgcv::gam]
#'
#' @param formula formula provided by the user
#' @param family set the family to mvn_scm (whose parameters are d - dimension of the outcome; nb - number of observations' blocks; param - the type of parametrisation: cuurently only mcd (default) and logm are available)
#' @param data  same arguments as in [mgcv::gam]
#' @param auxGam list of further arguments to be passed to [mgcv::gam]
#' @return An object of class "scm"
#' @name gam_scm
#' @export
#'
#' @importFrom Rcpp evalCpp
#' @importFrom mgcv gam summary.gam print.summary.gam
#'
#' @examples
gam_scm <- function(formula, family = mvn_scm(d = 2, nb = 1, param = NULL), optimizer = NULL, data = list(), auxGam = list()){ #passare nb come auxSCM
  d <- family$getd()
  q <- d + d * (d + 1)/2

  param <- family$getparam()
  if (is.null(optimizer) | param == "logm")    opt <- "efs" else opt <- optimizer

  # here we convert the formula provided by the user to a formula that could be used to the family
  # Check on the model formula provide by the users (especially check on the number of components specified and if there are variables in the dataset called like Th_<something>)
  idxTh <- which(grepl( "Th_", colnames(data), fixed = TRUE))
  if ( !(identical(idxTh,integer(0))) ) {
    Th_v<- rep(0, length(idxTh))
    for(j in 1 : length(idxTh)) Th_v[j] <- gregexpr("Th_", colnames(data[idxTh[j]]))[[1]][1]
    if(1 %in% Th_v) stop("One (or more) variable in the data frame is specified as Th_<something>: Please change the label of such variable")
  }

  # !!! In the new formulation Iwe need to use q+1 to deal with the case where there are formulae specified by using |; We need to avoid the possibility to specify two formulas for the same element as
  # Th_11 ~ dow
  # Th_11 ~ s(doy)
  #if(length(formula) > q)  stop("Incorrect specifiction of the model formula: more formula than elements to be modelled")

  # Processing the model formula to be used by the family mvn_scm
  foo <- get_foo(formula, d = d)

  # Fit the model
  obj <- do.call("gam", c(list("formula" = foo$foo_eval, "family" = family,
                              "data"= data, "optimizer" = opt), auxGam))
  obj$foo_print <- foo$foo_print
  obj$foo_summary <- foo$foo_summary
  class(obj) <- c("scm", class(obj))

  return(obj)
}
