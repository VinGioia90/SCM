gamlss.gH_scm <- function (X, jj, eta, y, w, z, t, Gm, # l3 = 0,
                           l1, l1_l, l2, l2_v, l2_l, l2_v_l, idx_b, idx_aux, param = 1, #by default param=1 for the mcd parameterisation
                           l3, l3_l, idx_l3, idx_neq0, idx_jkq,
                           a, a_l, d1H, # V, V_l,
                           l4 = 0, i4 = 0, d1b = 0, d2b = 0, deriv = 0, fh = NULL, D = NULL){
  K <- length(jj) # no_eta
  sparse <- discrete <- FALSE
  p <- ncol(X)
  n <- nrow(X)

  trHid2H <- d2H <- NULL
  if( deriv < 1 ) d1H <- NULL


  lpi <- lapply(jj, function(x) x - 1) # trick to allow the use of the c++functions
  #####################
  # Score w.r.t. beta #
  #####################
  lb <- rep(0, p)
  internal()$d1_beta(X, eta, y, lpi, K, lb, l1, l1_l, idx_b, z, w, Gm, param) #first derivatives w.r.t beta

  #######################
  # Hessian w.r.t. beta #
  #######################
  lbb <- if ( sparse )
    Matrix(0, p, p)
  else matrix(0, p, p)

  internal()$d2_beta(X, eta, y, lpi, K, lbb, l2, l2_v, l2_l, l2_v_l, idx_b, z, w, Gm, t, #second derivatives w.r.t beta
                     idx_aux$b1_eta, idx_aux$b1, idx_aux$b2, idx_aux$b3, idx_aux$idx_b1_eta, idx_aux$idx_b3,
                     idx_aux$l2_el, idx_aux$l2_el2 , param)

  #######################################
  # Derivative of the Hessian w.r.t eta #
  #######################################
  if ( deriv > 0 & param == 1) {
      if (discrete) stop("er... no discrete methods for higher derivatives")
      internal()$dHess_drho(X, eta, y, lapply(jj, function(x) x - 1), K, l3, l3_l, idx_b, d1b, #d1eta, d1eta_l, #V, V_l,
                                a, a_l, d1H, fh,  z, w, Gm, t,  idx_l3, idx_neq0, idx_jkq)
    }# end  if deriv > 0 & param == 1

  if ( deriv > 0 & param == 2) {
    warning("bfgs is not available for the logm parametrisation: the fitting is done by using efs optimizer")
    d1H <- NULL
  }


  list(lb = lb, lbb = lbb, d1H = d1H, d2H = d2H, trHid2H = trHid2H)
}
