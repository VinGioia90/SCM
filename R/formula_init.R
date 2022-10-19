formula_init <- function(d) {
  Theta_formula_int <-  list()
  for(ii in 1 : (d * (d + 1)/2)){
    Theta_formula_int[[ii]] <- as.formula("~ 1", env = globalenv())
  }
  return(Theta_formula_int)
}
