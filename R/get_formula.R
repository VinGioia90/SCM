#' Processing the model formula provided by the users
#'
#' @description This function take the list of formulas provided by the user, then they are  processed to be passed to the gam function and to be estimated by gam_scm model
#'
#' @param foo_user formula provided by the user
#' @param d dimension of the outcome
#' @return Three list of formulas:one need to be passed to the gam function and the other two are useful for the summary function
#' @export
#'
#' @importFrom stringr word
#' @importFrom BMisc rhs lhs.vars
#' @examples
get_foo <- function ( foo_user, d) {
  d_len <- nchar(d)
  foo_len <- length(foo_user)
  foo_bar_lhs <- list()
  foo_bar <-list()
  el_bar <- rep(0, foo_len)
  el_bar2 <- rep(0, foo_len)

  # It is possible to specify a formula for an element no more than on time in the common model formula statement
  # and no more than one time in the term specific model formula

  # Here we detect the common model formula specification:
  count<-1
  for(j in 1 : foo_len) {
    idx_bar <- gregexpr("\\|", foo_user[[j]])[[2]]
    if ( idx_bar[1] != -1 ) {
      el_bar[j] <- j
      idx_bar <- c(0, idx_bar, (nchar( word(as.character(foo_user[[j]]), 1, sep = "\\~")[2])) + 1)
      for(i in 1: (length(idx_bar) - 1)){
        foo_bar_lhs[[count]] <- substring(as.character(foo_user[[j]])[2], first = idx_bar[i] + 1, last = idx_bar[i + 1] -1)
        foo_bar[[count]] <- paste(foo_bar_lhs[[count]], "~", as.character(rhs(foo_user[[j]]))[2])
        count <- count + 1
      }
    }
  }

  if ( length(unique(unlist(foo_bar_lhs))) <  length(unlist(foo_bar_lhs)) ) stop ("It is not possible to specify two formulas for the same element in common formula statement specification")
  nsp_foonames <- lapply((1 : foo_len)[-el_bar], function(x) lhs.vars(foo_user[[x]]))
  if ( length(unique(unlist(nsp_foonames))) <  length(unlist(nsp_foonames)) ) stop ("It is not possible to specify two formulas for the same element in term specific formula statement specification")

  foob_len <- length(foo_bar)

  for(j in (1 : foo_len)[-el_bar]){
    for(k in 1 : foob_len) {
      if ( lhs.vars(foo_user[[j]]) ==  lhs.vars(formula(foo_bar[[k]])) ) {
        foo_bar[[k]] <- paste(foo_bar[[k]], as.character(rhs(foo_user[[j]]))[2], sep = "+")
        el_bar2[j] <- j
      }
    }
  }

  name_lhs <- rep(0, foob_len)
  for(j in 1 : foob_len) name_lhs[j] <- gsub("[[:space:]]", "",  word(as.character(foo_bar[[j]]), 1, sep = "\\~"))

  count <- foob_len + 1
  for(j in (1 : foo_len)[-c(el_bar, el_bar2)]){
    if ( !(as.character(lhs.vars(foo_user[[j]])) %in%   name_lhs) ) foo_bar[[count]] <- foo_user[[j]]
    count <- count + 1
  }


  foo_user2 <- lapply(foo_bar, formula, env = globalenv())

  foo_len <- length(foo_user2)

  #if ( foo_len > (d + d*(d + 1)/2) ) stop("More formulas than elements to be modelled")
  #library(stringr)
  # Check on the mean model formula
  lhs <- unlist(lapply(1 : foo_len, function(x) word(as.character(foo_user2[[x]]), 1, sep = "\\~")[2]))
  meanlhs <- lhs[which(!word(lhs, 1, sep = "\\Th_")=="")]
  if ( length(meanlhs) != d ) stop("More/Less mean components than d")
  if ( length(meanlhs) > length(unique(meanlhs)) ) stop("A mean component is modelled more than one time")


  # Build a list where we save:
  # [[1]]: the row and column of the components to be modelled (with a counter)
  # [[2]]: the list of formulas specified by the user
  # [[3]]:

  nrc <- list()
  nrc[[1]] <- matrix(0, foo_len - d, 2)
  nrc[[2]] <- list()
  nrc[[3]] <- list()
  for ( j in (d + 1) : foo_len ) {
    if ( nchar(word(as.character(foo_user2[[j]]), 1, sep = "\\~")[2]) == (3 + 2 * d_len) ) {
      nrc[[1]][j - d, 1] <- as.integer(substring(format(foo_user2[[j]]), first = 4, last = 4 + d_len - 1))
      if ( nrc[[1]][j - d, 1] == 0 | nrc[[1]][j - d, 1] > d ) stop("Wrong specification of the formula: the elements 0 and the lements greater than d  could not be specified")
      nrc[[1]][j - d, 2] <- as.integer(substring(format(foo_user2[[j]]), first = 4 + d_len, last = 4 + 2 * d_len - 1))
      if ( nrc[[1]][j - d, 2] == 0 | nrc[[1]][j - d, 2] > d ) stop("Wrong specification of the formula: the elements 0 and the lements greater than d  could not be specified")
      if ( nrc[[1]][j - d, 1] >  nrc[[1]][j - d, 2] ) {
        nrc[[1]][j - d,] <- nrc[[1]][j - d, c(2:1)]
      }
      nrc[[2]][[j - d]] <- substring(format(foo_user2[[j]]), first = 4 + 2 * d_len)
      nrc[[3]][[j - d]] <- format(foo_user2[[j]])
    } else {
      stop("Wrong formula initialization, if the dimension of the outcome is < 10 specify two digits after Th_,
                                          if it is > 9 and < 100 specify four digits after Th_, ...")
    }
  }

  if ( nrow(unique(nrc[[1]][,1:2])) <  ( foo_len -d) ) stop ("Two elements are associated to different formulas")

  # Facciamno un controllo sulla specificazione del mean model
  mean_foo <- lapply(1 : d, function(x) foo_user2[[x]])  # the first d components of the list lies on mean model specification
  cov_foo <- internal()$formula_init(d) # initialization of the covariance model
  cov_foop <- list() # initialization of the covariance model (to be printed)
  cov_foos <- list() # initialization of the covariance model (auxiliary for the summary)


  # Now, we build the index needed to the formula function to be passed to the gam function
  count <- d + 1
  idxCm <- matrix(0, d * (d + 1)/2, 3)
  for ( j in 1 : d ){
    idxCm[j, 1] <- idxCm[j, 2] <- idxCm[j, 3] <- j
    if ( j > 1 ) {
      for ( k in 1 : (j - 1) ) {
        idxCm[count, 1] <- k
        idxCm[count, 2] <- j
        idxCm[count, 3] <- count
        count <- count + 1
      }
    }
  }

  count <- 1
  for(i in 1 : (d * (d + 1)/2)){
    for( i2 in 1: nrow(nrc[[1]])){
      if ((idxCm[i, 1] %in% nrc[[1]][i2, 1]) & (idxCm[i, 2] %in% nrc[[1]][i2, 2]) ) {
        cov_foo[[idxCm[i, 3]]] <- formula(unlist(nrc[[2]][idxCm[i2, 3]]), env = globalenv())
        cov_foop[[count]] <- formula(unlist(nrc[[3]][idxCm[i2, 3]]), env = globalenv())
        cov_foos[[idxCm[i, 3]]] <- formula(unlist(nrc[[3]][idxCm[i2, 3]]), env = globalenv())
        count <- count + 1
      }
    }
  }
  foo_eval <- c(mean_foo, cov_foo)
  if(length(foo_eval) > (d + d*(d + 1)/2)) stop("Some problems in model formula")
  foo_print <- c(mean_foo,cov_foop)
  foo_summary <- c(mean_foo, cov_foos)

  return (list(foo_eval = foo_eval, foo_print = foo_print, foo_summary = foo_summary))

}



