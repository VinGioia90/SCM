#' This function take the list of formulas provide by the user, then processing it  to be passed to the gam function to estimate gam_scm model
#'
#' @description These is a function that fits a GAM model via covariance modelling using [mgcv::gam]
#'
#' @param foo_user formula provided by the user
#' @param d dimension of the outcome
#' @return A list of formulas
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
  el_bar <- rep(0, length(foo_user))
  el_bar2 <- rep(0, length(foo_user))
  count<-1
  for(j in 1:length(foo_user)){
    idx_bar <- gregexpr("\\|", foo_user[[j]])[[2]]
    if(idx_bar[1] != -1){
      el_bar[j] <- j
      idx_bar <- c(0,idx_bar,(nchar( word(as.character(foo_user[[j]]), 1, sep = "\\~")[2]))+1 )
      for(i in 1:(length(idx_bar)-1)){
        foo_bar_lhs[[count]] <- substring(as.character(foo_user[[j]])[2], first = idx_bar[i]+1, last = idx_bar[i + 1] -1)
        foo_bar[[count]] <- paste(foo_bar_lhs[[count]], "~", as.character(rhs(foo_user[[j]]))[2])
        count <- count + 1
      }
    }
  }


  for(j in (1:length(foo_user))[-el_bar]){
    for(k in 1 : length(foo_bar)){
      if(lhs.vars(foo_user[[j]]) ==  lhs.vars(formula(foo_bar[[k]]))){
        foo_bar[[k]] <- paste(foo_bar[[k]],as.character( rhs(foo_user[[j]]))[2], sep="+")
        el_bar2[j] <- j
      }
    }
  }
  name_lhs <- rep(0, length(foo_bar))
  for(j in 1:length(foo_bar)){
    name_lhs[j]<-gsub("[[:space:]]", "",  word(as.character(foo_bar[[j]]), 1, sep = "\\~"))
  }

  count <- length(foo_bar) + 1
  for(j in (1:length(foo_user))[-c(el_bar,el_bar2)]){
    if( !(as.character(lhs.vars(foo_user[[j]])) %in%   name_lhs)) foo_bar[[count]] <- foo_user[[j]]
    count <- count + 1
  }


  foo_user2 <- lapply(foo_bar, formula, env = globalenv())

  foo_len <- length(foo_user2)

  if ( foo_len > (d + d*(d + 1)/2) ) stop("More formulas than elements to be modelled")
  #library(stringr)
  # Check on the mean model formula
  lhs <- unlist(lapply(1:foo_len, function(x) word(as.character(foo_user2[[x]]), 1, sep = "\\~")[2]))
  meanlhs <- lhs[which(!word(lhs, 1, sep = "\\Th_")=="")]
  if ( length(meanlhs) != d ) stop("More/Less mean components than d")
  if ( length(meanlhs) > length(unique(meanlhs)) ) stop("A mean component is modelled more than one time")


  # Build a list where we save:
  # [[1]]: the row and column of the components to be modelled (with a counter)
  # [[2]]: the list of formulas specified by the user
  nrc <- list()
  nrc[[1]] <- matrix(0, foo_len-d, 2)
  nrc[[2]] <- list()
  nrc[[3]] <- list()
  for( j in (d + 1) : foo_len ) {
    if( nchar(word(as.character(foo_user2[[j]]), 1, sep = "\\~")[2]) == (3 + 2 * d_len) ) {
      nrc[[1]][j - d, 1] <- as.integer(substring(format(foo_user2[[j]]), first = 4, last = 4 + d_len - 1))
      if( nrc[[1]][j - d, 1] == 0 | nrc[[1]][j - d, 1] > d) stop("Wrong specification of the formula: the elements 0 and the lements greater than d  could not be specified")
      nrc[[1]][j - d, 2] <- as.integer(substring(format(foo_user2[[j]]), first = 4 + d_len, last = 4 + 2 * d_len - 1))
      if( nrc[[1]][j - d, 2] == 0 | nrc[[1]][j - d, 2] > d) stop("Wrong specification of the formula: the elements 0 and the lements greater than d  could not be specified")
      if ( nrc[[1]][j - d, 1] >  nrc[[1]][j - d, 2] ) {
        nrc[[1]][j - d,] <- nrc[[1]][j - d, c(2:1)]
      }
      nrc[[2]][[j - d]] <- substring(format(foo_user2[[j]]), first = 4 + 2 * d_len)
      nrc[[3]][[j - d]] <- format(foo_user2[[j]])
    } else {
      stop("Wrong formula initialization, if the dimension of the outcome is < 10 specify two digits after Th_,
                                          if it is > 9 and < 100 specify four digits after Th_,
            see the vignettes")
    }
  }

  if ( nrow(unique(nrc[[1]][,1:2])) <  ( foo_len -d) ) stop ("Two elements are associated to different formulas")

  # Facciamno un controllo sulla specificazione del mean model
  mean_foo <- lapply(1:d, function(x) foo_user2[[x]])  # the first d components of the list lies on mean model specification
  cov_foo <- internal()$formula_init(d) # initialization of the covariance model
  cov_foop <- list() # initialization of the covariance model (to be printed)
  cov_foos <- list() # initialization of the covariance model (auxiliary)


  # Now, we build the index needed to the formula function to be passed to the gam function
  count <- d + 1
  idxCm <- matrix(0, d * (d + 1)/2, 3)
  for ( j in 1 : d ){
    idxCm[j, 1] <- idxCm[j, 2] <- idxCm[j,3] <- j
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
  for(i in 1:(d*(d+1)/2)){
    for( i2 in 1: nrow(nrc[[1]])){
      if ((idxCm[i, 1] %in% nrc[[1]][i2, 1]) & (idxCm[i, 2] %in% nrc[[1]][i2, 2]) ) {
        cov_foo[[idxCm[i, 3]]] <- formula(unlist(nrc[[2]][idxCm[i2, 3]]), env = globalenv())
        cov_foop[[count]] <- formula(unlist(nrc[[3]][idxCm[i2, 3]]), env = globalenv())
        cov_foos[[idxCm[i, 3]]] <- formula(unlist(nrc[[3]][idxCm[i2, 3]]), env = globalenv())
        count <- count + 1
      }
    }
  }

  return (list(foo_eval = c(mean_foo, cov_foo), foo_print = c(mean_foo,cov_foop), foo_summary=c(mean_foo, cov_foos)))

}



