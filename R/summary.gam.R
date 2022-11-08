#' Processing the summary to be printed
#' @description  This function allows to modify the summary.gam function, in order to have a bespoke summary
#'
#' @param obj an object of class scm
#' @param intercept set to TRUE if you want print the inctercept
#'
#' @return The summary is printed
#' @export summary.scm
#' @export
#'
#' @importFrom stringr word
#' @importFrom mgcv summary.gam
#' @name summary.scm
#' @examples
summary.scm <- function(obj, intercept = FALSE){

  print.summary.gam <-function (x, parcoef, smooth, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), ...) {
    print(x$family)
    cat("Formula:\n")
    for (i in 1:length(x$foo_print)) if ( !(word(as.character(x$foo_print[[i]]), 1, sep = "\\~")[3]=="1")) print(x$foo_print[[i]])
    #for (i in 1:length(x$formula))  print(x$formula[[i]])
    #if (length(x$p.coeff) > 0) {
      cat("\nParametric coefficients:\n")
      printCoefmat(parcoef, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    #}
    cat("\n")
    if (length(x$smooth) > 0) {
      cat("Approximate significance of smooth terms:\n")
      printCoefmat(smooth, digits = digits, signif.stars = signif.stars, has.Pvalue = TRUE, na.print = "NA", cs.ind = 1, ...)
    }
    cat("\n")
    if (!is.null(x$rank) && x$rank < length(x$coefficients))
      cat("Rank: ", x$rank, "/", x$np, "\n",
          sep = "")
    if (!is.null(x$r.sq))
      cat("R-sq.(adj) = ", formatC(x$r.sq, digits = 3,
                                   width = 5), "  ")
    if (length(x$dev.expl) > 0)
      cat("Deviance explained = ", formatC(x$dev.expl *
                                             100, digits = 3, width = 4), "%", sep = "")
    cat("\n")
    if (!is.null(x$method) && !(x$method %in% c("PQL",
                                                "lme.ML", "lme.REML")))
      cat(x$method, " = ", formatC(x$gcv.ubre, digits = 5),
          sep = "")
    cat("  Scale est. = ", formatC(x$scale, digits = 5,
                                   width = 8, flag = "-"), "  n = ", nrow(x$model), "\n",
        sep = "")
    invisible(x)
  }

  sgam <- summary.gam(obj)

  sgam.pt <- sgam$p.table
  lab_pt<-rownames(sgam.pt)
  for(i in 1:length(lab_pt)){
    detlp<-gregexpr("\\.", lab_pt[i])

    if(detlp[[1]][1] == -1){
      # Detect the lhs of the formula of the first terms
      lhs_foo <- word(as.character(obj$foo_summary[[1]]), 1, sep = "\\~")[2]
      lab_pt[i] <- paste0(lab_pt[i],".", lhs_foo)
    } else {
      # Detect the lhs of the formula of the first terms
      lhs_foo <- word(as.character(obj$foo_summary[[ as.integer(substring(lab_pt[i], first = detlp[[1]][length(detlp[[1]])]+1)) + 1]]), 1, sep = "\\~")[2]
      if(!is.na(lhs_foo)){
        lab_pt[i] <- paste0(substring(lab_pt[i], first=1, last = detlp[[1]][length(detlp[[1]])]), lhs_foo)
      } else {
        lab_pt[i] <- paste0("(Intercept).",internal()$labTh(obj$family$getd(), as.integer(substring(lab_pt[i], first = detlp[[1]][length(detlp[[1]])]+1)) + 1))
      }
      }
  }
  rownames(sgam.pt) <- lab_pt

  sgam.st <- sgam$s.table
  lab_st<-rownames(sgam.st)
  for(i in 1:length(lab_st)){
    detlp<-gregexpr("\\.", lab_st[i])
    detlpbl<-gregexpr("\\(", lab_st[i])
    if(detlp[[1]][1] == -1){
      # Detect the lhs of the formula of the first terms
      lhs_foo <- word(as.character(obj$foo_summary[[1]]), 1, sep = "\\~")[2]
      lab_st[i] <- paste0(lab_st[i],".", lhs_foo)
    } else {
    # Detect the lhs of the formula of the first terms
    lhs_foo <- word(as.character(obj$foo_summary[[ as.integer(substring(lab_st[i], first = detlp[[1]][length(detlp[[1]])]+1,last = detlpbl[[1]][length(detlpbl[[1]])]-1)) + 1]]), 1, sep = "\\~")[2]
    lab_st[i] <- paste0(substring(lab_st[i], first=1, last=1), substring(lab_st[i], first = detlpbl[[1]][length(detlp[[1]])], last=nchar(lab_st[i])),".",lhs_foo)
    }
  }
  rownames(sgam.st) <- lab_st

  if ( !intercept ) sgam.pt <- sgam.pt[-which(grepl( "(Intercept)", row.names(sgam.pt), fixed = TRUE)),] #(intercepts are removed)
  class(sgam) <- "summary.scm"
  return(print.summary.gam(obj, sgam.pt, sgam.st))
}
