#' Processing the summary to be printed
#' @description  This function allows to modify the summary.gam function, in order to have a bespoke summary
#'
#' @param object an object of class scm
#' @param intercept set to TRUE if you want print the intercept
#' @param print print the summary
#' @param ... further arguments to be passed to gam
#'
#' @return The summary is printed
#' @export summary.scm
#' @export
#'
#' @importFrom stringr word
#' @importFrom mgcv summary.gam
#' @name summary.scm
#' @examples
summary.scm <- function(object, intercept = FALSE, print = TRUE, ...){

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

  sgam <- summary.gam(object)

  sgam.pt <- sgam$p.table
  lab_pt<-rownames(sgam.pt)
  for(i in 1:length(lab_pt)){
    detlp<-gregexpr("\\.", lab_pt[i])

    if(detlp[[1]][1] == -1){
      # Detect the lhs of the formula of the first terms
      lhs_foo <- word(as.character(object$foo_summary[[1]]), 1, sep = "\\~")[2]
      lab_pt[i] <- paste0(lab_pt[i],".", lhs_foo)
    } else {
      # Detect the lhs of the formula of the first terms
      lhs_foo <- word(as.character(object$foo_summary[[ as.integer(substring(lab_pt[i], first = detlp[[1]][length(detlp[[1]])]+1)) + 1]]), 1, sep = "\\~")[2]
      if(!is.na(lhs_foo)){
        lab_pt[i] <- paste0(substring(lab_pt[i], first=1, last = detlp[[1]][length(detlp[[1]])]), lhs_foo)
      } else {
        lab_pt[i] <- paste0("(Intercept).",internal()$labTh(object$family$getd(), as.integer(substring(lab_pt[i], first = detlp[[1]][length(detlp[[1]])]+1)) + 1))
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
      lhs_foo <- word(as.character(object$foo_summary[[1]]), 1, sep = "\\~")[2]
      lab_st[i] <- paste0(lab_st[i],".", lhs_foo)
    } else {
    # Detect the lhs of the formula of the first terms
    lhs_foo <- word(as.character(object$foo_summary[[ as.integer(substring(lab_st[i], first = detlp[[1]][length(detlp[[1]])]+1,last = detlpbl[[1]][length(detlpbl[[1]])]-1)) + 1]]), 1, sep = "\\~")[2]
    lab_st[i] <- paste0(substring(lab_st[i], first=1, last=1), substring(lab_st[i], first = detlpbl[[1]][length(detlp[[1]])], last=nchar(lab_st[i])),".",lhs_foo)
    }
  }
  rownames(sgam.st) <- lab_st

  if ( !intercept ) sgam.pt <- sgam.pt[-which(grepl( "(Intercept)", row.names(sgam.pt), fixed = TRUE)),] #(intercepts are removed)
  if(print){
    print.summary.gam(object, sgam.pt, sgam.st)
  }
  ret <- with(sgam, list(p.coeff = p.coeff, se = se, p.t = p.t, p.pv = p.pv,
                         residual.df = residual.df, m = m, chi.sq = chi.sq, s.pv = s.pv,
                         scale = dispersion, r.sq = r.sq, family = object$family,
                         formula = object$formula, n = nobs, dev.expl = dev.expl,
                         edf = edf, dispersion = dispersion, pTerms.pv = pTerms.pv,
                         pTerms.chi.sq = pTerms.chi.sq, pTerms.df = pTerms.df,
                         cov.unscaled = cov.unscaled, cov.scaled = cov.scaled,
                         p.table = sgam.pt, pTerms.table = pTerms.table, s.table = sgam.st,
                         method = object$method, sp.criterion = object$gcv.ubre,
                         rank = object$rank, np = length(object$coefficients)))
  class(ret) <- "summary.scm"
  return(invisible(ret))
}
