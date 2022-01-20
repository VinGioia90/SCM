#' Extract information from the boosting procedure
#' @description This function allows to obtain a summary of the effects chosen by the boosting
#'
#' @param data_boost an object obtained by ...
#'
#' @return The function
#' @export
#' @importFrom stringr word
#'
#' @examples
#'
boost_eff <- function(data_boost){

  l_mat <- function(d){
    mat <- matrix(NA, d, d)
    for(i in 1:d){
      for(j in 1:i){
        if(i == j){
          mat[i, i] <- paste0("logD2[", i, ",", i, "]")
        } else {
          mat[i, j] <- mat[j, i] <- paste0("T[", i, ",", j, "]")
        }
      }
    }
    return(mat = mat)
  }

  A <- l_mat(data_boost$d)

  step_upd <- data_boost$nstep
  max_imp <- rep(0,step_upd)
  elem <- rep(0,step_upd)
  eff <- rep(0,step_upd)
  row <- rep(0, step_upd)

  for(i in 1: step_upd){
    id_boost <- data_boost$idx[[i]]
    max_imp[i] <- data_boost$dll[[i]][id_boost]
    elem[i] <- c(diag(A), A[upper.tri(A)])[id_boost[1,1]]
    row[i] <- id_boost[1,1]
    eff[i] <- id_boost[1,2]
  }
  elab<-paste(elem, "-", eff)
  erow <- paste(row, "-", eff)

  uni_elab <- unique(elab)

  s_imp <- rep(0, length(uni_elab))
  for(i in 1:length(uni_elab)){
    s_imp[i] <- sum(max_imp[which(elab %in%  uni_elab[i])])
  }

  dsum_imp <-data.frame(sum_imp = s_imp,
                        uni_label =  uni_elab,
                        ele_rows = as.numeric(word(unique(erow),1)) ,
                        label  =  as.numeric(word(unique(elab),3)),
                        name_eff = data_boost$effects[as.numeric(word(uni_elab,3))])

  dsum_imp <- dsum_imp[order(dsum_imp$sum_imp, decreasing="TRUE"),]
  return(data_sum_imp = dsum_imp)
}


