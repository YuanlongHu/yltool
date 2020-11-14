#' Calculate the similarity
#'
#'
#' @title sim_expr
#' @param data data
#' @return a data.frame
#' @importFrom Hmisc rcorr
#' @importFrom magrittr %>%
#' @export
#' @author Yuanlong Hu



sim_expr <- function(data, method=c("pearson","spearman")){

  cormat <- rcorr(data, type=method[1])
  cormat <- convCorrMatrix(cormat$r,cormat$P)
  return(cormat)
}




#' Conversion CorrMatrix
#'
#'
#' @title convCorrMatrix
#' @param mat CorrMatrix
#' @param pmat p-value matrix
#' @return a data.frame
#' @author Yuanlong Hu
#' @noRd
convCorrMatrix <- function(mat, pmat) {
  ut <- upper.tri(mat)
  data.frame(
    row = rownames(mat)[row(mat)[ut]],
    column = rownames(mat)[col(mat)[ut]],
    cor  =(mat)[ut],
    p = pmat[ut]
  )
}
