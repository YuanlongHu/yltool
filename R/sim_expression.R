#' Calculate the similarity
#'
#'
#' @title sim_expr
#' @param data data
#' @importFrom Hmisc rcorr
#' @importFrom magrittr %>%
#' @return a data.frame
#' @export
#' @author Yuanlong Hu



sim_expr <- function(data, gene, method=c("pearson","spearman"), pcutoff=0.05, cut=2){

  data <- as.data.frame(data)
  seed <- as.numeric(data[gene,])

  if (cut > (nrow(data)/2)-1) {
    stop("cut > (nrow(data)/2)-1)")
  }

  data$spilt <- cut(1:nrow(data), cut, labels=F)
  data <- split(data,data$spilt)

  res <- lapply(data, function(x){
    x <- x[,-ncol(x)]
    x <- as.data.frame(t(x))
    x$seed <- seed
    res <- rcorr(as.matrix(x), type=method[1])
    res <- convCorrMatrix(res$r,res$P)
    res <- res[res$row=="seed",]
    res <- res[,-1]
    colnames(res) <- c("Gene","Cor","Pvalue")
    res <- res[res$Pvalue < pcutoff,]
    return(res)
  })

  res <- Reduce(rbind,res)
  return(res)
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
