#' Estimates GSVA enrichment scores.
#'
#'
#' @title pathway_gsva
#' @param data a matrix of expression values where rows correspond to genes and columns correspond to samples.
#' @param geneset a list or dataframe
#' @param gmt gmtfile
#' @param method one of "ssgsea","gsva", "zscore",and "plage"
#' @param rna_seq rna-seq data or microarray data
#' @importFrom immcp read_gmt
#' @importFrom immcp to_list
#' @importFrom GSVA gsva
#' @importFrom magrittr %>%
#' @return a data.frame
#' @export
#' @author Yuanlong Hu

pathway_gsva <- function(data,geneset,gmt=NULL, method=c("ssgsea","gsva", "zscore", "plage"),rna_seq=FALSE){

  if(is.null(gmt)){
    if (!is.list(geneset)) {
      geneset <- immcp::to_list(geneset)
    }
  }else{
    geneset <- immcp::read_gmt(gmt, out_type = "list")
  }
  res <- gsva(expr=as.matrix(data), gset.idx.list=geneset,
              method=method[1],
              kcdf=ifelse(rna_seq,"Poisson","Gaussian"),
              abs.ranking=FALSE,
              min.sz=1, max.sz=Inf,
              parallel.sz=1L,
              mx.diff=TRUE,
              ssgsea.norm=TRUE,verbose=TRUE)
  res <- as.data.frame(res)
  return(res)
}
