##' Calculate differentially expressed genes using limma method.
##'
##'
##' @title findDEG
##' @param data A matrix of expression values where rows correspond to genes and columns correspond to samples.
##' @param pdata A character vector of phenotype.
##' @param contrasts character vector specifying contrasts
##' @return A list
##' @importFrom stats model.matrix
##' @importFrom limma makeContrasts
##' @importFrom limma lmFit
##' @importFrom limma contrasts.fit
##' @importFrom limma eBayes
##' @importFrom limma topTable
##' @importFrom magrittr %>%
##' @export
##' @author Yuanlong Hu

findDEG <- function(data, pdata, contrasts){
  Group <- factor(pdata)
  design <- model.matrix(~0 + Group)
  colnames(design) <- gsub("Group","", colnames(design))

  # Construct Matrix of Custom Contrasts
  contrast_matrix <- makeContrasts(contrasts = contrasts,
                                   levels = design)
  fit <- lmFit(data, design) %>%
    contrasts.fit(contrast_matrix) %>%
    eBayes()

  res_DEG0 <- NULL
  for (i in 1:length(contrasts)) {
    res_DEG <-  topTable(fit, adjust.method = "fdr", number = Inf, coef = i) %>%
      list()
    names(res_DEG) <- contrasts[i]
    res_DEG0 <- c(res_DEG0, res_DEG)
  }
  return(res_DEG0)
}


