##' Calculate differentially expressed genes using limma method.
##'
##'
##' @title findDEG
##' @param data A matrix of expression values where rows correspond to genes and columns correspond to samples.
##' @param pdata A character vector of phenotype.
##' @param contrasts list of specifying contrasts
##' @return A list
##' @importFrom stats model.matrix
##' @importFrom limma makeContrasts
##' @importFrom limma lmFit
##' @importFrom limma contrasts.fit
##' @importFrom limma eBayes
##' @importFrom limma topTable
##' @importFrom DESeq2 DESeqDataSetFromMatrix
##' @importFrom DESeq2 DESeq
##' @importFrom DESeq2 results
##' @importFrom magrittr %>%
##' @export
##' @author Yuanlong Hu

findDEG <- function(data, pdata, contrasts, rnaseq=F){

  if(rnaseq){
    coldata <- data.frame(sample=colnames(data),
                          condition=pdata)
    rownames(coldata) <- coldata$sample

    dds <- DESeqDataSetFromMatrix(countData = data,
                                  colData = coldata,
                                  design = ~condition)
    dds <- DESeq(dds)

    res_DEG0 <- lapply(contrasts, function(x){
      res <- results(dds, contrast=c("condition",x[1],x[2]))
      res <- res[order(res$log2FoldChange, decreasing = T),]
      colnames(res) <- c("baseMean", "logFC", "lfcSE", "stat",
                         "P.Value", "adj.P.Val")
      return(res)
    })

    n <- lapply(contrasts, function(x){
      paste0(x[1], "_vs_",x[2])
    }) %>%
      unlist()
    names(res_DEG0) <- n

  }else{

  Group <- factor(pdata)
  design <- model.matrix(~0 + Group)
  colnames(design) <- gsub("Group","", colnames(design))

  # Construct Matrix of Custom Contrasts

  if(is.list(contrasts)){
    contrasts <- lapply(contrasts, function(x){
    paste0(x[1], "-",x[2])
  }) %>%
    unlist()
  }

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

  }

  return(res_DEG0)
}

##' Heatmap showing the summary of DEGs
##'
##'
##' @title plotSummaryDEGHeatmap
##' @param reslist the list of results
##' @param select selected genes
##' @param pvalueCutoff pvalue Cutoff
##' @param adjustpCutoff adjust p Cutoff
##' @param marker P marker
##' @return ggplot2
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_tile
##' @importFrom ggplot2 scale_fill_gradient2
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 geom_text
##' @importFrom scales muted
##' @export
##' @author Yuanlong Hu

plotSummaryDEGHeatmap <- function(reslist, select,
                                  pvalueCutoff=1,
                                  adjustpCutoff=1,
                                  marker=F){
  reslist <- lapply(reslist, function(x){
    x <- na.omit(x[select,])
    x$SYMBOL <- rownames(x)
    return(x)
  })

  for (i in 1:length(reslist)) {
    reslist[[i]]$group <- names(reslist)[i]
  }

  reslist <- Reduce(rbind, reslist)

  reslist <- reslist[reslist$P.Value<pvalueCutoff & reslist$adj.P.Val<adjustpCutoff,]
  reslist$SYMBOL <- with(reslist, reorder(SYMBOL, logFC, mean))
  reslist$marker <- ifelse(reslist$adj.P.Val<0.001,"***",
                           ifelse(reslist$adj.P.Val<0.01,"**",
                                  ifelse(reslist$adj.P.Val<0.05,"*",NA)
                           ))

  p <- ggplot(reslist, aes(x=SYMBOL, y=group))+
    geom_tile(aes(fill=logFC),color="gray",alpha=0.8)+
    scale_fill_gradient2(low = scales::muted("blue"),
                         high = scales::muted("red"))+

    theme_minimal()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.x = element_text(angle=45,hjust=1, vjust=1))

  if(marker){
    p <- p+geom_text(aes(label=marker))
  }

  return(p)
}
