#' Fit single lm methods
#'
#'
#' @title select_lm_single
#' @param expr a matrix of expression values where rows correspond to genes and columns correspond to samples.
#' @param pdata a vector.
#' @importFrom pbapply pblapply
#' @return a data.frame
#' @export
#' @author Yuanlong Hu

select_lm_single <- function(expr, pdata){

  expr <- as.data.frame(t(expr))
  if(!is.vector(pdata)) stop("The pdata must be a vector.")
  f <- paste0("y","~","`",colnames(expr),"`")

  expr$y <- pdata

  res_list <- pblapply(f, function(x){
    x <- as.formula(x)
    fit <- lm(x,expr)
    fit_s <- summary(fit)
    df <- data.frame(y=as.character(fit_s$terms)[3],
                     Beta=fit_s$coefficients[2,1],
                     Pvalue=fit_s$coefficients[2,4],
                     CI1=confint(fit)[2,1],
                     CI2=confint(fit)[2,2]
    )
    return(df)
  })
  res_list <- Reduce(rbind, res_list)
  return(res_list)
}

#' Select by Boruta
#'
#'
#' @title select_Boruta
#' @param expr a matrix of expression values where rows correspond to genes and columns correspond to samples.
#' @param pdata a vector.
#' @importFrom Boruta Boruta
#' @return a Boruta object
#' @export
#' @author Yuanlong Hu

select_Boruta <- function(expr, pdata){
  expr <- as.data.frame(t(expr))
  expr$pdata <- pdata
  res <- Boruta::Boruta(pdata~., expr, doTrace=1)
  return(res)
}


#' plot Boruta ImpHistory
#'
#'
#' @title plotBorutaImpHistory
#' @param res the result of Boruta
#' @param select a vector.
#' @importFrom Boruta Boruta
#' @return a Boruta object
#' @export
#' @author Yuanlong Hu

plotBorutaImpHistory <- function(res,
                                 select = NULL,
                                 select2 = c("Confirmed", "Rejected", "Tentative")){
  if (is.null(select)) {
    select <- names(res$finalDecision)
  }

  ImpHistory <- res$ImpHistory[,select]
  ImpHistory <- as.data.frame(t(ImpHistory))
  ImpHistory$Group <- res$finalDecision[select]
  ImpHistory$gene <- rownames(ImpHistory)
  ImpHistory <- reshape2::melt(ImpHistory, id.vars=c("gene","Group"))
  ImpHistory <- ImpHistory[ImpHistory$Group %in% select2,]
  ImpHistory$gene <- with(ImpHistory, reorder(gene, value, median))

  p <- ggplot(ImpHistory, aes(x=value, y=gene, fill=Group))+
    geom_boxplot()+
    scale_fill_jco()+
    theme_minimal()+
    labs(x="Importance",y="Attributes",fill="")+
    theme(legend.justification=c(1,0), legend.position=c(1,0))
  return(p)
}


#' plot ExprBox
#'
#'
#' @title plotExprBox
#' @param expr the result of Boruta
#' @param select a vector.
#' @param pdata pdata
#' @param comparisons a list
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_violin
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggsci scale_fill_jco
#' @importFrom ggpubr stat_compare_means
#' @return a Boruta object
#' @export
#' @author Yuanlong Hu

plotExprBox <- function(expr, select, pdata, comparisons=list(c("C1","C2"))){

  expr <- data.frame(object=as.numeric(expr[select,]),
                     group=pdata)
  p <- ggplot(expr, aes(x=group, y=object, fill=group))+
    geom_violin()+
    geom_boxplot(width=0.2, fill="white")+
    scale_fill_jco()+
    stat_compare_means(comparisons = comparisons)+
    theme_minimal()
  return(p)

}
