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

#' Fit single cox methods
#'
#'
#' @title select_cox_single
#' @param expr a matrix of expression values where rows correspond to genes and columns correspond to samples.
#' @param time a vector.
#' @param status status
#' @param digits digits
#' @importFrom pbapply pblapply
#' @return a data.frame
#' @export
#' @author Yuanlong Hu

select_cox_single <- function(expr, time, status, digits=2){

  formulas <- sapply(rownames(expr),
                          function(x) as.formula(paste0('Surv(time, status)~', "`",x,"`")))

  expr <- as.data.frame(t(expr))
  expr$time <- as.numeric(time)
  expr$status <- as.numeric(status)
  message("*** Fitting Cox Models ***")
  models <- pblapply(formulas, function(x){coxph(x, data = expr)})

  message("*** Summary Cox Models Results ***")
  res <- pblapply(models,function(x){
                   x <- summary(x)
                   p.value<-signif(x$wald["pvalue"], digits=digits)
                   wald.test<-signif(x$wald["test"], digits=digits)
                   beta <- signif(x$coef[1], digits=digits);#coeficient beta
                  HR <- signif(x$coef[2], digits=digits);#exp(beta)
                  HR.confint.lower <- signif(x$conf.int[,"lower .95"],digits)
                  HR.confint.upper <- signif(x$conf.int[,"upper .95"],digits)
                  # HR <- paste0(HR, " (",
                  #                       HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR, HR.confint.lower, HR.confint.upper,wald.test, p.value)
                           names(res)<-c("beta", "HR", "CI_lower", "CI_upper", "wald.test", "p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                         })
  res <- t(as.data.frame(res, check.names = FALSE))
  res <- as.data.frame(res)
  return(res)
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
#' @param expr the expr data
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
#' @return a ggplot2 object
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

#' plot ExprVolcano
#'
#'
#' @title plotExprVolcano
#' @param res the result of limma
#' @param selectlabels a vector.
#' @param logFCcutoff logFC cutoff
#' @param xlab a charact
#' @param ylab a charact
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggrepel geom_text_repel
#' @return a ggplot2 object
#' @export
#' @author Yuanlong Hu


plotExprVolcano <- function(res, selectlabels=NULL,logFCcutoff=1,
                 xlab="log FoldChange", ylab="-log10 adjusted P-value"){

  res$SYMBOL <- rownames(res)
  res$group <- ifelse(res$adj.P.Val >=0.05|abs(res$logFC)< logFCcutoff,"none",
                      ifelse(res$logFC>=logFCcutoff,"up","down"))
  res$labels <- ifelse(res$SYMBOL %in% selectlabels, res$SYMBOL, NA)

  p <- ggplot(res, aes(x=logFC, y=-log10(adj.P.Val), color=group))+
    geom_point(size=2,alpha=0.5)+
    geom_hline(aes(yintercept=-log10(0.05)), color="red",alpha=0.5,linetype=2)+
    geom_vline(aes_(xintercept=-logFCcutoff), color="red", color="red",alpha=0.5,linetype=2)+
    geom_vline(aes_(xintercept=logFCcutoff), color="red", color="red",alpha=0.5,linetype=2)+
    theme_minimal()+
    theme(legend.position = "none")+
    labs(x=xlab, y=ylab, color="Change")+
    scale_color_manual(values=c("down"="#0073C2FF",
                                "none"="#868686FF",
                                "up"="#EFC000FF"))+
    geom_text_repel(aes(label=labels),size = 2.25,
                    segment.color = "black",
                    box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.1, "lines"))


  return(p)
}


#' plot GroupBar
#'
#'
#' @title plotGroupBar
#' @param pdata pdata
#' @param x a vector.
#' @param fill logFC cutoff
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggsci scale_fill_jco
#' @return a ggplot2 object
#' @export
#' @author Yuanlong Hu


plotGroupBar <- function(pdata, x, fill){
  pdata <- pdata[,c(x,fill)]
  names(pdata) <- c("x","fill")
  ggplot(pdata, aes(x=x, fill=fill))+
    geom_bar(stat="count",width=1,position='fill',color="white")+
    geom_text(stat='count',aes(label=scales::percent(..count../sum(..count..))),
              color="black", size=3.5,position=position_fill(0.5))+
    theme_minimal()+
    scale_fill_jco()
}

#' plot PCA and t-sne
#'
#'
#' @title plotExprDIM
#' @param expr expr
#' @param feature feature
#' @param pdata a vector.
#' @param methed "PCA" or "tSNE"
#' @param addEllipses TRUE or FALSE
#' @param ellipse_type "convex" or "confidence"
#' @param ellipse_level 0.95
#' @param tsne_perplexity perplexity
#' @param ... other
#' @importFrom FactoMineR PCA
#' @importFrom factoextra fviz_pca_ind
#' @importFrom Rtsne Rtsne
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 stat_ellipse
#' @importFrom ggplot2 labs
#' @importFrom ggsci scale_color_jco
#' @importFrom ggsci scale_fill_jco
#' @importFrom ggplot2 theme_minimal
#' @return a ggplot2 object
#' @export
#' @author Yuanlong Hu

plotExprDIM <- function(expr, feature, pdata=NULL,
                        methed=c("PCA","tSNE"),
                        addEllipses=TRUE,
                         ellipse_type="confidence",
                         ellipse_level=0.95,
                        tsne_perplexity=30,
                         ...){

  expr <- expr[feature,]
  expr <- na.omit(expr)
  message(paste("** A total of",nrow(expr), "features. **"))
  expr <- as.data.frame(t(expr))


  if (method[1]=="PCA") {


  res_pca <- PCA(expr, graph = FALSE)
  if(is.null(pdata)){
    p <- fviz_pca_ind(res_pca,
                                  col.ind = "blue",
                                  #palette = "jco",
                                  #addEllipses = addEllipses,
                                  label = "none",pointsize = 2.5,
                                  #col.var = "black",#alpha.ind = 0.5,
                                  alpha.var=0.7,
                                  #ellipse.level = ellipse.level,
                                  #gradient.cols = "RdYlBu",col.var = "black",
                                  repel = F,title = "",legend.title = "Group",...) +
      #theme(legend.position = "right")+
      theme_minimal()
  }else{



    if(ellipse_type=="confidence"){
      p <- fviz_pca_ind(res_pca,
                                    col.ind = factor(pdata),
                                    palette = "jco",
                                    addEllipses = addEllipses,
                                    label = "none",pointsize = 2.5,
                                    #col.var = "black",#alpha.ind = 0.5,
                                    alpha.var=0.7,ellipse.level = ellipse_level,
                                    #gradient.cols = "RdYlBu",col.var = "black",
                                    repel = F,title = "",legend.title = "Group",...) +
        #theme(legend.position = "right")+
        theme_minimal()
    }else{
      p <- fviz_pca_ind(res_pca,
                                    col.ind = factor(pdata),
                                    palette = "jco",
                                    addEllipses = addEllipses,
                                    label = "none",pointsize = 2.5,
                                    #col.var = "black",#alpha.ind = 0.5,
                                    alpha.var=0.7,ellipse.type=ellipse_type,
                                    ellipse.level=ellipse_level,
                                    #gradient.cols = "RdYlBu",col.var = "black",
                                    repel = F,title = "",legend.title = "Group",...) +
        #theme(legend.position = "right")+
        theme_minimal()
    }
  }
  }


  if(method[1]=="tSNE"){
    res_tsne <- Rtsne(expr,perplexity=tsne_perplexity)
    res_tsne <- data.frame(x = res_tsne$Y[,1], y = res_tsne$Y[,2], col = pdata)
    p <- ggplot(res_tsne,
                aes(x=x, y=y, color=col, fill=color)) +
      geom_point()+
      theme_minimal()+
      labs(x="Dim1",y="Dim2")+
      scale_fill_jco()+
      scale_color_jco()

    if(addEllipses){
      p <- p+stat_ellipse(level = ellipse_level,alpha=0.8)
    }
  }

  return(p)
}
