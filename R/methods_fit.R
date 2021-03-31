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


#' Fit single logistic methods
#'
#'
#' @title select_logistic_single
#' @param expr a matrix of expression values where rows correspond to genes and columns correspond to samples.
#' @param feature feature
#' @param pdata a vector.
#' @importFrom pbapply pblapply
#' @importFrom ROCR prediction
#' @importFrom ROCR performance
#' @return a data.frame
#' @export
#' @author Yuanlong Hu

select_logistic_single <- function(expr, feature=NULL, pdata){
  if(is.null(feature)){
    feature <- rownames(expr)
  }
  expr <- as.data.frame(t(expr[feature,]))
  expr$pdata <- as.factor(pdata)
  formula <- paste0("pdata~",feature)
  formula <- as.list(formula)
  res <- pblapply(formula, function(x){
    fit <- suppressMessages(glm(as.formula(x), data = expr, family = binomial()))
    auc <- predict(fit, type = "response")
    auc <- ROCR::prediction(auc, expr$pdata)
    auc <- ROCR::performance(auc, "auc")@y.values[[1]]
    CI <- confint(fit)
    fit <- as.data.frame(summary(fit)$coefficients)
    res <- c(exp(fit[2,1]),fit[2,4],exp(CI[2,1]), exp(CI[2,2]),auc)
    names(res) <- c("OR","pvalue","CI_lower","CI_upper","AUC")
    res <- signif(res,2)

    return(res)
  })
  names(res) <- feature

  res <- as.data.frame(t(as.data.frame(res)))
  return(res)

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
#' @importFrom survival coxph
#' @importFrom survival Surv
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
#' @param pValue confidence level
#' @param maxRuns maximal number of importance source runs.
#' @importFrom Boruta Boruta
#' @return a Boruta object
#' @export
#' @author Yuanlong Hu

select_Boruta <- function(expr, pdata, pValue=0.01, maxRuns=100){
  expr <- as.data.frame(t(expr))
  if (is.character(pdata)) {
    pdata <- as.factor(pdata)
  }
  expr$pdata <- pdata
  res <- Boruta::Boruta(pdata~., expr, doTrace=1,
                        pValue = pValue, maxRuns = maxRuns)
  return(res)
}


#' plot Boruta ImpHistory
#'
#'
#' @title plotBorutaImpHistory
#' @param res the result of Boruta
#' @param select a vector.
#' @importFrom Boruta Boruta
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_fill_manual
#' @return a Boruta object
#' @export
#' @author Yuanlong Hu

plotBorutaImpHistory <- function(res,
                                 select = NULL,
                                 select2 = c("Confirmed", "Rejected", "Tentative")){

  ImpHistory <- as.data.frame(t(res$ImpHistory))
  ImpHistory$Group <- res$finalDecision[rownames(ImpHistory)]
  ImpHistory$gene <- gsub("`","",rownames(ImpHistory))

  ImpHistory <- reshape2::melt(ImpHistory, id.vars=c("gene","Group"))

  ImpHistory$gene <- gsub("`","",ImpHistory$gene)
  if (is.null(select)) {
    ImpHistory <- ImpHistory[ImpHistory$Group %in% select2,]
  }else{
    ImpHistory <- ImpHistory[ImpHistory$gene %in% select & ImpHistory$Group %in% select2,]
  }

  #ImpHistory$gene <- gsub("."," ",ImpHistory$gene)
  ImpHistory$gene <- with(ImpHistory, reorder(gene, value, median))

  p <- ggplot(ImpHistory, aes(x=value, y=gene, fill=Group))+
    geom_boxplot()+
    theme_minimal()+
    labs(x="Importance",y="Attributes",fill="")+
    theme(legend.justification=c(1,0), legend.position=c(1,0))+
    scale_fill_manual(values = c("Confirmed"=="#EFC000FF",
                                 "Tentative"="#0073C2FF",
                                 "Rejected"="#868686FF"))
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

plotExprBox <- function(expr, select, pdata,
                        comparisons=list(c("C1","C2")),
                        label=c("p.signif","p.format")){

  expr <- data.frame(object=as.numeric(expr[select,]),
                     group=pdata)
  p <- ggplot(expr, aes(x=group, y=object, fill=group))+
    geom_violin()+
    geom_boxplot(width=0.2, fill="white")+
    scale_fill_jco()+
    stat_compare_means(label=label[1], comparisons = comparisons)+
    theme_minimal()+
    labs(y=select, x="")
  return(p)

}

#' plot ExprBox2
#'
#'
#' @title plotExprBox2
#' @param expr the expr data
#' @param select a vector.
#' @param pdata pdata
#' @param comparisons a list
#' @param label "p.signif" or "p.format"
#' @param style "A" or "B"
#' @param ncol the number of col
#' @param nrow the number of now
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_violin
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 labs
#' @importFrom ggsci scale_fill_jco
#' @importFrom ggsci scale_color_jco
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggpubr stat_compare_means
#' @return a ggplot2 object
#' @export
#' @author Yuanlong Hu
plotExprBox2 <- function(expr, select, pdata,
                         comparisons = list(c("S1", "S2")),
                         label=c("p.signif","p.format"),
                         style="A",
                         ncol=2, nrow = 2
){
  select <- intersect(rownames(expr),select)
  expr <- expr[select,]
  expr <- as.data.frame(t(expr))
  expr$group <- pdata

  expr <- reshape2::melt(expr, id.vars="group")

  if(style=="A"){
    p <- ggplot(expr, aes(x=group, y=value, fill=group))+
    geom_violin()+
    geom_boxplot(width=0.2, fill="white")+
    scale_fill_jco()+
    stat_compare_means(label=label[1], comparisons = comparisons)+
    theme_minimal()+
    labs(x="",y="")+
    facet_wrap(variable~.,scales = "free", ncol=ncol, nrow=nrow)
  }

  if(style=="B"){
  expr$variable <- with(expr, reorder(variable, value, mean))
  p <- ggplot(expr, aes(x=variable, y=value,
                        fill=group))+
    geom_boxplot(width=0.8, alpha=0.6)+
    stat_compare_means(aes(group=group), label = label[1])+
    scale_fill_jco()+
    #scale_color_jco()+
    theme_minimal()+
    labs(x="",y="")+
    coord_flip()
  }
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
#' @importFrom ggplot2 aes_
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 unit
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
                    color="black",
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
#' @importFrom ggplot2 position_fill
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
#' @param method "PCA" or "tSNE"
#' @param addEllipses TRUE or FALSE
#' @param ellipse_type "convex" or "confidence"
#' @param ellipse_level 0.95
#' @param tsne_perplexity perplexity
#' @param alpha 0.7
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
                        method,
                        addEllipses=TRUE,
                        ellipse_type="confidence",
                        ellipse_level=0.95,
                        tsne_perplexity=30,
                        alpha=0.7,
                        ...){

  expr <- expr[feature,]
  expr <- na.omit(expr)
  message(paste("** A total of",nrow(expr), "features. **"))
  expr <- as.data.frame(t(expr))


  if (method[1]=="PCA") {
    message("** Run PCA **")

  res_pca <- PCA(expr, graph = FALSE)
  if(is.null(pdata)){
    p <- fviz_pca_ind(res_pca,
                      col.ind = "blue",
                      #palette = "jco",
                      #addEllipses = addEllipses,
                      label = "none",pointsize = 2.5,
                       #col.var = "black",#alpha.ind = 0.5,
                      alpha.var=alpha,
                                  #ellipse.level = ellipse.level,
                                  #gradient.cols = "RdYlBu",col.var = "black",
                       repel = F, title = "", legend.title = "Group",...) +
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
                                    alpha.var=alpha,ellipse.level = ellipse_level,
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
                                    alpha.var=alpha,ellipse.type=ellipse_type,
                                    ellipse.level=ellipse_level,
                                    #gradient.cols = "RdYlBu",col.var = "black",
                                    repel = F,title = "",legend.title = "Group",...) +
        #theme(legend.position = "right")+
        theme_minimal()
    }
  }
  }


  if(method=="tSNE"){
    message("** Run t-SNE **")
    res_tsne <- Rtsne(expr,perplexity=tsne_perplexity)
    res_tsne <- data.frame(x = res_tsne$Y[,1],
                           y = res_tsne$Y[,2],
                           col = pdata)
    p <- ggplot(res_tsne,
                aes(x=x, y=y, color=col)) +
      geom_point()+
      theme_minimal()+
      labs(x="Dim1",y="Dim2")


    if(addEllipses){
      p <- p+
        stat_ellipse(aes(fill = col), geom = "polygon",
                     alpha = alpha, level = ellipse_level)+
        scale_fill_jco()+
        scale_color_jco()
    }
  }

  return(p)
}


