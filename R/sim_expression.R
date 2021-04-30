#' Calculate the similarity
#'
#'
#' @title sim_expr
#' @param data data
#' @param gene gene
#' @param method one of pearson and spearman.
#' @param MIC_pvalue MIC p-value
#' @param pcutoff <0.05
#' @param cut 2
#' @param n 10
#' @importFrom Hmisc rcorr
#' @importFrom pbapply pblapply
#' @importFrom minerva mine
#' @return a data.frame
#' @export
#' @author Yuanlong Hu



sim_expr <- function(data, gene, method=c("pearson","spearman","MIC"), MIC_pvalue=FALSE, pcutoff=0.05, cut=2, n=100){

  data <- as.data.frame(data)
  seed <- as.numeric(data[gene,])

  if (cut > (nrow(data)/2)-1) {
    stop("cut > (nrow(data)/2)-1)")
  }

  data$spilt <- cut(1:nrow(data), cut, labels=F)
  data <- split(data,data$spilt)

  if(method %in% c("pearson","spearman")){

    res <- pblapply(data, function(x){
      x <- x[,-ncol(x)]
      x <- as.data.frame(t(x))
      x$seed <- seed
      res <- rcorr(as.matrix(x), type=method[1])
      res <- convCorrMatrix(res$r,res$P)
      return(res)
  })

    res <- Reduce(rbind,res)
    res <- res[res$column=="seed",]
    res <- res[,-2]
    colnames(res) <- c("Object","Cor","Pvalue")
    res <- res[res$Pvalue < pcutoff,]
  }

  if(method =="MIC"){

    res <- pblapply(data, function(x){
      x <- x[,-ncol(x)]
      x <- as.data.frame(t(x))
      x$seed <- seed
      res <- mine(x,normalization = F, n.cores = 4)$MIC
      if(MIC_pvalue){
        p_value <- MIC_pvalue(x=x, res, n=n)
      }else{
        p_value <- res
        p_value <- as.matrix(p_value)
        p_value <- 0.00001
      }

      res <- convCorrMatrix(res,p_value)
      return(res)
    })

    res <- Reduce(rbind,res)
    res <- res[res$column=="seed",]
    res <- res[,-2]

    if(!MIC_pvalue){
      colnames(res) <- c("Object","MIC", "pvalue")
    }else{
      res <- res[,-3]
      colnames(res) <- c("Object","MIC")
    }

  }

  return(res)
}

#' Calculate the similarity
#'
#'
#' @title sim_expr
#' @param data data
#' @param gene gene
#' @param method one of pearson and spearman.
#' @param MIC_pvalue MIC p-value
#' @param pcutoff <0.05
#' @param cut 2
#' @param n 10
#' @importFrom Hmisc rcorr
#' @importFrom pbapply pblapply
#' @importFrom minerva mine
#' @return a data.frame
#' @export
#' @author Yuanlong Hu


sim_expr2 <- function(...,
                      expr=NULL, pdata=NULL,
                      method=c("pearson", "spearman")){
  if (is.null(expr)) {
    data <- list(...)
    data <- Reduce(rbind, data)
    data <- t(data)
  }else{
    data <- expr
  }


  if (!is.null(pdata)) data$pdata <- data$pdata

  res <- Hmisc::rcorr(data, type = method[1])
  return(res)
}



#' Conversion CorrMatrix
#'
#'
#' @title convCorrMatrix
#' @param mat CorrMatrix
#' @param pmat p-value matrix
#' @param pvalueCutoff p-value cutoff
#' @param corCutoff cor cutoff
#' @return a data.frame
#' @author Yuanlong Hu
#' @noRd

convCorrMatrix <- function(mat, pmat, pvalueCutoff=0.05, corCutoff=0.5) {
  ut <- upper.tri(mat)
  mat <- data.frame(
    row = rownames(mat)[row(mat)[ut]],
    column = rownames(mat)[base::col(mat)[ut]],
    cor  =(mat)[ut],
    p = pmat[ut]
  )
  mat <- mat[mat$p<pvalueCutoff & abs(mat$cor)>=corCutoff,]
  return(mat)
}


#' plot cor
#'
#'
#' @title plotCorHeatmap
#' @param res Hmis
#' @param x x
#' @param y y
#' @param pvalueCutoff p-value cutoff
#' @param mark_var mark "*"
#' @param style "A" or "B"
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 theme
#' @importFrom scales muted
#' @importFrom ggplot2 element_text
#' @return ggplot object
#' @author Yuanlong Hu
#' @export


plotCorHeatmap <- function(res, x, y,
                           pvalueCutoff=0.05){
  cormat <- res$r
  pmat <- res$P

  cormat <- as.data.frame(cormat[x,y])
  cormat$row <- rownames(cormat)
  pmat <- as.data.frame(pmat[x,y])
  pmat$row <- rownames(pmat)
  cormat <- reshape2::melt(cormat, id.vars="row")
  pmat <- reshape2::melt(pmat, id.vars="row")

  data <- cormat
  data$pvalue <- pmat$value

  data$marker <- ifelse(data$pvalue<pvalueCutoff,data$value ,"×")
  #data$marker <- ifelse(abs(data$value)>= mark_var,"*","" )

  data$row <- with(data, reorder(row, value, mean))
  data$variable <- with(data, reorder(variable, value, mean))

  ggplot(data, aes(x=row, y=variable))+
    geom_tile(aes(fill=value), color="gray")+
    geom_text(aes(label=marker),size=4)+
    labs(x="",y="")+
    scale_fill_gradient2(low = scales::muted("blue"),
                         high = scales::muted("red"))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle=45,
                                     hjust=1,
                                     vjust=1))
}



#' plot cor
#'
#'
#' @title plotExprCor
#' @param expr expr data
#' @param x x
#' @param y y
#' @param group group
#' @param point_group logical
#' @param geom_smooth logical
#' @param method c("lm","glm","gam","loess")
#' @param smooth_group logical
#' @param Cormethod "pearson" or "spearman"
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 aes
#' @importFrom ggpubr stat_cor
#' @return ggplot object
#' @author Yuanlong Hu
#' @export


plotExprCor <- function(expr, x, y,
                       group=NULL, point_group=TRUE,
                       geom_smooth=TRUE, method=c("lm","glm","gam","loess"),
                       smooth_group=TRUE, Cormethod="pearson"){


  df <- data.frame(x=as.numeric(expr[x,]),
                   y=as.numeric(expr[y,])
                   )
  df$group <- group
  p <- ggplot(df, aes(x=x, y=y))

  if (point_group) {
    p <- p+geom_point(aes(color=group))
  }else{
    p <- p+geom_point()
  }

  if(geom_smooth){

    if (smooth_group) {
      p <- p+geom_smooth(aes(color=group), method=method[1])+
        stat_cor(method=Cormethod)
    }else{
      p <- p+geom_smooth(method=method[1])
    }

  }

  p <- p+labs(x=x, y=y)

  # res_cor <- rcorr(as.matrix(df[,1:2]), type = "pearson")
  # message(paste0("Pearson: ",round(res_cor[[1]][1,2],3), "\n",
  #                "P-value: ",signif(res_cor[[3]][1,2],3), "\n",
  #                "N: ", res_cor[[2]][1,2]))
  return(p)
}

#' plot cor
#'
#'
#' @title MIC_pvalue
#' @param x data matrix
#' @param MIC_mat MIC matrix
#' @param n n
#' @importFrom minerva mine
#' @return p-value matrix
#' @author Yuanlong Hu
#' @export

MIC_pvalue <- function(x,MIC_mat, n){
  data <- x
  p_num <- MIC_mat
  p_num[abs(p_num)>0] <- 1

  set.seed(1234)
  for (i in 1:n) {
    random <- apply(data,2,sample)
    micN <- mine(random, normalization = F, n.cores = 4)$MIC
    micN[abs(micN)>=abs(MIC_mat)] <- 1
    micN[abs(micN)<abs(MIC_mat)] <- 0
    p_num <- p_num + micN
  }

  p <- p_num/(n+1)
}
