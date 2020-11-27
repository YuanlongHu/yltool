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


#' plot cor
#'
#'
#' @title plot_point
#' @param data expr data
#' @param x x
#' @param y y
#' @param group group
#' @param point_group logical
#' @param geom_smooth logical
#' @param method c("lm","glm","gam","loess")
#' @param smooth_group logical
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 aes
#' @return ggplot object
#' @author Yuanlong Hu
#' @export


plot_point <- function(data, x, y,
                       group=NULL, point_group=TRUE,
                       geom_smooth=TRUE, method=c("lm","glm","gam","loess"), smooth_group=TRUE){


  df <- data.frame(x=as.numeric(data[x,]),
                   y=as.numeric(data[y,])
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
      p <- p+geom_smooth(aes(color=group), method=method[1])
    }else{
      p <- p+geom_smooth(method=method[1])
    }

  }
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
