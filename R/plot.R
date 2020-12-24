#' plot point
#'
#'
#' @title plot_cor
#' @param expr a matrix of expression values where rows correspond to genes and columns correspond to samples.
#' @param gene.x x
#' @param gene.y y
#' @param geom_smooth rna-seq data or microarray data
#' @param method one of "lm", "loess", "glm"
#' @param group group
#' @param group_point group_point
#' @param group_smooth group_smooth
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 labs
#' @importFrom Hmisc rcorr
#' @return a ggplot
#' @export
#' @author Yuanlong Hu


plot_cor <- function(expr,gene.x, gene.y,
                     geom_smooth=TRUE, method="lm",
                     group=NULL,
                     group_point=TRUE, group_smooth=TRUE
                     ){


  data <- data.frame(x=as.numeric(expr[gene.x,]),
                     y=as.numeric(expr[gene.y,]))

  if (!is.null(group)) {
    data$group <- group
  }
  p <- ggplot(data, aes(x=x, y=y))

  if(group_point){
    p <- p + geom_point(aes(color=group))
  }else{
    p <- p + geom_point()
  }

  p <- p + labs(x=gene.x, y=gene.y)+
       theme_minimal()

  if (group_smooth) {
    p <- p + geom_smooth(method = method)
  }else{
    p <- p + geom_smooth(aes(color=group),method = method)
  }

  res_cor <- rcorr(data[,1:2], type = "pearson")
  message(paste0(res_cor[[1]][1,2], "\n",res_cor[[2]][1,2]))
  return(p)
}



