#' Convert probe ID to gene symbol
#'
#'
#' @title ann_expr
#' @param data data
#' @param ann_db The name of the annotation package
#' @param ann_df gene symbol data frame
#' @param Stat "max","mean","median","min","IQR"
#' @importFrom Biobase featureNames
#' @importFrom genefilter findLargest
#' @importFrom Biobase exprs
#' @return a data.frame
#' @export
#' @author Yuanlong Hu


ann_expr <- function(data, ann_db, ann_df,
                     Stat=c("max","mean","median","min","IQR")){

  if(Stat[1]=="max") testStat <- apply(exprs(data), 1, max)
  if(Stat[1]=="mean") testStat <- apply(exprs(data), 1, mean)
  if(Stat[1]=="median") testStat <- apply(exprs(data), 1, median)
  if(Stat[1]=="min") testStat <- apply(exprs(data), 1, min)
  if(Stat[1]=="IQR") testStat <- apply(exprs(data), 1, IQR)

  prbs <- findLargest(featureNames(data),
                      testStat = testStat,
                      data = ann_db)

  data <- data[prbs,]
  datExpr <- exprs(data)
  datExpr <- as.data.frame(datExpr)
  datExpr$probe_id <- rownames(datExpr)

  datExpr <- merge(toTable(ann_df), datExpr, by = "probe_id")
  rownames(datExpr) <- datExpr$symbol
  datExpr <- datExpr[,-c(1:2)]
  return(datExpr)

}


#' Add a row
#'
#'
#' @title add_row
#' @param data data
#' @param add_row a new row
#' @param rowname rowname
#' @return a data.frame
#' @export
#' @author Yuanlong Hu

add_row <- function(data, add_row, rowname){
  data <- rbind(add_row, data)
  rownames(data)[1] <- rowname
}
