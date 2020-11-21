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

#' Convert probe ID to gene symbol using GPL
#'
#'
#' @title ann_expr2
#' @param data data
#' @param GPL GPL file
#' @param probe_symbol probe id and gene symbol
#' @param Stat "max","mean","median","min","IQR"
#' @importFrom Biobase exprs
#' @return a data.frame
#' @export
#' @author Yuanlong Hu

ann_expr2 <- function(data, GPL, probe_symbol=c("NAME","GENE_SYMBOL"),
                      Stat=c("max","mean","median","min","IQR")){

  # read gpl file
  gpl <- read.table(GPL,
                 header = TRUE, fill = T,sep = "\t",
                 comment.char = "#",
                 stringsAsFactors = FALSE,
                 quote = "")
  gpl <- gpl[,probe_symbol]
  colnames(gpl) <- c('probe_id','symbol')



  datExpr <- exprs(data)
  datExpr <- as.data.frame(datExpr)
  datExpr$probe_id <- rownames(datExpr)

  datExpr <- merge(gpl, datExpr, by = "probe_id")

  datExpr <- datExpr[datExpr$symbol != "",]

  if(Stat[1]=="max") datExpr <- aggregate(datExpr, by=list(datExpr$symbol), max)
  if(Stat[1]=="mean") datExpr <- aggregate(datExpr, by=list(datExpr$symbol), mean)
  if(Stat[1]=="median") datExpr <- aggregate(datExpr, by=list(datExpr$symbol), median)
  if(Stat[1]=="min") datExpr <- aggregate(datExpr, by=list(datExpr$symbol), min)
  if(Stat[1]=="IQR") datExpr <- aggregate(datExpr, by=list(datExpr$symbol), IQR)


  #colnames(data_GSE47460_GPL14550)

  rownames(datExpr) <- datExpr$symbol
  datExpr <- datExpr[,-c(1:3)]
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
