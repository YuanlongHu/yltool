#' Scoring ZhengSu
#'
#'
#' @title ZSscore
#' @param data data
#' @param power power
#' @return a data.frame
#' @export
#' @author Yuanlong Hu

ZSscore <- function(data,power){
  ZS2 <- data

  message("---- Scoring ----")
  dat0 <- list()
  for (i in 1:ncol(ZS2)) {
    dat <- data.frame(zz=rownames(ZS2),value= ZS2[,i])
    dat <- merge(dat,power,by="zz")
    dat$score <- dat$value*dat$power2
    dat <- dat[,-c(1:2,4:5)]
    dat <- aggregate(dat[,2], by=list(dat$zs), sum)
    names(dat) <- c("zs",colnames(ZS2)[i])
    dat0 <- c(dat0, list(dat))
  }
  dat1 <- Reduce(merge,dat0, list(by="zs"))
  rownames(dat1) <- dat1$zs
  dat1 <- dat1[,-c(1:2)]
  dat1 <- as.data.frame(t(dat1))

  message("---- Conversion ----")
  dat2 <- apply(dat1, 2, function(x){
    x <- ifelse(x<70,0,
                ifelse(x<100,1,
                       ifelse(x<150,2,3))
    )
    return(x)
  })

  return(dat2)
}
