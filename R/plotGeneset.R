#' plot geneset heatmap
#'
#'
#' @title plotExprGenesetHeatmap
#' @param res result
#' @param select selected geneset name
#' @param genesetlist a geneset list
#' @param logFCCutoff abs logFC cutoff
#' @param pvalueCutoff p-value cutoff
#' @param adjpvalueCutoff adjusted p-value cutoff
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom scales muted
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 theme_minimal
#' @return a ggplot2 object
#' @export
#' @author Yuanlong Hu

plotExprGenesetHeatmap <- function(res, select, selectgene=NULL,genesetlist,
                                   logFCCutoff=0.2,
                                   pvalueCutoff=0.05,adjpvalueCutoff=0.05){
  res$Gene <- rownames(res)
  select <- as.list(select)
  genesetlist <- lapply(select, function(x){
    data.frame(Geneset=rep(x,length(genesetlist[[x]])),
               Gene=genesetlist[[x]])
  })

  genesetlist <- Reduce(rbind, genesetlist)
  genesetlist <- merge(genesetlist,res[,c("logFC","P.Value","adj.P.Val","Gene")],
                       by = "Gene")

  genesetlist <- genesetlist[abs(genesetlist$logFC)>= logFCCutoff,]
  genesetlist <- genesetlist[genesetlist$P.Value<pvalueCutoff,]
  genesetlist <- genesetlist[genesetlist$adj.P.Val<adjpvalueCutoff,]

  if(!is.null(selectgene)){
    genesetlist <- genesetlist[genesetlist$Gene %in% selectgene,]
  }

  genesetlist$Gene <- with(genesetlist, reorder(Gene, logFC,mean))
  genesetlist$Geneset <- with(genesetlist, reorder(Geneset, Gene,length))

  p <- ggplot(genesetlist,aes(x=Gene,y=Geneset))+
    geom_tile(aes(fill=logFC))+
    scale_fill_gradient2(low = muted("blue"), high = muted("red"))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle=45,hjust=1, vjust=1))

  return(p)
}

#' Heatmap displaying gene expression in GSEA
#'
#'
#' @title plotExprGSEAHeatmap
#' @param res_gsea GSEA result
#' @param res_deg DEGs result
#' @param select selected ID
#' @param sort sort method
#' @param logFCCutoff abs logFC cutoff
#' @param pvalueCutoff p-value cutoff
#' @param adjpvalueCutoff adjusted p-value cutoff
#' @param show_NES Whether to show GSEA NES
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom scales muted
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 theme_minimal
#' @return a ggplot2 object
#' @export
#' @author Yuanlong Hu

plotExprGSEAHeatmap <- function(res_gsea, res_deg, select, sort="logFC",
                                logFCCutoff=0.1, pvalueCutoff=0.05, adjpvalueCutoff=0.1,
                                show_NES=TRUE){

  res_gsea <- Reduce(rbind,lapply(res_gsea, function(x) as.data.frame(x)))

  res_gsea <- res_gsea[select,c("ID","Description","NES","core_enrichment")]
  res_deg$gene <- rownames(res_deg)
  res_gsea1 <- list()
  for (i in 1:nrow(res_gsea)) {
    core_enrichment <- unlist(stringr::str_split(res_gsea$core_enrichment[i],"/"))
    res_gsea0 <- res_deg[intersect(core_enrichment,rownames(res_deg)),c("gene","logFC","P.Value","adj.P.Val")]
    res_gsea0$pathway <- res_gsea$Description[i]
    res_gsea1 <- c(res_gsea1, list(res_gsea0))
  }
  res_gsea1 <- Reduce(rbind,res_gsea1)

  res_gsea1 <- res_gsea1[res_gsea1$P.Value < pvalueCutoff & res_gsea1$adj.P.Val <adjpvalueCutoff,]
  res_gsea1 <- res_gsea1[abs(res_gsea1$logFC) >= logFCCutoff,]

  if (sort=="logFC") {
    res_gsea1$gene <- with(res_gsea1, reorder(gene,logFC, mean))
  }

  p1 <- ggplot(res_gsea1, aes(x=gene, y=pathway))+
    geom_tile(aes(fill=logFC),color="gray")+
    scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle=45,hjust=1, vjust=1),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())+
    labs(x="",y="")

  if(show_NES){
    if (sort=="NES") {
      res_gsea$Description <- with(res_gsea,reorder(Description,NES, mean))
    }
    res_gsea$type <- "NES"
    p2 <- ggplot(res_gsea, aes(x=type, y=Description))+
      geom_tile(aes(fill=NES),color="gray")+
      theme_minimal()+
      scale_fill_gradient2(low = scales::muted("blue"),
                           high = scales::muted("red"))+
      labs(y="",x="")+
      theme(axis.text.y = element_blank(),
            axis.ticks=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())

    p3 <- aplot::insert_right(p1,p2, width=0.05)
    return(p3)
  }
}

#' Summary NES in GSEA results
#'
#'
#' @title plotSummaryGSEAHeatmap
#' @param reslist GSEA result list
#' @param selectID selected ID
#' @param show_top top
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom scales muted
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 theme_minimal
#' @return a ggplot2 object
#' @export
#' @author Yuanlong Hu


plotSummaryGSEAHeatmap <- function(reslist, selectID=NULL, show_top=20){


  reslist <- lapply(reslist, function(x){
    if (!is.data.frame(x)) {
      x <- lapply(x, as.data.frame)
      x <- Reduce(rbind, x)
    }
    return(x)
  })

  for(i in 1:length(reslist)){
    reslist[[i]]$Type <- names(reslist)[i]
  }
  reslist <- Reduce(rbind, reslist)
  reslist <- reslist[order(abs(reslist$NES), decreasing = T),]

  if(is.null(selectID)) selectID <- unique(reslist$ID)[1:show_top]

  reslist <- reslist[reslist$ID %in% selectID,]

  reslist$Description <- with(reslist, reorder(Description,NES,max))
  ggplot(reslist, aes(x=Type, y=Description))+
    geom_tile(aes(fill=NES),color="gray",alpha=0.8)+
    scale_fill_gradient2(low = scales::muted("blue"),
                         high = scales::muted("red"))+
    theme_minimal()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
}
