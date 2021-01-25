#' gsea
#'
#'
#' @title enrich_gsea
#' @param res res
#' @param pvalueCutoff p-value cutoff
#' @param kegg_internal_data Use kegg internal data
#' @importFrom clusterProfiler bitr
#' @importFrom clusterProfiler gseKEGG
#' @importFrom clusterProfiler gseGO
#' @importFrom ReactomePA gsePathway
#' @importFrom clusterProfiler setReadable
#' @return a ggplot2 object
#' @export
#' @author Yuanlong Hu

enrich_gsea <- function(res, pvalueCutoff=0.05, kegg_internal_data=FALSE){


  message("** Biological Id translation **")
  res$SYMBOL <- rownames(res)
  deg_SYMBOL_ENTREZID <- clusterProfiler::bitr(res$SYMBOL,fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  res <- merge(res, deg_SYMBOL_ENTREZID, by="SYMBOL")

  genelist <- res$logFC
  names(genelist) <- res$ENTREZID
  genelist <- sort(genelist, decreasing = T)

  message("** KEGG GSEA **")
  kegg <- clusterProfiler::gseKEGG(geneList=genelist,organism="hsa",
                                   minGSSize=5, maxGSSize=500,
                                   nPerm = 10000, pAdjustMethod = "fdr",
                                   pvalueCutoff = pvalueCutoff,
                                   verbose = TRUE,
                                   use_internal_data = kegg_internal_data)
  message("** GO-BP GSEA **")
  ego_BP <- clusterProfiler::gseGO(geneList= genelist,
                                   OrgDb=org.Hs.eg.db,
                                   ont= "BP",nPerm=10000,
                                   minGSSize=5,maxGSSize = 500,
                                   pvalueCutoff = pvalueCutoff,
                                   verbose= FALSE)
  message("** Reactome GSEA **")
  Reactome <- ReactomePA::gsePathway(geneList=genelist, organism="human",
                                     exponent= 1, nPerm=10000,
                                     minGSSize= 5, maxGSSize= 500,
                                     pvalueCutoff=pvalueCutoff,
                                     pAdjustMethod= "fdr",
                                     verbose=TRUE, seed = FALSE,
                                     by = "fgsea")

  message("** Biological Id translation2 **")
  kegg <- clusterProfiler::setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  ego_BP <- clusterProfiler::setReadable(ego_BP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  Reactome <- clusterProfiler::setReadable(Reactome, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

  message("** Summary Result **")
  res <- list(kegg=kegg,
              ego_BP=ego_BP,
              Reactome=Reactome)
  return(res)
}
