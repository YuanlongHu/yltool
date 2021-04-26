#' gsea
#'
#'
#' @title enrich_gsea
#' @param res res
#' @param pvalueCutoff p-value cutoff
#' @param kegg_internal_data Use kegg internal data
#' @param GMTset gmt file
#' @param useGMTset use geneset from gmt file
#' @param IDtoNAME ID to name
#' @param na.omit Remove NA
#' @importFrom clusterProfiler bitr
#' @importFrom clusterProfiler gseKEGG
#' @importFrom clusterProfiler gseGO
#' @importFrom clusterProfiler gseMKEGG
#' @importFrom ReactomePA gsePathway
#' @importFrom clusterProfiler setReadable
#' @importFrom clusterProfiler GSEA
#' @importFrom immcp read_gmt
#' @return a gseaResult object
#' @export
#' @author Yuanlong Hu

enrich_gsea <- function(res, pvalueCutoff=0.05,
                        kegg_internal_data=FALSE,
                        GMTset=NULL,
                        useGMTset=FALSE,
                        IDtoNAME=NA,
                        na.omit=TRUE){
  if(na.omit){
    res <- na.omit(res)
    res <- res[res$logFC != Inf,]
    res <- res[res$logFC != -Inf,]
    res <- res[res$logFC != " ",]
    res <- res[res$logFC != "",]
  }

  if(useGMTset){

    message("** Build Genelist **")
    genelist <- res$logFC
    names(genelist) <- rownames(res)
    genelist <- sort(genelist, decreasing = T)

    message("** Read GMT file **")
    if(is.data.frame(GMTset)) gmt <- GMTset[,1:2]
    if(is.character(GMTset)) gmt <- clusterProfiler::read.gmt(GMTset)
    if(is.list(GMTset)) gmt <- immcp::to_df(GMTset)

    message("** Run GSEA **")
    res <- clusterProfiler::GSEA(
      geneList=genelist,
      pvalueCutoff=0.05,
      minGSSize = 5,maxGSSize = 500,
      TERM2GENE=gmt,
      TERM2NAME = IDtoNAME)

  }else{

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

  message("** KEGG Module GSEA **")
  mkegg <- clusterProfiler::gseMKEGG(geneList = genelist,
                                     organism = 'hsa')


  message("** GO-BP GSEA **")
  ego_BP <- clusterProfiler::gseGO(geneList= genelist,
                                   OrgDb=org.Hs.eg.db,
                                   ont= "BP",nPerm=10000,
                                   minGSSize=5,maxGSSize = 500,
                                   pvalueCutoff = pvalueCutoff,
                                   verbose= FALSE)
  message("** GO-CC GSEA **")
  ego_CC <- clusterProfiler::gseGO(geneList= genelist,
                                   OrgDb=org.Hs.eg.db,
                                   ont= "CC",nPerm=10000,
                                   minGSSize=5,maxGSSize = 500,
                                   pvalueCutoff = pvalueCutoff,
                                   verbose= FALSE)
  message("** GO-MF GSEA **")
  ego_MF <- clusterProfiler::gseGO(geneList= genelist,
                                   OrgDb=org.Hs.eg.db,
                                   ont= "MF",nPerm=10000,
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

  message("** Wikipathways GSEA **")
  genesetlist <- prepareGeneset("wikipathways")
  Wikipathways <- GSEA(geneList = genelist,
                       TERM2GENE = genesetlist$TERM2GENE,
                       TERM2NAME = genesetlist$TERM2NAME,
                       verbose = FALSE)

  message("** Summary Result **")
  res <- list(kegg=kegg,
              mkegg=mkegg,
              go_CC=ego_CC,
              go_MF=ego_MF,
              go_BP=ego_BP,
              Reactome=Reactome,
              Wikipathways=Wikipathways)
  }

  res <- lapply(res, function(x){
    x <- clusterProfiler::setReadable(x,
                                      OrgDb = org.Hs.eg.db,
                                      keyType="ENTREZID")
  })
  return(res)
}

#' prepare Geneset
#'
#'
#' @title prepareGeneset
#' @param geneset Name
#' @importFrom clusterProfiler read.gmt
#' @importFrom tidyr separate
#' @importFrom dplyr select
#' @importFrom vroom vroom
#' @importFrom tidyr unite
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @return a list
#' @export
#' @author Yuanlong Hu

prepareGeneset <- function(geneset){

  if(geneset == "wikipathways"){
    wp2gene <- clusterProfiler::read.gmt("https://wikipathways-data.wmcloud.org/current/gmt/wikipathways-20210410-gmt-Homo_sapiens.gmt")
    wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
    wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
    wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

  genesetlist <- list(TERM2GENE=wpid2gene,
                   TERM2NAME=wpid2name)
  }

  if(geneset=="CellMarker"){
    cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
      tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>%
      dplyr::select(cellMarker, geneID) %>%
      dplyr::mutate(geneID = strsplit(geneID, ', '))

    genesetlist=list(TERM2GENE=cell_markers)

  }
  return(genesetlist)

}
