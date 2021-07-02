#' graph node analysis
#'
#'
#' @title graph_node
#' @param data data
#' @param input "df","adjmat"
#' @param weighted weighted
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph E
#' @importFrom igraph V
#' @importFrom igraph degree
#' @importFrom igraph closeness
#' @importFrom igraph betweenness
#' @importFrom igraph evcent
#' @importFrom igraph strength
#' @return a data frame
#' @export
#' @author Yuanlong Hu


graph_node <- function(data, input=c("df","adjmat"), weighted=FALSE){

  if(input[1]=="adjmat"){
    igraph <- graph_from_adjacency_matrix(as.matrix(data),
                                          mode = 'undirected',
                                          weighted = ifelse(weighted,TRUE, NULL),
                                          diag = FALSE)
    if(weighted){
      E(igraph)$corr <- E(igraph)$weight
      E(igraph)$weight <- abs(E(igraph)$weight)
    }

  }else{
    if(input[1]=="df" & weighted) stop("Must be adjacency matrix")

    igraph <- graph_from_data_frame(d=as.data.frame(data),
                                    directed = FALSE,
                                    vertices = NULL)

  }

  node_list <- data.frame(
    node_id = V(igraph)$name,
    degree = degree(igraph),
    #weight_degree = V(igraph)$weight_degree,
    closeness_centrality = closeness(igraph),
    betweenness_centrality = betweenness(igraph),
    eigenvector_centrality = evcent(igraph)$vector
  )

  if(weighted) node_list$weight_degree <- strength(igraph)

  return(node_list)
}
