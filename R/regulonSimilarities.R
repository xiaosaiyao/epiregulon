#' Find interaction partners of a transcription factor of interest
#' @param graph a igraph object from `buildGraph` or `buildDiffGraph`
#' @param focal_tf character string indicating the name of the transcription factors to find
#' interaction partners of
#' @export
#' @return A list with elements corresponding to each transcription factor apart from
#' the focal one. Each list element is represented as a data frame with columns containing
#' names of all target genes shared with focal transcription factor, weights of edges
#' connecting transcription factor with target genes, equivalent weights for focal transcription
#' factor and the element wise product of both weight columns.
findPartners <- function(graph, focal_tf){
  checkmate::assert_character(focal_tf, len = 1)
  checkmate::assert_class(graph, "igraph")
  if(length(unique(get.vertex.attribute(graph, name = "type"))) != 2)
    stop("Graph should contain vertices of two types. Input graph should not be created with the use of tripartite mode")
  focal_vertex <- V(graph)[V(graph)$name == focal_tf & V(graph)$type == "transcription factor"]
  if(length(focal_vertex) == 0) stop(sprintf("Focal transcription factor (%s) is not present in the input graph"), focal_tf)
  focal_tf_targets <- igraph::ego(graph, node = focal_vertex, mode = "out", mindist = 1)[[1]]
  # keep other transcription factors
  redundant_vertices <- setdiff(as.numeric(V(graph)), as.numeric(V(graph)[V(graph)$type=="transcription factor"]))
  # keep focal tf's target genes
  redundant_vertices <- setdiff(redundant_vertices, as.numeric(focal_tf_targets))
  # prune network removing target genes which do not belong to the focal tf's regulon
  graph <- igraph::delete.vertices(graph, redundant_vertices)
  if(any(duplicated(igraph::ends(graph, es=E(graph))))) stop("Duplicated edges")
  # find focal tf in the new graph
  focal_vertex <- V(graph)[V(graph)$name == focal_tf & V(graph)$type == "transcription factor"]
  other_tfs <- V(graph)[V(graph)$type == "transcription factor" & V(graph)$name != focal_tf]
  all_tfs <- c(focal_vertex, other_tfs)   # focal tf first on the list
  focal_tf_targets <- V(graph)[V(graph)$type == "target gene"]
  # create regulons using focal tf's target genes
  tf_targets_subgraphs <- lapply(all_tfs, function(tf, targets) igraph::subgraph(graph, vids = c(tf, intersection(targets, subcomponent(graph, tf, mode ="out")))), focal_tf_targets)
  res_list <- list()
  focal_weights <- igraph::get.edge.attribute(tf_targets_subgraphs[[1]], index =  E(tf_targets_subgraphs[[1]]), name = "weight")
  for(i in seq_along(all_tfs)){
    # skip focal tf
    if(i == 1) next
    tf <- V(tf_targets_subgraphs[[i]])[V(tf_targets_subgraphs[[i]])$type=="transcription factor" & V(tf_targets_subgraphs[[i]])$name==all_tfs[i]$name]
    common_targets <- V(tf_targets_subgraphs[[i]])[V(tf_targets_subgraphs[[i]])$type=="target gene"]
    edges_to_tf <- shortest_paths(tf_targets_subgraphs[[i]], from=tf, to=common_targets, output="epath")
    edges_to_tf <- do.call(c, edges_to_tf)
    # find indices of common targets in all targets of focal tf
    common_targets_ind <- match(common_targets$name, focal_tf_targets$name)
    res_list[[tf$name]] <- data.frame(target = common_targets$name,
                                      focal_weight = focal_weights[common_targets_ind],
                                      other_tf_weight = igraph::get.edge.attribute(tf_targets_subgraphs[[i]],
                                                                                   index = edges_to_tf,
                                                                                   name = "weight"))
    res_list[[tf$name]]$weight_product <- res_list[[tf$name]]$focal_weight * res_list[[tf$name]]$other_tf_weight
  }
  res_list
}

#' Calculate Jaccard Similarity between regulons of all transcription factors
#' @param graph a igraph object from `buildGraph` or `buildDiffGraph`
#' @export
#' @return A matrix with Jaccard similarity between all pairs of transcription factors.
calculateJaccardSimilarity <- function(graph){
  all_tfs <- V(graph)[V(graph)$type == "transcription factor"]
  res <- similarity(graph, vids = all_tfs, method = "jaccard", mode = "out")
  rownames(res) <- colnames(res) <- V(graph)[all_tfs]$name
  res
}


#' Calculate similarity score from permuted graphs to estimate background similarity
#' @param graph an igraph object from `buildGraph` or `buildDiffGraph`
#' @param focal_tf character string indicating the name of the transcription factors to
#' calculate similarity score
#' @param n an integer indicating the number of permutations
#' @param p a scalar indicating the probability of rewiring the graphs
#' @return A matrix with Jaccard similarity between the focal transcription factor and all pairs of transcription factors
#' for n permuted graphs
#' @export
permuteGraph <- function(graph, focal_tf, n=100, p=1){
  if(!focal_tf %in% names(V(graph)[V(graph)$type == "transcription factor"])) stop(focal_tf, " vertex shoud be present in the graph")
  all_tfs <- names(V(graph)[V(graph)$type == "transcription factor"])
  permute_matrix <- matrix(rep(NA, length(all_tfs)*n), nrow=length(all_tfs))
  rownames(permute_matrix) <-  all_tfs

  for ( iteration in seq_len(n)) {
    diff_graph_permute <- rewire(graph, each_edge(prob = p))
    similarity_score <- calculateJaccardSimilarity(diff_graph_permute)
    permute_matrix[, iteration ] <- similarity_score[focal_tf, all_tfs]
  }
  permute_matrix

}
