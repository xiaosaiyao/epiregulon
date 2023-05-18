#' @importFrom igraph E V
#' @export
findPartners <- function(graph, focal_tf){
  checkmate::assert_character(focal_tf, len = 1)
  checkmate::assert_class(graph, "igraph")
  if(length(unique(get.vertex.attribute(graph, name = "type"))) != 2)
    stop("Graph should contain vertices of two types. Input graph should not be created with the use of tripartite mode")
  target_type <- setdiff(unique(get.vertex.attribute(graph, name = "type")), "transcription factor")
  focal_vertex <- V(graph)[V(graph)$name == focal_tf & V(graph)$type == "transcription factor"]
  if(length(focal_vertex) == 0) stop(sprintf("Focal transcription factor (%s) is not present in the input graph"), focal_tf)
  focal_tf_targets <- igraph::ego(graph, node = focal_vertex, mode = "out")[[1]]
  focal_tf_targets <- focal_tf_targets[focal_tf_targets$type == target_type]
  redundant_vertices <- setdiff(as.numeric(V(graph)), as.numeric(V(graph)[V(graph)$type=="transcription factor"]))
  redundant_vertices <- setdiff(redundant_vertices, as.numeric(focal_tf_targets))
  graph <- igraph::delete.vertices(graph, redundant_vertices)
  if(any(duplicated(igraph::ends(graph, es=E(graph))))) stop("Duplicated edges")
  # find focal tf in the new graph
  focal_vertex <- V(graph)[V(graph)$name == focal_tf & V(graph)$type == "transcription factor"]
  focal_tf_targets <- igraph::ego(graph, node = focal_vertex, mode = "out")[[1]]
  focal_tf_targets <- focal_tf_targets[focal_tf_targets$type == target_type]
  other_tfs <- V(graph)[V(graph)$type == "transcription factor" & V(graph)$name != focal_tf]
  all_tfs <- c(focal_vertex, other_tfs)
  tf_targets <- ego(graph, nodes = all_tfs, mindist = 1, mode = "out")
  tf_targets_subgraphs <- mapply(function(x,y) subgraph(graph, vids = c(x, y)), all_tfs, tf_targets)
  res_list <- list()
  focal_weights <- get.edge.attribute(tf_targets_subgraphs[[1]], index =  E(tf_targets_subgraphs[[1]]), name = "weight")
  for(i in seq_along(other_tfs)){
    # skip focal tf
    if(i == 1) next
    tf <- other_tfs[i]
    common_targets <- ends(tf_targets_subgraphs[[i]], E(tf_targets_subgraphs[[i]]))[,2]
    # find indices of common targets in all targets of focal tf
    common_targets_ind <- match(common_targets, ends(tf_targets_subgraphs[[1]], E(tf_targets_subgraphs[[1]]))[,2])
    res_list[[tf$name]] <- data.frame(target = common_targets,
                                      focal_weigth = focal_weights[common_targets_ind],
                                      other_tf_weight = get.edge.attribute(tf_targets_subgraphs[[i]], index = E(tf_targets_subgraphs[[i]]), name = "weight"))
    res_list[[tf$name]]$weight_product <- res_list[[tf$name]]$focal_weigth * res_list[[tf$name]]$other_tf_weight

  }
  res_list
}


