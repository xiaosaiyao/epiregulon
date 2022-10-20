#' Creating graphs and related operations
#'
#' @description
#' \code{build_graph} function creates a directed graph based on the output of
#' the \code{getRegulon} function.
#'
#' \code{build_difference_graph} a graph difference by subtracting the edges of \code{graph_obj_2}
#' from those of the \code{graph_obj_1}. If \code{weighted} is set to \code{TRUE} then for each
#' ordered pair of vertices (nodes) the difference in number of edges between \code{graph_obj_1}
#' and \code{graph_obj_1} is calculated. The result is used to set the number of
#' corresponding edges in output graph. Note that unless \code{abs_diff} is set to
#' \code{TRUE} any non-positive difference will translate into lack of the edges
#' for a corresponding ordered pair of vertices in the output graph (equivalent
#' to 0 value in the respective position in adjacency matrix). In case of
#' weighted graphs, the weight of the output graph is calculated as a difference
#' of the corresponding weights between input graphs.
#'
#' \code{add_centrality_degree} calculates degree centrality for each vertex using
#' \code{igraph::strength}.
#'
#' With \code{normalize_centrality} function the normalized values of centrality
#' are calculated from the original ones divided by
#' \code{FUN}(total number of non-zero edges associated with each node).
#'
#' \code{rank_tfs} assign ranks to transcription factors according to degree
#' centrality of their vertices
#'
#' @param regulon an object returned by the getRegulon or addWeights function
#' @param mode a character specifying whch type of graph will be built. In \code{'tg'} mode
#' a bipartite graph is built connecting transcription factors directly to the target genes
#' and ignoring information about mediating regulatory elements; in \code{'pairs'} mode
#' transcription factors are connected to unique target gene-regulatory element pairs;
#' in \code{'tripartite'} mode the network is build of three type of vertices (nodes):
#' trascription factors, regulatory elements and target genes; here the path from
#' target gene to regulatory element always contains a regulatory element; in
#' \code{'re'} mode data in the target genes is dropped and only connections are
#' between transcription factors and regulatory elements.
#' @param graph,graph_obj_1,graph_obj_2  an igraph object
#' @param weights a character specifying which variable should be used to assign
#' weights to edges. If set to 'NA' then unweighted graph is built.
#' @param aggregation_function a function used to collapse duplicated edges
#' @param weighted a logical indicating whether weighted graphs are used
#' @param abs_diff a logical indicating whether absulute difference in the number
#' edges or their weights will be calculated
#' @param FUN a function used for normalization. The input to this
#' function is be the number of edges connected with each node (incident edges).
#' @param type_attr a character corresponding to the name of the vertex attribute
#' which indicate the type of vertex
#' @return an igraph object
#' \code{rank_tfs} returns a data.frame with transcription factors sorted according
#' to the value of the \code{centrality} attribute
#' @examples
#' # create an artificial getRegulon output
#' set.seed(1234)
#' tf_set <- apply(expand.grid(LETTERS[1:10], LETTERS[1:10]),1,  paste, collapse = "")
#' regulon <- data.frame(tf = sample(tf_set, 5e3, replace = TRUE))
#' gene_set <- expand.grid(LETTERS[1:10], LETTERS[1:10], LETTERS[1:10])
#' gene_set <- apply(gene_set,1,function(x) paste0(x,collapse=""))
#' regulon$target <- sample(gene_set, 5e3, replace = TRUE)
#' regulon$idxATAC <- 1:5e3
#' regulon$corr <- runif(5e3)*0.5+0.5
#' graph_tripartite <- build_graph(regulon, mode = "tripartite")
#' # build bipartite graph using regulatory element-target gene pairs
#' graph_pairs_1 <- build_graph(regulon, mode = "pairs")
#' regulon$corr <- runif(5e3)*0.5+0.5
#' graph_pairs_2 <- build_graph(regulon, mode = "pairs")
#' graph_diff <- build_difference_graph(graph_pairs_1, graph_pairs_2)
#' graph_diff <- add_centrality_degree(graph_diff)
#' graph_diff <- normalize_centrality(graph_diff)
#' tf_ranking <- rank_tfs(graph_diff)
#' @importFrom igraph graph_from_data_frame V V<- E E<- vcount strength incident_edges
#' list.edge.attributes graph_from_adjacency_matrix get.adjacency list.vertex.attributes
#' vertex_attr delete.edges delete_vertices
#' @export
build_graph <-function(regulon, mode = "tripartite", weights = "corr",
                       aggregation_function = mean){
    stopifnot(mode %in% c("tg", "re", "tripartite", "pairs"))

    # give names to the peaks and target genes which will be easy to extract
    regulon$idxATAC <- paste0(as.character(regulon$idxATAC), "_peak")
    regulon$target <- paste0(regulon$target, "_gene")
    if (mode == "re")
        vertex_columns <- c("tf", "idxATAC")

    else if (mode == "pairs")
        vertex_columns <- c("tf", "idxATAC", "target")

    else
        vertex_columns <- c("tf", "target")

    graph_data <- regulon[,stats::na.omit(c(vertex_columns, weights))]
    if (mode =="tripartite"){
        # add tf-re data
        colnames(graph_data) <- stats::na.omit(c("from", "to", weights))
        graph_data_tf_re <- data.frame(from = regulon$tf, to = regulon$idxATAC)

        # weights between tf and re are set to 1 for weighted graph
        if (!is.na(weights)) graph_data_tf_re[,weights] <- 1
        graph_data <- rbind(graph_data, graph_data_tf_re)
        rm(graph_data_tf_re)
        vertex_columns <- c("from", "to")
    }
    if (mode == "pairs"){
        # create node names corresponding to re-tg pairs
        graph_data$target <- paste(graph_data$idxATAC, graph_data$target, sep ="@")
        graph_data <- graph_data[, c("tf", "target", weights)]
        vertex_columns <- c("tf", "target")
    }

    if (is.na(weights))
        graph_data <- unique(graph_data)
    else{
        colnames(graph_data)[colnames(graph_data) == weights] <- "weight"
        grouping_factors <- paste(vertex_columns, collapse="+")
        aggregation_formula <- eval(parse(text=paste0("weight~", grouping_factors)))
        graph_data <- stats::aggregate(graph_data, aggregation_formula, aggregation_function)
    }
    epiregulon_graph <- graph_from_data_frame(graph_data)
    if (mode == "tripartite"){
        layer_numb <- rep(1, vcount(epiregulon_graph))
        layer_numb[grepl("_gene$", V(epiregulon_graph)$name)] <- 2
        layer_numb[grepl("_peak$", V(epiregulon_graph)$name)] <- 3
        V(epiregulon_graph)$layer <- layer_numb
    }
    # set 'type' attribute for vertices required by bipartite graphs
    vertex_type <- rep("transcription factor", vcount(epiregulon_graph))
    vertex_type[grepl("_peak$", V(epiregulon_graph)$name)] <- "regulatory element"
    vertex_type[grepl("_gene$", V(epiregulon_graph)$name)] <- "target gene"
    V(epiregulon_graph)$type <- vertex_type
    # transform character constants to numeric values for later use by graphics functions
    V(epiregulon_graph)$type.num <- match(V(epiregulon_graph)$type,
                                                  c("transcription factor", "peak", "target gene"))

    # restore original names
    V(epiregulon_graph)$name <- gsub("_gene|_peak", "", V(epiregulon_graph)$name)
    epiregulon_graph
}

#' @rdname build_graph
#' @export
build_difference_graph <- function(graph_obj_1, graph_obj_2, weighted = TRUE,
                                   abs_diff = TRUE){
    checkmate::assertClass(graph_obj_1, "igraph")
    checkmate::assertClass(graph_obj_2, "igraph")
    if(!identical(V(graph_obj_1)$name, V(graph_obj_2)$name)) {
        stop("The nodes should be the same in both graphs")}
    transformation_function <- ifelse(abs_diff, abs, identity)
    if(weighted) {
        res <- graph_from_adjacency_matrix(transformation_function(get.adjacency(graph_obj_1, attr = "weight") -
                                                 get.adjacency(graph_obj_2, attr = "weight")),
                                           weighted = TRUE)
        # remove zero-weight edges
        res <- delete.edges(res, E(res)[E(res)$weight == 0])
   } else {
        res <- graph_from_adjacency_matrix(abs(get.adjacency(graph_obj_1) -
                                                 get.adjacency(graph_obj_2)),
                                           weighted = FALSE)
    }

    if (!identical(igraph::V(graph_obj_1)$type, igraph::V(graph_obj_2)$type)) {
        warning("Types of nodes differ between graphs. Only those from the first graph are used.")
    }
  igraph::V(res)$type <- igraph::V(graph_obj_1)$type
  igraph::V(res)$type.num <- igraph::V(graph_obj_1)$type.num

  # remove nodes with no edges
  edge_numbers <- vapply(incident_edges(res, V(res), mode = "all"), length,
                         FUN.VALUE = numeric(1))
  res <- delete_vertices(res, V(res)[edge_numbers == 0])
  res
}

#' @rdname build_graph
#' @export
add_centrality_degree <- function(graph){
    checkmate::assertClass(graph, "igraph")
    V(graph)$centrality <- strength(graph)
    graph
}

#' @rdname build_graph
#' @export
normalize_centrality <- function(graph, FUN = sqrt, weighted = TRUE){
  checkmate::assertClass(graph, "igraph")
  if (!"centrality" %in% list.vertex.attributes(graph)) stop("Vertices do not have 'centrality' attribute")
  if (!"weight" %in% list.edge.attributes(graph) & weighted) stop("Set 'weight' attribute to edges or use with 'weighted = FALSE'")

  # calculate number of edges for each node
  edge_numbers <- vapply(incident_edges(graph, V(graph), mode = "all"), length,
                         FUN.VALUE = numeric(1))

  V(graph)$centrality <- V(graph)$centrality/FUN(edge_numbers)
  graph
}

#' @rdname build_graph
#' @export
rank_tfs <- function(graph, type_attr = "type"){
    checkmate::assertClass(graph, "igraph")
    rank_df <- data.frame(tf = V(graph)$name[order(V(graph)$centrality[vertex_attr(graph, type_attr) == "transcription factor"], decreasing = TRUE)],
               centrality = sort(V(graph)$centrality[vertex_attr(graph, type_attr) == "transcription factor"], decreasing = TRUE))
    rank_df$rank <- base::rank(-rank_df$centrality)
    rank_df
}

#' Plot a graph build based on \code{getRegulon} output
#'
#' This function takes an input an igraph object created by any of the following:
#' \code{build_graph}, \code{add_centrality_degree}, \code{igraph::strength}, \code{normalize_centrality}.
#' It makes a force-directed layout plot to visualize it at a high level.
#'
#' @param graph an igraph object
#' @param layout a layout specification. Any values that are valid for
#' \link[ggraph]{ggraph} or \link[ggraph]{create_layout} will work. Defaults to
#' "stress". Consider also trying "mds", "nicely", and "fr" while you experiment.
#' @param label_size an integer indicating how large the labels of highlighted
#' transcription factors should be
#' @param tfs_to_highlight a character vector specifying which TFs in the plot
#' should be highlighted. Defaults to NULL (no labels).
#' @param edge_alpha a numeric value between 0 and 1 indicating the level of
#' transparency to use for the edge links in the force-directed layout. Defaults
#' to 0.02.
#' @param point_size a numeric value indicating the size of nodes in the force-directed layout
#' @param point_border_size a numeric value indicating the size of point
#' borders for nodes in the force-directed layout
#' @param label_alpha a numeric value between 0 and 1 indicating the level of
#' transparency to use for the labels of highlighted nodes
#' @param label_nudge_x a numeric value indicating the shift of the labels
#' along the x axis that should be used in the force-directed layout
#' @param label_nudge_y A numeric value indicating the shift of the labels
#' along the y axis that should be used in the force-directed layout.
#' @param ... optional additional arguments to pass to \link[ggraph]{create_layout}
#' @return a ggraph object
#' @author Timothy Keyes, Tomasz Wlodarczyk
#' @examples
#' # create an artificial getRegulon output
#' set.seed(1234)
#' tf_set <- apply(expand.grid(LETTERS[seq_len(5)], LETTERS[seq_len(5)]),1,  paste, collapse = "")
#' regulon <- data.frame(tf = sample(tf_set, 5e2, replace = TRUE))
#' gene_set <- expand.grid(LETTERS[seq_len(5)], LETTERS[seq_len(5)], LETTERS[seq_len(5)])
#' gene_set <- apply(gene_set,1,function(x) paste0(x,collapse=""))
#' regulon$target <- sample(gene_set, 5e2, replace = TRUE)
#' regulon$idxATAC <- seq_len(5e2)
#' regulon$corr <- runif(5e2)*0.5+0.5
#' #create igraph object
#' graph_tripartite <- build_graph(regulon, mode = "tripartite")
#' plot_epiregulon_network(graph_tripartite, tfs_to_highlight = sample(unique(tf_set),3),
#' edge_alpha = 0.2)
#' @export
plot_epiregulon_network <-
    function(
        graph,
        layout = "stress",
        label_size = 3,
        tfs_to_highlight = NULL,
        edge_alpha = 0.02,
        point_size = 1,
        point_border_size = 0.5,
        label_alpha = 0.8,
        label_nudge_x = 0.2,
        label_nudge_y = 0.2,
        ...
    ) {
        checkmate::assertClass(graph, "igraph")
        my_layout <- ggraph::create_layout(graph, layout = layout, ...)
        highlighted <- my_layout[my_layout$name  %in% tfs_to_highlight, ]
        my_plot <-
            ggraph::ggraph(graph = my_layout) +
            #ggraph::ggraph(graph = graph, layout  = my_layout) +
            ggraph::geom_edge_link(alpha = edge_alpha) +
            ggraph::geom_node_point(
                ggplot2::aes_string(fill = "type"),
                shape = 21,
                size = point_size,
                stroke = point_border_size
            ) +
            ggraph::geom_node_label(
                ggplot2::aes_string(label = "name"),
                data = highlighted,
                alpha = label_alpha,
                nudge_x = label_nudge_x,
                nudge_y = label_nudge_y,
                size = label_size
            ) +
            ggraph::geom_node_point(
                ggplot2::aes_string(fill = "type"),
                shape = 21,
                data = highlighted,
                size = 3,
                stroke = point_border_size
            ) +
            ggplot2::theme_void() +
            ggplot2::labs(fill = NULL)

        return(my_plot)
    }

#' Plot graph according to grouping factor
#'
#' Plot graph with separate weights for different levels of the grouping factor
#'
#' @param regulon an object returned by the getRegulon or addWeights function
#' @param cutoff a numeric used to select values of the variables passed in `groups`
#' parameter. Values greater than `cutoff` are retained and used as
#' graph edge weights.
#' @param tf a character vector storing the names of transcription factors to be
#' included in the graph
#' @param groups a character indicating levels of the grouping factor; should
#' correspond to columns of `regulon` to be used as `weight`
#' variable in the long format of `epiregulon` data frame.
#' @param layout a layout specification. Any values that are valid for
#' \link[ggraph]{ggraph} or \link[ggraph]{create_layout} will work.
#' @author Xiaodai Yao, Tomasz Wlodarczyk
#' @return a ggraph object
#' @examples
#' #' # create an artificial getRegulon output
#' set.seed(1234)
#' tf_set <- apply(expand.grid(LETTERS[1:10], LETTERS[1:10]),1,  paste, collapse = "")
#' regulon <- data.frame(tf = sample(tf_set, 5e3, replace = TRUE))
#' gene_set <- expand.grid(LETTERS[1:10], LETTERS[1:10], LETTERS[1:10])
#' gene_set <- apply(gene_set,1,function(x) paste0(x,collapse=""))
#' regulon$target <- sample(gene_set, 5e3, replace = TRUE)
#' regulon$idxATAC <- 1:5e3
#' regulon <- cbind(regulon, data.frame(C1 = runif(5e3), C2 = runif(5e3),
#' C3 = runif(5e3)))
#' plot_difference_network(regulon, tf = unique(tf_set)[1:3],
#' groups = c("C1", "C2", "C3"), cutoff = 0.2)
#' @export

plot_difference_network <- function(regulon,
                                    cutoff = 0.01,
                                    tf = NULL,
                                    groups  = NULL,
                                    layout = "stress"){
  regulon.tf <- list()
  for (group in groups) {
    regulon_group <- regulon[regulon$tf %in% tf, c("tf","target", group)]

    #apply cutoff
    regulon_group <- regulon_group[regulon_group[, group] > cutoff, ]

    #rename colnames as weight to be consistent across all groups
    colnames(regulon_group)[colnames(regulon_group) == group] <- "weight"

    #rename tf to be tf_group
    regulon_group[,"tf"] <- paste0(regulon_group[,"tf"], "_", group)
    regulon.tf[[group]] <- regulon_group
  }

  combined.regulon <- do.call("rbind", regulon.tf)

  combined.graph <- build_graph(combined.regulon, mode = "tg", weights = "weight")

  plot_epiregulon_network(combined.graph,
                          layout = layout,
                          tfs_to_highlight = unique(combined.regulon$tf),
                          label_nudge_x = 0.1,
                          label_nudge_y = 0.1)
}


