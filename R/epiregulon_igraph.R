#' Building igraph directed graph object based on the output of the \code{getRegulon} function
#'
#' @param regulon An object returned by the getRegulon or addWeights function
#' @param mode A character specifying whch type of graph will be built. In \code{'tg'} mode
#' a bipartite graph is built connecting transcription factors directly to the target genes
#' and ignoring information about mediating regulatory elements; in \code{'pairs'} mode
#' transcription factors are connected to unique target gene-regulatory element pairs;
#' in \code{'tripartite'} mode the network is build of three type of vertices (nodes):
#' trascription factors, regulatory elements and target genes; here the path from
#' target gene to regulatory element always contains a regulatory element; in
#' \code{'re'} mode data in the target genes is dropped and only connections are
#' between transcription factors and regulatory elements.
#' @param weights A character specifying which variable should be used to assign
#' weights to edges. If set to 'NA' then unweighted graph is built.
#' @return A regulatory network graph
#' @importFrom igraph graph_from_data_frame V vcount
#' @return An igraph object
#' @examples
#' # create an artificial getRegulon output
#' set.seed(1234)
#' tf_set <- apply(expand.grid(LETTERS[1:10], LETTERS[1:10]),1,  paste, collapse = "")
#' regulon <- data.frame(tf = sample(tf_set, 5e3, replace = TRUE))
#' gene_set <- expand.grid(LETTERS[1:10], LETTERS[1:10], LETTERS[1:10])
#' regulon$target <- sample(gene_set, 5e3, replace = TRUE)
#' regulon$idxATAC <- 1:5e3
#' regulon$corr <- runif(5e3)*0.5+0.5
#' # build bipartite graph using regulatory element-target gena pairs
#' gr1 <- build_graph(regulon, mode = "pairs")
#' # bulid tripartite graph
#' gr2 <- build_graph(regulon, mode = "tripartite")

build_graph <- function(regulon, mode = "tripartite", weights = "corr"){
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
        graph_data <- stats::aggregate(graph_data, aggregation_formula, mean)
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

#' Build a graph difference
#'
#' Build a graph based on the edge difference of two input graphs
#'
#' Function building a graph difference by subtracting the edges of \code{graph_obj_2}
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
#' @param graph_obj_1 an igraph object from which \code{graph_obj_2} will be
#' subtracted
#' @param graph_obj_2 an igraph object used for being subtracted from \code(graph_obj_1)
#' @param weighted a logical indicating whether weighted graphs are used
#' @param abs_diff a logical indicating whether absulute difference in the number
#' edges or their weights will be calculated
#' @return an igraph object
#' @importFrom igraph get.adjacency V graph_from_adjacency_matrix
build_difference_graph <- function(graph_obj_1, graph_obj_2, weighted = TRUE,
                                   abs_diff = TRUE){
    if(!identical(V(graph_obj_1)$name, V(graph_obj_2)$name)) {
        stop("The nodes should be the same in both graphs")}
    transformation_function <- ifelse(abs_diff, abs, identity)
    if(weighted) {
        res <- graph_from_adjacency_matrix(transformation_function(get.adjacency(graph_obj_1, attr = "weight") -
                                                 get.adjacency(graph_obj_2, attr = "weight")),
                                           weighted = TRUE)
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
    res
}

#' Calculate degree centrality
#'
#' Calculate degree centrality for each vertex
#'
#' @param graph an igraph object
#' @return an igraph object with attribute \code{centrality} added to vertices
#' @importFrom igraph V
add_centrality_degree <- function(graph){
    V(graph)$centrality <- strength(graph)
    graph
}


#' Rank transcription factors
#'
#' Rank transcription factors according to degree centrality of their vertices
#'
#' @param graph an igraph object with \code{centrality} attribute added to vertices
#' @return a data.frame with transcription factors sorted according to the value of the
#' \code{centrality} attribute
#' @importFrom igraph V vcount
rank_tfs <- function(graph){
    rank_df <- data.frame(tf = V(graph)$name[order(V(graph)$centrality, decreasing = TRUE)],
               centrality = sort(V(graph)$centrality, decreasing = TRUE))
    rank_df$rank <- base::rank(-rank_df$centrality)
    rank_df
}

plot_epiregulon_network <-
    function(
        graph,
        layout = "stress",
        label_size = 3,
        tfs_to_label = NULL,
        edge_alpha = 0.02,
        tf_point_size = 3,
        re_and_tg_point_size = 1,
        tf_point_border_size = 0.5,
        re_and_tg_point_border_size = 0.5,
        label_alpha = 0.8,
        label_nudge_x = 0.5,
        label_nudge_y = 0.5,
        ...
    ) {
        my_layout <- ggraph::create_layout(graph, layout = layout)
        highlighted <- my_layout[my_layout$name  %in% tfs_to_label, ]
        my_plot <-
            ggraph::ggraph(graph = my_layout) +
            #ggraph::ggraph(graph = graph, layout  = my_layout) +
            ggraph::geom_edge_link(alpha = 0.02, aes_string(color = "from")) +
            ggraph::geom_node_point(
                ggplot2::aes_string(fill = "type"),
                shape = 21,
                size = re_and_tg_point_size,
                stroke = re_and_tg_point_border_size
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
                stroke = tf_point_border_size
            ) +
            ggplot2::theme_void() +
            ggplot2::labs(fill = NULL)

        return(my_plot)
    }

plot_difference_network <- function(regulon,
                                    cutoff = 0.01,
                                    tf = NULL,
                                    groups  = NULL,
                                    layout = "stress"){
  regulon.tf <- list()
  for (group in groups) {
    regulon.tf[[group]] <- regulon[which(regulon$tf %in% tf), c("tf","target", group)]

    #apply cutoff
    regulon.tf[[group]] <- regulon.tf[[group]][which(regulon.tf[[group]][, group] > cutoff),]

    #rename colnames as weight to be consistent across all groups
    colnames(regulon.tf[[group]])[which(colnames(regulon.tf[[group]]) == group)] <- "weight"

    #rename tf to be tf_group
    regulon.tf[[group]][,"tf"] <- paste0(regulon.tf[[group]][,"tf"], "_", group)
  }

  combined.regulon <- do.call("rbind", regulon.tf)

  combined.graph <- build_graph(combined.regulon, mode = "tg", weights = "weight")


  plot_epiregulon_network(combined.graph,
                          layout = layout,
                          tfs_to_label = unique(combined.regulon$tf),
                          label_nudge_x = 0.1,
                          label_nudge_y = 0.1)
}
