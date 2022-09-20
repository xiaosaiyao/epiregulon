#' Building igraph directed graph object based on the output of the getRegulon function
#'
#' @param regulon An object returned by the getRegulon or addWeights function
#' @param mode A character specifying whch type of graph will be built. In 'tg' mode
#' a bipartite graph is built connecting transcription factors directly to target genes
#' and ignoring information about mediating regulatory elements; in 'pairs' mode
#' transcription factors are connected to unique target gene-regulatory element pairs;
#' in tripartite mode the network is build of three type of vertices (nodes):
#' trascription factors, regulatory elements and target genes; here the path from
#' target gene to regulatory element always contains a regulatory elements
#' @param weights A character specifying which variable should be used to assign
#' weights to edges. If set to 'NA' then unweighted graph is built.
#' @importFrom igraph graph_from_data_frame
#'

build_graph <-function(regulon, mode = "tripartite", weights = "corr"){
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

    graph_data <- regulon[,na.omit(c(vertex_columns, weights))]
    if (mode =="tripartite"){
        # add tf-re data
        colnames(graph_data) <- na.omit(c("from", "to", weights))
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
        graph_data <- aggregate(graph_data, aggregation_formula, mean)
    }
    epiregulon_graph <- igraph::graph_from_data_frame(graph_data)
    if (mode == "tripartite"){
        layer_numb <- rep(1, igraph::vcount(epiregulon_graph))
        layer_numb[grepl("_gene$", igraph::V(epiregulon_graph)$name)] <- 2
        layer_numb[grepl("_peak$", igraph::V(epiregulon_graph)$name)] <- 3
        igraph::V(epiregulon_graph)$layer <- layer_numb
    }
    # set 'type' attribute for vertices required by bipartite graphs
    vertex_type <- rep("transcription factor", igraph::vcount(epiregulon_graph))
    vertex_type[grepl("_peak$", igraph::V(epiregulon_graph)$name)] <- "regulatory element"
    vertex_type[grepl("_gene$", igraph::V(epiregulon_graph)$name)] <- "target gene"
    igraph::V(epiregulon_graph)$type <- vertex_type
    # transform character constants to numeric values for later use by graphics functions
    igraph::V(epiregulon_graph)$type.num <- match(igraph::V(epiregulon_graph)$type,
                                                  c("transcription factor", "peak", "target gene"))

    # restore original names
    igraph::V(epiregulon_graph)$name <- gsub("_gene|_peak", "", igraph::V(epiregulon_graph)$name)
    epiregulon_graph
}



build_difference_graph <- function(graph_obj_1, graph_obj_2, weighted = TRUE){
    res <- igraph::graph_from_adjacency_matrix(abs(igraph::get.adjacency(graph_obj_1, attr = "weight") -
                                                       igraph::get.adjacency(graph_obj_2, attr = "weight")), weighted = weighted)
    V(res)$type <- V(graph_obj_1)$type
    V(res)$type.num <- V(graph_obj_1)$type.num
    res
}

calculate_centrality_degree <- function(graph){
    igraph::V(graph)$centrality <- igraph::strength(graph)
    graph
}

rank_tfs <- function(graph){
    data.frame(tf = V(graph)$name[order(V(graph)$centrality, decreasing = TRUE)],
               rank = 1:vcount(graph),
               centrality = sort(V(graph)$centrality, decreasing = TRUE))
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
        highlighted <- my_layout[my_layout$label_col  %in% tfs_to_label, ]
        my_plot <-
            ggraph::ggraph(graph = my_layout) +
            ggraph::geom_edge_link(alpha = 0.02) +
            ggraph::geom_node_point(
                ggplot2::aes(fill = type),
                shape = 21,
                size = re_and_tg_point_size,
                stroke = re_and_tg_point_border_size
            ) +
            ggraph::geom_node_label(
                ggplot2::aes(label = name),
                data = highlighted,
                alpha = label_alpha,
                nudge_x = label_nudge_x,
                nudge_y = label_nudge_y,
                size = label_size
            ) +
            ggraph::geom_node_point(
                ggplot2::aes(fill = type),
                shape = 21,
                data = highlighted,
                size = 3,
                stroke = tf_point_border_size
            ) +
            ggplot2::theme_void() +
            ggplot2::labs(fill = NULL)

        return(my_plot)
    }

