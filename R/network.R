
#' Plot networks graph of significant genesets from regulonEnrich results
#'
#' @param tf A vector of gene names to be plotted. They should be present in enrichresults
#' @param enrichresults Output from regulonEnrich that computes enriched genesets from user-specified regulons
#' of interest
#' @param ntop_pathways An integer indicating the number of top pathways to be included in the graph
#' @param p.adj_cutoff A scalar indicating the p.adjusted cutoff for pathways to be included in the graph. Default value is 0.05
#' @param layout String indicating layout option from igraph
#' @param tf_label String indicating the name of the tf label
#' @param gset_label String indicating the name of the geneset label
#' @param tf_color String indicating the color of the tf label
#' @param gset_color String indicating the color of the geneset label
#'
#' @return an igraph plot of interconnected pathways through TFs
#' @export
#' @examples
#' \dontrun{
#' AR <- data.frame(ID = c("ANDROGEN RESPONSE","PROLIFERATION","MAPK"),
#' p.adjust = c(0.001, 0.01, 0.04))
#' GATA6 <- data.frame(ID = c("STK33","PROLIFERATION","MAPK"),
#' p.adjust = c(0.001, 0.01, 0.04))
#' enrichresults <- list(AR = AR, GATA6 = GATA6)
#' plotGseaNetwork(tf = names(enrichresults), enrichresults = enrichresults)
#' }

#' @author Phoebe Guo, Xiaosai Yao

plotGseaNetwork <- function(tf,
                            enrichresults,
                            ntop_pathways = 10,
                            p.adj_cutoff = 0.05,
                            layout = "sugiyama",
                            tf_label = "tf",
                            gset_label = "ID",
                            tf_color = "tomato",
                            gset_color = "grey") {

  #filter non-significant pathways by p.adj
  enrichresults.filter <- lapply(stats::setNames(names(enrichresults),names(enrichresults)), function (x) {
    enrichresults[[x]][which(as.numeric(enrichresults[[x]][,"p.adjust"]) < p.adj_cutoff), ]
  })

  # extract top pathways from enrichr results
  pathway.list <- lapply(stats::setNames(tf,tf), function (x) {
    head(enrichresults.filter[[x]][order(as.numeric(enrichresults.filter[[x]][, "p.adjust"])), ], ntop_pathways)
  })


  #filter out empty elements
  pathway.list <- pathway.list[lapply(pathway.list, nrow)>0]


  # bind rows of pathway list
  pathway.df <- do.call(rbind, lapply(names(pathway.list), function(x) {
    data.frame(tf = x, ID = pathway.list[[x]][, gset_label])
  }))

  # create node object for igraph
  nodes.list <- lapply(names(pathway.list), function(x) {
    data.frame(rbind(data.frame(name = x, type = tf_label),
          data.frame(name = pathway.list[[x]][, gset_label], type = gset_label)))
  })

  nodes <- unique(do.call(rbind, nodes.list))


  # create the network object
  p <- igraph::graph_from_data_frame(d = pathway.df,
                             vertices = nodes,
                             directed = F)

  # Create a vector of color
  col  <- c(gset_color, tf_color)
  pal <- col[as.numeric(as.factor(igraph::V(p)$type))]
  igraph::V(p)$type.num <- as.numeric(factor(igraph::V(p)$type, levels = c(tf_label, gset_label)))

  igraph::V(p)$color <- pal

  # Create sugiyama layout
  l <- igraph::layout_with_sugiyama(p, attributes = "all",
                            layers = igraph::V(p)$type.num)

  # plot
  gseaplot <-
    ggraph::ggraph(l$extd_graph,
                   layout =  "sugiyama",
                   layers = l$layout[, 2]) +
    ggraph::geom_edge_link(alpha = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(values = pal) +
    ggraph::theme_graph(fg_text_colour = 'black') +
    ggraph::geom_node_point(aes_string(color = "type"), size = 5)  +
    ggraph::geom_node_text(aes(label = name, filter = type == gset_label),
                            nudge_y = 0.1, hjust = 0) +
    ggraph::geom_node_text(aes(label = name, filter = type == tf_label),
                           nudge_y = -0.1, hjust = 1) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(1, 4)))
  print(gseaplot)

}


