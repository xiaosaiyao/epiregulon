
#' plot networks graph of genesets from regulonEnrich results
#'
#' @param tf A vector of gene names to be plotted. They should be present in enrichresults
#' @param enrichresults Output from regulonEnrich that computes enriched genesets from user-specified regulons
#' of interest
#' @param ntop_pathways An integer indicating the number of top pathways to be included in the graph
#' @param p.adj_cutoff p.adjusted cutoff for pathways to be included in the graph. Default value is 0.05
#' @param layout layout option from igraph
#'
#' @return an igraph plot of interconnected pathways through TFs
#' @import igraph
#' @export
#'
#' @examples
#' \dontrun{plotGseaNetwork(tf = names(enrichresults), enrichresults, p.adj_cutoff = 0.1,
#' ntop_pathways = 10)}

plotGseaNetwork <- function(tf, enrichresults, ntop_pathways = 10, p.adj_cutoff=0.05,
                            layout = layout.fruchterman.reingold) {

  #filter non-significant pathways by p.adj
  enrichresults.filter <- lapply(setNames(names(enrichresults),names(enrichresults)), function (x) {
    enrichresults[[x]][which(as.numeric(enrichresults[[x]][,"p.adjust"]) < p.adj_cutoff), ]
  })

  # extract top pathways from enrichr results
  pathway.list <- lapply(setNames(tf,tf), function (x) {
    head(enrichresults.filter[[x]][order(as.numeric(enrichresults.filter[[x]][, "p.adjust"])), ], ntop_pathways)
  })


  #filter out empty elements
  pathway.list <- pathway.list[lapply(pathway.list, nrow)>0]


  # bind rows of pathway list
  pathway.df = do.call(rbind, lapply(names(pathway.list), function(x) {
    data.frame(tf = x, ID = pathway.list[[x]][, "ID"])
  }))

  # create node object for igraph
  nodes.list = lapply(names(pathway.list), function(x) {
    data.frame(rbind(data.frame(name = x, type = "tf"),
          data.frame(name = pathway.list[[x]][, "ID"], type = "ID")))
  })

  nodes = unique(do.call(rbind, nodes.list))


  # create the network object
  p <- igraph::graph_from_data_frame(d = pathway.df,
                             vertices = nodes,
                             directed = F)

  # Create a vector of color
  col  <- c("grey", "tomato")
  pal <- col[as.numeric(as.factor(V(p)$type))]

  plot(
    p,
    vertex.size = 15,
    vertex.color = pal,
    vertex.frame.color = "white",
    vertex.label.cex = 0.5,
    vertex.label.color = "black",
    layout = layout
  )
}
