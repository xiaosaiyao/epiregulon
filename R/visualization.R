

plotActivityDim_ <- function(sce,
                             activity_matrix,
                             tf,
                             dimtype,
                             label,
                             legend.label,
                             colors,
                             limit,
                             ...){

  tf.activity <- as.numeric(activity_matrix[tf,])
  sce$activity <- tf.activity

  g <- scater::plotReducedDim(sce, dimred = dimtype, colour_by="activity",
                                text_by = label,...)

  g <- g + scale_color_gradient(low = colors[1], high = colors[2], limit = limit, oob=scales::squish) + ggtitle(tf) +
    labs(color=legend.label) +
    theme_classic(base_size = 12) + theme(plot.title = element_text(hjust = 0.5))

  return(g)

}


#' Plot cell-level reduced dimension results stored in a SingleCellExperiment object, colored by activities for a list of TFs
#'
#' @param sce  A SingleCellExperiment object containing dimensionality reduction coordinates
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param tf A character vector indicating the names of the transcription factors to be plotted
#' @param dimtype String indicating the name of dimensionality reduction matrix to be extracted from the SingleCellExperiment
#' @param label String corresponding to the field in the colData of sce for annotation on plot
#' @param ncol A integer to specify the number of columns in the combined plot, if combine == TRUE
#' @param nrow A integer to specify the number of rows in the combined plot, if combine == TRUE
#' @param title A string to specify the name of the combined plot
#' @param combine logical to indicate whether to combine and visualize the plots in one panel
#' @param legend.label String indicating the name of variable to be plotted on the legend
#' @param colors A vector of 2 colors for the intensity, with the first element refering to the lower value and
#' the second elment refering to the higher value. Default is c("blue","yellow").
#' @param limit A vector of lower and upper bounds for the color scale. The default option is NULL and will adjust
#' to minimal and maximal values
#' @param ... Additional arguments from scater::plotReducedDim
#'
#' @return A combined ggplot object or a list of ggplots if combine == FALSE
#' @import SingleCellExperiment
#' @export
#'
#' @examples
#' # create a mock singleCellExperiment object for gene expression matrix
#' example_sce <- scuttle::mockSCE()
#' example_sce <- scuttle::logNormCounts(example_sce)
#' example_sce <- scater::runPCA(example_sce)
#' example_sce <- scater::runUMAP(example_sce)
#' example_sce$cluster <- sample(LETTERS[1:5], ncol(example_sce), replace = TRUE)
#' plotActivityDim(sce = example_sce, activity = logcounts(example_sce),
#' tf = c("Gene_0001","Gene_0002"),  label = "cluster")

#' @author Xiaosai Yao, Shang-yang Chen
#'
plotActivityDim <- function(sce = NULL,
                            activity_matrix,
                            tf,
                            dimtype ="UMAP",
                            label = NULL,
                            ncol = NULL,
                            nrow = NULL,
                            title = NULL,
                            combine = TRUE,
                            legend.label = "activity",
                            colors = c("blue","yellow"),
                            limit = NULL,
                            ...){

  # give warning for genes absent in tf list
  missing <- tf[which(! tf %in% rownames(activity_matrix))]
  if (!identical(missing, character(0))) {
    writeLines(paste0(missing, " not found in activity matrix. Excluded from plots"))
  }

  tf <- tf[which(tf %in% rownames(activity_matrix))]


  gs <- lapply(tf, function(x) {
    suppressMessages(return(plotActivityDim_(sce, activity_matrix, x, dimtype, label, legend.label, colors, limit, ...)))
  })

  if (combine == TRUE) {

    gs <- patchwork::wrap_plots(gs, ncol = ncol, nrow=nrow) +
      patchwork::plot_annotation(title = title)

    return(gs)

  } else {

    return(gs)

  }

}



plotActivityViolin_ <- function(activity_matrix,
                                tf,
                                clusters,
                                legend.label,
                                colors){

  tf.activity <- as.numeric(activity_matrix[tf,])
  df <- data.frame(activity = tf.activity, clusters = clusters)

  g <- ggplot2::ggplot(df, aes_string(x = "clusters",
                                      y = "activity",
                                      fill = "clusters")) +
    geom_violin() +
    theme_classic(base_size = 12) +
    ggtitle(tf) + ylab(legend.label) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45,hjust = 1))

  if (!is.null(colors)){
    g <- g +  scale_fill_manual(values = colors)
  }

  return(g)

}

#' Generate violin plots of inferred activities for a list of TFs grouped by cluster/group labels
#'
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param tf A character vector indicating the names of the transcription factors to be plotted
#' @param clusters A vector of cluster or group labels for single cells
#' @param ncol A integer to indicate the number of columns in the combined plot, if combine == TRUE
#' @param combine logical to indicate whether to combine and visualize the plots in one panel
#' @param legend.label String indicating the name of variable to be plotted on the legend
#' @param colors  A character vector representing the names of colors
#'
#' @return A combined ggplot object or a list of ggplots if combine == FALSE
#' @export
#' @import ggplot2
#'
#' @examples
#' # create a mock singleCellExperiment object for gene expression matrix
#' example_sce <- scuttle::mockSCE()
#' example_sce <- scuttle::logNormCounts(example_sce)
#' example_sce$cluster <- sample(LETTERS[1:5], ncol(example_sce), replace = TRUE)
#' plotActivityViolin(activity_matrix = logcounts(example_sce),
#' tf = c("Gene_0001","Gene_0002"),  clusters = example_sce$cluster)
#'
#' @author Xiaosai Yao, Shang-yang Chen
plotActivityViolin <- function(activity_matrix,
                               tf,
                               clusters,
                               ncol = NULL,
                               combine = TRUE,
                               legend.label = "activity",
                               colors = NULL){

  # give warning for genes absent in tf list
  missing <- tf[which(! tf %in% rownames(activity_matrix))]
  if (!identical(missing, character(0))) {
    writeLines(paste0(missing, " not found in activity matrix. Excluded from plots"))
  }
  tf <- tf[which(tf %in% rownames(activity_matrix))]

  gs <- lapply(tf, function(x) {
    return(plotActivityViolin_(activity_matrix, x, clusters, legend.label, colors))
  })

  if (combine == TRUE) {

    gs <- patchwork::wrap_plots(gs, ncol = ncol)

    return(gs)

  } else {

    return(gs)

  }

}

#' Generate bubble plots of relative activities across cluster/group labels for a list of TFs
#'
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param tf A character vector indicating the names of the transcription factors to be plotted
#' @param clusters A character or integer vector of cluster or group labels for single cells
#' @param bubblesize String indicating the variable from findDifferentialActivity output to scale size of bubbles by. Default is FDR
#'
#' @return A ggplot object
#' @export
#' @import ggplot2
#'
#' @examples
#' example_sce <- scuttle::mockSCE()
#' example_sce <- scuttle::logNormCounts(example_sce)
#' example_sce$cluster <- sample(LETTERS[1:5], ncol(example_sce), replace = TRUE)
#' plotBubble(activity_matrix = logcounts(example_sce),
#' tf = c("Gene_0001","Gene_0002"),  clusters = example_sce$cluster)
#' @author Shang-yang Chen
plotBubble <- function (activity_matrix,
                        tf,
                        clusters,
                        bubblesize = "FDR"){
  # give warning for genes absent in tf list
  missing <- tf[which(! tf %in% rownames(activity_matrix))]
  if (!identical(missing, character(0))) {
    writeLines(paste0(missing, " not found in activity matrix. Excluded from plots"))
  }
  tf <- tf[which(tf %in% rownames(activity_matrix))]

  #find logFC and FDR of TFs
  markers <- findDifferentialActivity(activity_matrix, clusters,
                                      pval.type = "some", direction = "up", test.type = "t")
  markers <- suppressMessages(getSigGenes(markers, fdr_cutoff = 1.5,
                                          logFC_cutoff = -100))
  markers <- markers[which(markers$tf %in% tf), ]
  levels <- make.names(unique(tf[tf %in% markers$tf]))
  markers$tf <- make.names(markers$tf)

  # z normalize activity and compute mean by cluster
  tf.activity <- activity_matrix[tf, ,drop=FALSE]
  df <- data.frame(clusters = clusters, t(as.matrix(tf.activity)))
  df.mean <- stats::aggregate(. ~ clusters, df, mean)
  zscores <- apply(df.mean[, -1], 2, scale)
  df.mean <- data.frame(clusters = df.mean[, 1], as.data.frame(zscores))
  df.plot <- suppressMessages(reshape2::melt(df.mean, id.variable="clusters", variable.name = "tf", value.name = "relative_activity"))

  # merge logFC, FDR and mean activity
  df.plot <- merge(df.plot, markers)
  df.plot$tf <- factor(as.character(df.plot$tf), levels = levels )

  # generate bubble plots
  if (bubblesize == "FDR") {
    logpval <- -log10(df.plot$FDR)
    max.logpval <- max(logpval[is.finite(logpval)])
    logpval <- replace(logpval, is.infinite(logpval), max.logpval)
    g <- ggplot2::ggplot(df.plot, aes_string("clusters", "tf",
                                             color = "relative_activity")) +
      geom_point(stat = "identity", aes(size = logpval)) +
      scale_color_viridis_c() +
      scale_size_continuous(range = c(0, 7)) + theme_classic(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  else if (bubblesize == "summary.logFC") {
    g <- ggplot2::ggplot(df.plot, aes_string("clusters", "tf",
                                             color = "relative_activity")) + geom_point() + scale_color_viridis_c() +
      scale_size_continuous(range = c(0, 7)) + theme_classic(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  return(g)
}


enrichPlot_ <- function(results,
                        title,
                        top) {
  results$logP.adj <- -log10(results$p.adjust)
  ggplot(results[seq_len(top), ] , aes_string(y = "logP.adj",
                                       x = "Description",
                                       color = "GeneRatio")) +
    scale_colour_gradient(high = "red", low = "blue") +
    geom_point(stat = 'identity', aes_string(size = "Odds.Ratio")) +
    coord_flip() +
    theme_bw() + ggtitle (title) +
    theme(
      text = element_text(size = 10),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.5, "lines")
    )
}

#' Plot results of regulonEnrich
#'
#' @param results  Output from regulonEnrich
#' @param top An integer to indicate the number of pathways to plot ranked by significance. Default is 15.
#' @param ncol An integer to indciate the number of columns in the combined plot, if combine == TRUE. Default is 3.
#' @param combine logical to indicate whether to combine and visualize the plots in one panel
#'
#' @return A combined ggplot object or a list of ggplots if combine == FALSE
#' @export
#' @import ggplot2

#' @examples
#' #retrieve genesets
#' H <- EnrichmentBrowser::getGenesets(org = "hsa", db = "msigdb",
#' cat = "H", gene.id.type = "SYMBOL" )
#' C6 <- EnrichmentBrowser::getGenesets(org = "hsa", db = "msigdb",
#' cat = "C6", gene.id.type = "SYMBOL" )
#'
#' #combine genesets and convert genesets to be compatible with enricher
#' gs <- c(H,C6)
#' gs.list <- do.call(rbind,lapply(names(gs), function(x) {
#' data.frame(gs=x, genes=gs[[x]])}))
#'
#' head(gs.list)
#'
#' #get regulon
#' library(dorothea)
#' data(dorothea_hs, package = "dorothea")
#' regulon <- dorothea_hs
#' enrichment_results <- regulonEnrich(c("ESR1","AR"), regulon = regulon, corr = "mor",
#' genesets = gs.list)
#'
#' # plot graph
#' enrichPlot(results = enrichment_results )
#'
#' @author Xiaosai Yao
enrichPlot <- function(results,
                       top = 15,
                       ncol = 3,
                       combine = TRUE) {

  gs <- lapply(names(results), function(x) {
    enrichPlot_(results[[x]], x, top)
  })

  if (combine == TRUE) {

    gs <- patchwork::wrap_plots(gs, ncol = ncol)

    return(gs)

  } else {

    return(gs)

  }

}
