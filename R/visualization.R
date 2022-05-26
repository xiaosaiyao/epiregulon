#' A function to plot cell-level reduced dimension results stored in a SingleCellExperiment object, colored by activities for a specific TF
#'
#' @param sce A SingleCellExperiment object containing dimensionality reduction coordinates
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param tf A character vector indicating the names of the transcription factors to be plotted
#' @param dimtype String indicating the name of dimensionality reduction matrix to be extracted from the SingleCellExperiment
#' @param label logical to determine whether the cluster/group labels should be annotated on plot
#' @param legend.label String indicating the name of variable to be plotted on the legend
#' @param colors A vector of 2 colors for the intensity, with the first element refering to the lower value and
#' the second elment refering to the higher value. Default is c("blue","yellow").
#' @param limit A vector of lower and upper bounds for the color scale. The default option is NULL and will adjust
#' to minimal and maximal values
#' @param ... Additional arguments from scater::plotReducedDim
#'
#' @return A ggplot object
#' @export
plotActivityDim_ <- function(sce, activity_matrix, tf, dimtype, label, legend.label, colors, limit, ...){

  tf.activity <- as.numeric(activity_matrix[tf,])
  sce$activity <- tf.activity

  g <- scater::plotReducedDim(sce, dimred = dimtype, colour_by="activity",
                                text_by = label, text_size = 3, text_colour = "black",...)

  g <- g + scale_color_gradient(low=colors[1], high=colors[2], limit = limit, oob=scales::squish) + ggtitle(tf) +
    labs(color=legend.label) +
    theme_classic(base_size = 12) + theme(plot.title = element_text(hjust = 0.5))

  return(g)

}


#' A function to plot cell-level reduced dimension results stored in a SingleCellExperiment object, colored by activities for a list of TFs
#'
#' @param sce  A SingleCellExperiment object containing dimensionality reduction coordinates
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param tf A character vector indicating the names of the transcription factors to be plotted
#' @param dimtype String indicating the name of dimensionality reduction matrix to be extracted from the SingleCellExperiment
#' @param label logical to determine whether the cluster/group labels should be annotated on plot
#' @param ncol A integer to specify the number of columns in the combined plot, if combine == TRUE
#' @param combine logical to indicate whether to combine and visualize the plots in one panel
#' @param legend.label String indicating the name of variable to be plotted on the legend
#' @param colors A vector of 2 colors for the intensity, with the first element refering to the lower value and
#' the second elment refering to the higher value. Default is c("blue","yellow").
#' @param limit A vector of lower and upper bounds for the color scale. The default option is NULL and will adjust
#' to minimal and maximal values
#' @param ... Additional arguments from scater::plotReducedDim
#'
#' @return A combined ggplot object or a list of ggplots if combine == FALSE
#' @export
#'
#' @examples
#' \dontrun{
#' plotActivityDim(sce, score.combine, c("FOXA1","GATA3","SOX9","SPI1"), "TSNE", point_size = 0.25)
#'}
#'
#'
plotActivityDim <- function(sce, activity_matrix, tf, dimtype="UMAP", label = NULL, ncol = NULL, combine = TRUE,
                            legend.label = "activity", colors = c("blue","yellow"), limit=NULL, ...){

  # give warning for genes absent in tf list
  missing = tf[which(! tf %in% rownames(activity_matrix))]
  if (!identical(missing, character(0))) {
    message(paste0(missing, " not found in activity matrix. Excluded from plots"))
  }

  tf = tf[which(tf %in% rownames(activity_matrix))]


  gs <- lapply(tf, function(x) {
    suppressMessages(return(plotActivityDim_(sce, activity_matrix, x, dimtype, label, legend.label, colors, limit, ...)))
  })

  if (combine == TRUE) {

    gs <- patchwork::wrap_plots(gs, ncol = ncol)

    return(gs)

  } else {

    return(gs)

  }

}


#' A function to draw a violin plot of inferred activities for a specific TF grouped by cluster/group labels
#'
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param tf A character vector indicating the names of the transcription factors to be plotted
#' @param class A character or integer vector of cluster or group labels for single cells
#' @param legend.label String indicating the name of variable to be plotted on the legend
#'
#' @return A ggplot object
#' @export
plotActivityViolin_ <- function(activity_matrix, tf, class, legend.label){

  tf.activity <- as.numeric(activity_matrix[tf,])
  df <- data.frame(activity = tf.activity, class = class)

  g <- ggplot2::ggplot(df, aes(class, activity, fill=class)) + geom_violin() + theme_classic(base_size = 12) +
    ggtitle(tf) + ylab(legend.label) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
                        axis.text.x=element_text(angle=45,hjust=1))

  return(g)

}

#' A function to draw violin plots of inferred activities for a list of TFs grouped by cluster/group labels
#'
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param tf A character vector indicating the names of the transcription factors to be plotted
#' @param class A vector of cluster or group labels for single cells
#' @param ncol A integer to indicate the number of columns in the combined plot, if combine == TRUE
#' @param combine logical to indicate whether to combine and visualize the plots in one panel
#' @param legend.label String indicating the name of variable to be plotted on the legend
#'
#' @return A combined ggplot object or a list of ggplots if combine == FALSE
#' @export
#'
#' @examples
#' \dontrun{
#' plotActivityViolin(score.combine, c("FOXA1","GATA3","SOX9","SPI1"), sce$BioClassification)
#' }

plotActivityViolin <- function(activity_matrix, tf, class, ncol = NULL, combine = TRUE, legend.label = "activity"){
  # give warning for genes absent in tf list
  missing = tf[which(! tf %in% rownames(activity_matrix))]
  if (!identical(missing, character(0))) {
    message(paste0(missing, " not found in activity matrix. Excluded from plots"))
  }
  tf = tf[which(tf %in% rownames(activity_matrix))]

  gs <- lapply(tf, function(x) {
    return(plotActivityViolin_(activity_matrix, x, class, legend.label))
  })

  if (combine == TRUE) {

    gs <- patchwork::wrap_plots(gs, ncol = ncol)

    return(gs)

  } else {

    return(gs)

  }

}

#' A function to draw bubble plot of relative activities across cluster/group labels for a list of TFs
#'
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param tf A character vector indicating the names of the transcription factors to be plotted
#' @param class A character or integer vector of cluster or group labels for single cells
#' @param bubblesize String indicating the variable from findDifferentialActivity output to scale size of bubbles by. Default is FDR
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' plotBubble(score.combine, c("SPI1", "ZNF623", "IRF4","SOX4"), sce$BioClassification)
#'}

plotBubble <- function (activity_matrix, tf, class, bubblesize = "FDR"){
  # give warning for genes absent in tf list
  missing = tf[which(! tf %in% rownames(activity_matrix))]
  if (!identical(missing, character(0))) {
    message(paste0(missing, " not found in activity matrix. Excluded from plots"))
  }
  tf = tf[which(tf %in% rownames(activity_matrix))]

  #find logFC and FDR of TFs
  markers <- findDifferentialActivity(activity_matrix, class,
                                      pval.type = "some", direction = "up", test.type = "t")
  markers <- suppressMessages(getSigGenes(markers, fdr_cutoff = 1.5,
                                          logFC_cutoff = -100))
  markers <- markers[which(markers$tf %in% tf), ]
  levels=make.names(unique(tf[tf %in% markers$tf]))
  markers$tf = make.names(markers$tf)

  # z normalize activity and compute mean by cluster
  tf.activity <- activity_matrix[tf, ,drop=FALSE]
  df <- data.frame(class = class, t(as.matrix(tf.activity)))
  df.mean <- aggregate(. ~ class, df, mean)
  zscores <- apply(df.mean[, -1], 2, scale)
  df.mean <- data.frame(class = df.mean[, 1], as.data.frame(zscores))
  df.plot <- suppressMessages(reshape2::melt(df.mean, id.variable="class", variable.name = "tf", value.name = "relative_activity"))

  # merge logFC, FDR and mean activity
  df.plot <- merge(df.plot, markers)
  df.plot$tf = factor(as.character(df.plot$tf), levels = levels )

  # generate bubble plots
  if (bubblesize == "FDR") {
    logpval <- -log10(df.plot$FDR)
    max.logpval <- max(logpval[is.finite(logpval)])
    logpval <- replace(logpval, is.infinite(logpval), max.logpval)
    g <- ggplot2::ggplot(df.plot, aes_string("class", "tf",
                                             color = "relative_activity")) + geom_point(stat = "identity",
                                                                                        aes(size = logpval)) + scale_color_viridis_c() +
      scale_size_continuous(range = c(0, 7)) + theme_classic(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  else if (bubblesize == "summary.logFC") {
    g <- ggplot2::ggplot(df.plot, aes_string("class", "tf",
                                             color = "relative_activity")) + geom_point() + scale_color_viridis_c() +
      scale_size_continuous(range = c(0, 7)) + theme_classic(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  return(g)
}


enrichPlot_ <- function(results, title, top) {
  ggplot(results[1:top, ] , aes(y = -log10(p.adjust), x = Description, color = GeneRatio)) +
    scale_colour_gradient(high = "red", low = "blue") +
    geom_point(stat = 'identity', aes(size = Odds.Ratio)) +
    coord_flip() +
    theme_bw() + ggtitle (title) +
    theme(
      text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.margin = unit(0.5, "lines")
    )
}

#' A function to plot results of regulonEnrich
#'
#' @param results  Output from regulonEnrich
#' @param top An integer to indicate the number of pathways to plot ranked by significance. Default is 15.
#' @param ncol An integer to indciate the number of columns in the combined plot, if combine == TRUE. Default is 3.
#' @param combine logical to indicate whether to combine and visualize the plots in one panel
#'
#' @return A combined ggplot object or a list of ggplots if combine == FALSE
#' @export
#'
#' @examples
#' \dontrun{
#' enrichplot(results=enrichment_results )
#' }
#'

enrichPlot <- function(results, top = 15, ncol=3, combine=TRUE) {

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
