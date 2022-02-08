#' A function to plot cell-level reduced dimension results stored in a SingleCellExperiment object, colored by activities for a specific TF
#'
#' @param sce A SingleCellExperiment object containing dimensional reduction coordinates (UMAP or tSNE)
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param tf The name of the transcription factor to be plotted
#' @param dimtype Name of dimensional reduction to be plotted (UMAP or tSNE)
#' @param label A boolean value to determine whether the cluser/group labels should be annotated on plot
#'
#' @return A ggplot object
#' @export
#'
#' @examples 1+1 = 2
plotActivityDim <- function(sce, activity_matrix, tf, dimtype="UMAP", label = NULL,...){

  tf.activity <- as.numeric(subset(activity_matrix, rownames(activity_matrix)==tf))
  sce$activity <- tf.activity

  if (!is.null(label)){


    g <- scater::plotReducedDim(sce, dimred = dimtype, colour_by="activity",
                                text_by = label, text_size = 3, text_colour = "black",...)


  } else {

    g <- scater::plotReducedDim(sce, dimred = dimtype, colour_by="activity",
                                text_by = label, text_size = 3, text_colour = "black",...)

  }

  return(g + scale_color_gradient(low="blue", high="yellow") + ggtitle(tf) +
           theme_classic(base_size = 12) + theme(plot.title = element_text(hjust = 0.5)))

}

#' A function to draw a violin plot of inferred activities for a specific TFs grouped by cluster/group labels
#'
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param tf The name of the transcription factor to be plotted
#' @param class A vector of cluster or group labels for single cells
#'
#' @return A ggplot object
#' @export
#'
#' @examples 1+1 = 2
plotActivityViolin <- function(activity_matrix, tf, class){

  tf.activity <- as.numeric(subset(activity_matrix, rownames(activity_matrix)==tf))
  df <- data.frame(activity = tf.activity, class = class)

  g <- ggplot2::ggplot(df, aes(class, activity, fill=class)) + geom_violin() + theme_classic(base_size = 12) +
    ggtitle(tf) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
                        axis.text.x=element_text(angle=45,hjust=1))

  return(g)

}

#' A function to draw bubble plot of relative activities across cluster/group labels for a list of TFs
#'
#' @param activity_matrix The name of the transcription factor to be plotted
#' @param tf.list The list of the transcription factors to be plotted
#' @param class A vector of cluster or group labels for single cells
#'
#' @return A ggplot object
#' @export
#'
#' @examples 1+1 = 2
plotBubble <- function(activity_matrix, tf.list, class){

  tf.activity <- t(subset(activity_matrix, rownames(activity_matrix) %in% tf.list))
  df <- data.frame(class = class, tf.activity)

  markers <- findDifferentialActivity(score.combine, class, pval.type="some", direction="up", test.type= "t")
  markers <- getSigGenes(markers, fdr_cutoff = 1.5)
  markers <- subset(markers, tf %in% tf.list)

  ### aggregate by class and calculate z scores
  df.mean <- aggregate(.~class, df, mean)
  zscores <- apply(df.mean[,-1], 2, scale)
  df.mean <- data.frame(class = df.mean[,1], as.data.frame(zscores))

  ### plot bubble
  df.plot <- tidyr::pivot_longer(df.mean, -c(.data$class), names_to = "tf", values_to = "relative_activity")
  df.plot <- merge(df.plot, markers)
  g <- ggplot2::ggplot(df.plot, aes(class, tf, color=relative_activity, size=summary.logFC)) +
    geom_point() + scale_color_viridis_c()+ scale_size_continuous(range=c(0,7))+
    theme_classic(base_size = 12) + theme(axis.text.x=element_text(angle=45,hjust=1))

  return(g)

}
