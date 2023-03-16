aggregateMatrix <- function(regulon, mode, FUN,...){
  regulon$tf <- as.factor(regulon$tf)
  regulon$target <- as.factor(regulon$target)
  groupings <- interaction(regulon$tf,regulon$target, sep = '_')
  index <- order(groupings)
  regulon <- regulon[index,]
  breaks <- which(!duplicated(groupings[index]))
  aggregated <- lapply(seq_len(length(breaks)-1), function(i){
    FUN(as.matrix(regulon[breaks[i]:(breaks[i+1]-1), mode, drop=FALSE]), ...)})

  aggregated[[length(breaks)]] <-
    FUN(as.matrix(regulon[breaks[length(breaks)]:nrow(regulon), mode, drop=FALSE]))
  aggregated <- do.call(rbind, aggregated)
  aggregated <- S4Vectors::DataFrame(tf=regulon$tf[breaks],
                                     target = regulon$target[breaks],
                                     weight = I(aggregated))
  colnames(aggregated)[3] <- mode
  aggregated
}



renameCluster <- function(clusters) {
  if (!is.null(clusters)){
    clusters[clusters == "all"] <- "clusters_all"}
  clusters
}



initiateMatCluster <- function(clusters, nrow, value=NA) {
  unique_clusters <- sort(unique(clusters))
  cluster_mat <- matrix(value, nrow = nrow, ncol = length(unique_clusters) + 1)
  colnames(cluster_mat) <- c("all", unique_clusters)
  cluster_mat
}



combineSCE <- function(expMatrix, exp_assay, peakMatrix, peak_assay, reducedDim, useDim) {

  # convert expMatrix and peakMatrix in case they weren't already so
  expMatrix <- as(expMatrix,"SingleCellExperiment")
  peakMatrix <- as(peakMatrix,"SingleCellExperiment")

  # assay of peakMatrix needs to be named as counts
  names(assays(peakMatrix)[peak_assay]) <- "counts"

  sce <- SingleCellExperiment(list(counts = assay(expMatrix, exp_assay)),
                              altExps = list(peakMatrix = peakMatrix))

  rowRanges(sce) <- rowRanges(expMatrix)
  rowRanges(altExp(sce)) <- rowRanges(peakMatrix)

  # add reduced dimension information to sce object
  reducedDim(sce, useDim) <- reducedDim

  sce

}

#' @import ggplot2
plotDiagnostic <- function(idx,regulon, expMatrix,exp_assay, exp_cutoff=1, peakMatrix, peak_assay, peak_cutoff=0, clusters){

  target <- regulon$target[idx]
  tf <- regulon$tf[idx]
  peak <- regulon$idxATAC[idx]

  tg_plot <- list()
  for (cluster in unique(clusters)){
    target_exp <- SummarizedExperiment::assay(expMatrix, exp_assay)[target, clusters == cluster]
    tf_exp <- SummarizedExperiment::assay(expMatrix, exp_assay)[tf, clusters == cluster]
    peak_accessibility <- SummarizedExperiment::assay(peakMatrix, peak_assay)[peak, clusters == cluster]

    if (is.null(exp_cutoff)){
      exp_cutoff <- mean(tf_exp)
    }

    if (is.null(peak_cutoff)){
      peak_cutoff <- mean(peak_accessibility)
    }
    tf_exp.bi <- binarize(tf_exp, cutoff=exp_cutoff)
    peak.bi <- binarize(peak_accessibility, cutoff=peak_cutoff)
    tf_re.bi <- tf_exp.bi*peak.bi
    tf_re.bi <- factor(tf_re.bi, levels = c(0,1))

    require(ggplot2)
    require(ggbeeswarm)
    target_group <- data.frame(groups=as.vector(tf_re.bi), target=as.vector(target_exp))
    tg_plot[[cluster]] <- ggplot(target_group, aes(groups, target)) +  geom_boxplot() + geom_point(position = "jitter") +
      ggtitle(paste0("target:", target, " tf:", tf, " peak:", peak, " in ", cluster, "\n",
                      " GRN corr:", round(regulon$corr[idx,cluster],2),
                      " GRN pval:", round(regulon$pval[idx,cluster],2),
                      " weight:", round(regulon$weight[idx,cluster],2)))

  }

  tg_plot


}

binarize <- function(input_vector, cutoff) {
  filtered <- as.numeric(input_vector > cutoff)
}
