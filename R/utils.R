
#' Convert ArchR Gene Expression Matrix to SingleCellExperiment
#'
#' @param rse a RangedSummarizedExperiment object obtained from `ArchR::getMatrixFromProject`
#'
#' @return A SingleCellExperiment object
#' @export
#' @details When gene expression matrices (either GeneIntegrationMatrix or GeneExpressionMatrix)
#' are retrieved from ArchR project, they are returned as `RangedSummarizedExperiment`. This function
#' converts the class of the object from `RangedSummarizedExperiment` to `SingleCellExperiment`.
#' During the conversion, the genomic location of the genes are transferred from `rowData` in the
#' `RangedSummarizedExperiment` to the `rowRanges` of the `SingleCellExperiment`.
#'
ArchRMatrix2SCE <- function(rse) {

  rowData_temp <- rowData(rse)
  rse <- as(rse, "SingleCellExperiment")
  rowData_temp$strand[rowData_temp$strand == 1] = "+"
  rowData_temp$strand[rowData_temp$strand == 2] = "-"
  rowData_temp$strand[rowData_temp$strand == 3] = "*"

  # change order
  swap_rows <- which(rowData_temp$start > rowData_temp$end)
  rowData_temp_start <- rowData_temp$start
  rowData_temp$start[swap_rows] <- rowData_temp$end[swap_rows]
  rowData_temp$end[swap_rows] <- rowData_temp_start[swap_rows]

  rowRanges(rse) <- makeGRangesFromDataFrame(rowData_temp)
  rowData(rse) <- rowData_temp[, -which(colnames(rowData_temp) %in% c("seqnames","start","end","strand"))]
  count_names <- switch(names(assays(rse)),
                        GeneExpressionMatrix = "logcounts",
                        GeneIntegrationMatrix = "logcounts",
                        PeakMatrix = "counts")
  names(assays(rse))[1] <- count_names
  rse
}

aggregateMatrix <- function(regulon, mode, FUN){

  groupings <- paste(regulon$tf,regulon$target, sep = '#')
  groupings <- factor(groupings, levels = unique(groupings))

  weights <- regulon[, mode]
  agg.weights <- rowsum(weights, groupings, reorder = FALSE)

  if (FUN == "mean") {
    num <- table(groupings)
    agg.weights <- agg.weights / as.integer(num[rownames(agg.weights)])
  }


  rownames.split <- do.call(rbind,strsplit(rownames(agg.weights), "#"))
  aggregated <- S4Vectors::DataFrame(tf=rownames.split[,1],
                                     target = rownames.split[,2],
                                     weight = I(agg.weights))
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





#' @import ggplot2 ggbeeswarm
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
