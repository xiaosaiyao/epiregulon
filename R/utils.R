aggregateMatrix <- function(regulon, mode, FUN) {

    groupings <- paste(regulon$tf, regulon$target, sep = "#")
    groupings <- factor(groupings, levels = unique(groupings))

    weights <- regulon[, mode]
    agg.weights <- rowsum(weights, groupings, reorder = FALSE)

    if (FUN == "mean") {
        num <- table(groupings)
        agg.weights <- agg.weights/as.integer(num[rownames(agg.weights)])
    }


    rownames.split <- do.call(rbind, strsplit(rownames(agg.weights),
        "#"))
    aggregated <- S4Vectors::DataFrame(tf = rownames.split[, 1],
        target = rownames.split[, 2], weight = I(agg.weights))
    colnames(aggregated)[3] <- mode
    aggregated
}

renameCluster <- function(clusters) {
    if (!is.null(clusters)) {
        clusters[clusters == "all"] <- "clusters_all"
    }
    clusters
}

initiateMatCluster <- function(clusters, nrow, value = NA) {
    unique_clusters <- sort(unique(clusters))
    cluster_mat <- matrix(value, nrow = nrow, ncol = length(unique_clusters) +
        1)
    colnames(cluster_mat) <- c("all", unique_clusters)
    cluster_mat
}


#' @importFrom SummarizedExperiment colData<-
#' @importFrom SingleCellExperiment applySCE
.aggregateCells <- function(cellNum, expMatrix, peakMatrix, caller_env, useDim,
                            exp_assay, peak_assay, BPPARAM, clusters=NULL){
  message("performing pseudobulk using an average of ", cellNum, " cells")
  barcodes <- list()
  kclusters <- list()
  if (!is.null(clusters)) {

    #add cluster info to expMatrix
    colData(expMatrix)[, "cluster_for_pseudobulk"] <- clusters

    # K-means clustering
    kclusters <- list()
    for (cluster in unique(clusters)) {
      sce <- expMatrix[, which(clusters == cluster)]
      kNum <- trunc(ncol(sce)/cellNum)
      kclusters[[cluster]] <- scran::clusterCells(sce, use.dimred = useDim,
                                                  BLUSPARAM = bluster::KmeansParam(centers = kNum, iter.max = 5000))
      barcodes[[cluster]] <- names(kclusters[[cluster]])
      kclusters[[cluster]] <- paste(cluster, kclusters[[cluster]], sep = "_")
    }
    kclusters <- unlist(kclusters)
    barcodes <- unlist(barcodes)
    names(kclusters) <- barcodes


  } else {
    kclusters <- scran::clusterCells(expMatrix, use.dimred = useDim, BLUSPARAM = bluster::KmeansParam(centers = trunc(ncol(peakMatrix)/cellNum),
                                                                                                      iter.max = 5000))
  }

  kclusters <- kclusters[colnames(expMatrix)]

  #replace clusters with clusters of pseudobulked samples


  expMatrix <- applySCE(expMatrix, scuttle::aggregateAcrossCells, WHICH = NULL,
                        ids = kclusters, statistics = "sum", use.assay.type = exp_assay, BPPARAM = BPPARAM)

  peakMatrix <- applySCE(peakMatrix, scuttle::aggregateAcrossCells, WHICH = NULL,
                         ids = kclusters, statistics = "sum", use.assay.type = peak_assay, BPPARAM = BPPARAM)
  if (!is.null(clusters))
    clusters <- colData(expMatrix)[, "cluster_for_pseudobulk"]
  assign("expMatrix",expMatrix,envir = caller_env)
  assign("peakMatrix",peakMatrix,envir = caller_env)
  assign("clusters",clusters,envir = caller_env)
}



