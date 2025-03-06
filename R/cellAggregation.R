#' Aggregate expression across cells
#'
#' Aggregate expression values across cells based on one or more grouping factors.
#' This is primarily used to create pseudo-bulk profiles for each cluster/sample combination.
#'
#' @param x Any matrix-like object.
#' Expression values are typically expected to be counts.
#' @param factors A list or data frame containing one or more grouping factors.
#' Each entry should be a factor of the same length as the number of cells in \code{x}.
#' @param num.threads Integer specifying the number of threads to be used for aggregation.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{sums}, a numeric matrix where each row corresponds to a gene and each column corresponds to a unique combination of grouping levels.
#' Each entry contains the summed expression across all cells with that combination.
#' \item \code{detected}, an integer matrix where each row corresponds to a gene and each column corresponds to a unique combination of grouping levels.
#' Each entry contains the number of cells with detected expression in that combination.
#' \item \code{combinations}, a data frame describing the levels for each unique combination.
#' Rows of this data frame correspond to columns of \code{sums} and \code{detected},
#' while columns correspond to the factors in \code{factors}.
#' \item \code{counts}, the number of cells associated with each combination.
#' Each entry corresponds to a row of \code{combinations}.
#' \item \code{index}, an integer vector of length equal to the number of cells in \code{x}.
#' This specifies the combination in \code{combinations} to which each cell was assigned.
#' }
#'
#' @author Aaron Lun
#' @examples
#' # Mocking a matrix:
#' library(Matrix)
#' y <- round(abs(rsparsematrix(1000, 100, 0.1) * 100))
#'
#' # Simple aggregation:
#' clusters <- sample(LETTERS, 100, replace=TRUE)
#' agg <- aggregateAcrossCells(y, list(cluster=clusters))
#' str(agg)
#'
#' # Multi-factor aggregation
#' samples <- sample(1:5, 100, replace=TRUE)
#' agg2 <- aggregateAcrossCells(y, list(cluster=clusters, sample=samples))
#' str(agg2)
#'
#' @export
aggregateAcrossCells <- function(x, factors, num.threads = 1) {
  f <- lapply(factors, factor)
  f0 <- lapply(f, as.integer)

  x <- beachmat::initializeCpp(x)
  output <- aggregate_across_cells(x, f0, num.threads)
  rownames(output$sums) <- rownames(output$detected) <- rownames(x)

  for (i in seq_along(output$combinations)) {
    current <- output$combinations[[i]]
    output$combinations[[i]] <- levels(f[[i]])[current]
  }

  names(output$combinations) <- names(factors)
  output$combinations <- data.frame(output$combinations)
  output$index <- output$index + 1L

  output
}


#' @importFrom SummarizedExperiment colData<-
#' @importFrom SingleCellExperiment applySCE
.aggregateCells <- function(cellNum, expMatrix, peakMatrix, caller_env, useDim,
                            exp_assay, peak_assay, BPPARAM, clusters=NULL){
  message("performing pseudobulk using an average of ", cellNum, " cells")

  if (!is.null(clusters)) {
    kclusters <- list()
    barcodes <- list()

    #add cluster info to expMatrix
    colData(expMatrix)[, "cluster_for_pseudobulk"] <- clusters

    # K-means clustering
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

  expMatrix <- aggregateAcrossCellsFast(expMatrix, clusters = kclusters, fun_name="sum",
                                         assay.name = exp_assay)
  peakMatrix <- aggregateAcrossCellsFast(peakMatrix, clusters = kclusters, fun_name="sum",
                                          assay.name = peak_assay)

  if (!is.null(clusters))
    clusters <- colData(expMatrix)[, "cluster_for_pseudobulk"]
  assign("expMatrix",expMatrix,envir = caller_env)
  assign("peakMatrix",peakMatrix,envir = caller_env)
  assign("clusters",clusters,envir = caller_env)
}


#' Aggregate cells in SingleCellExperiment
#'
#' Aggregate expression values across cells in SingleCellExperiment based on a
#' grouping factor. This is primarily used to create pseudo-bulk profiles
#' for each cluster/sample combination. It is wrapped around `aggregateAcrossCells`,
#' which relies on the C++ code.
#'
#' @param sce A SingleCellExperiment, SummarizedExperiment or RangedSummarizedExperiment object
#' @param clusters A vector used as a grouping variable. The length should be equal to
#' the number of cells.
#' @param assay.name A character indicating the name of the assay containing the
#' values to be aggregated.
#' @param fun_name A character indicating the function used to aggregate data. The
#' selection is restricted to "mean" or "sum".
#' @param num.threads Integer specifying the number of threads to be used for aggregation.
#' @param aggregateColData A logical specifying if the columns in the `colData`
#' should be included in the output object. Only those columns are selected which
#' can be decomposed by grouping variable into the vectors whose all elements
#' are the same.
#' @return A SingleCellExperiment object containing aggregated cells.
#'
#' @importFrom stats setNames
#' @importFrom SingleCellExperiment altExps altExps<- altExpNames
#' @examples
#' # create a mock singleCellExperiment object for gene expression matrix
#' set.seed(1000)
#' example_sce <- scuttle::mockSCE()
#' ids <- sample(LETTERS[1:5], ncol(example_sce), replace=TRUE)
#' out <- aggregateAcrossCellsFast(example_sce, ids)
#' @export
aggregateAcrossCellsFast <- function(sce, clusters, assay.name="counts", fun_name=c("mean", "sum"),
                                     num.threads=1, aggregateColData = TRUE) {
  .validate_input_sce(SCE=sce, assay_name=assay.name,
                      accepted_classes = c("SingleCellExperiment", "SummarizedExperiment", "RangedSummarizedExperiment"))
  .validate_clusters(clusters, sce)
  clusters <- as.vector(clusters)
  fun_name <- match.arg(fun_name, several.ok = FALSE)
  # aggregate counts in assay
  if(!is.null(assay.name)) {
      x <- setNames(assays(sce)[assay.name], assay.name)
  }
  else {
      x <- setNames(assays(sce), names(assays(sce)))
  }
  aggr.counts <- lapply(x, aggregateAcrossCells, factors = list(clusters), num.threads=num.threads)
  if(fun_name=="sum") {
      assay_matrices <- setNames(lapply(aggr.counts, "[[", "sums"), names(x))
  }
  else {
      assay_matrices <- setNames(lapply(aggr.counts, function(x) t(t(x$sums)/x$counts)), names(x)) #mean
  altExps_list <- NULL
  }
  if(is(sce, "SingleCellExperiment") && length(altExps(sce))>0){
    altExps_list <- lapply(altExps(sce), aggregateAcrossCellsFast, clusters, NULL, fun_name,
                           FALSE)
    names(altExps_list) <- altExpNames(sce)
  }

  # reassemble the singleCellExperiment object
  sce.bulk <- SingleCellExperiment(assay_matrices,
                                   rowData = rowData(sce))
  rownames(colData(sce.bulk)) <- colData(sce.bulk)$idx <- aggr.counts[[1]]$combinations[,1]
  colData(sce.bulk)$ncells <- aggr.counts[[1]]$counts
  if(is(sce, "SingleCellExperiment") && length(altExps(sce))>0){
      altExps(sce.bulk) <- altExps_list
  }
  if(aggregateColData){
    colData.sce.consistent <- .select_consistent_columns(colData(sce), clusters)
    if(!is.null(colData.sce.consistent)){
        duplicated_colnames <- intersect(colnames(colData(sce.bulk)), colnames(colData.sce.consistent))
        if(length(duplicated_colnames)>0){
            stop(spritnf("The following columns are already present in the colData: %s", paste(duplicated_colnames,collapse=", ")))
        }
        # one row per cluster
        unique_clusters_idx <- match(aggr.counts[[1]]$combinations[,1], clusters)
        colData(sce.bulk) <- cbind(colData(sce.bulk), colData.sce.consistent[unique_clusters_idx,,drop=FALSE])
    }
  }
  rowRanges(sce.bulk) <- rowRanges(sce)
  sce.bulk
}

.select_consistent_columns <- function(df, ids){
    if(is(df,"DataFrame")){ # handle DataFrame separately to preserve its hierarchical structure
        current_col = 0
        for(i in seq_len(ncol(df))){
            current_col = current_col+1
            if(length(dim(df[,current_col]))<2) { # select vectors and 1-dim arrays
                if(!.is_consistent(df[,current_col],ids)) {
                    df[[current_col]] <- NULL
                    current_col = current_col-1
                }
            }
            else{
                consistent_data <- .select_consistent_columns(df[,current_col], ids)
                df[[current_col]] <- consistent_data
                if(is.null(consistent_data)){
                    current_col = current_col-1
                }
            }
        }
    }
    else{
        preserve_columns <- integer(0)
        for(i in seq_len(ncol(df))){
            if(.is_consistent(df[,i],ids)){
                preserve_columns <- c(preserve_columns, i)
            }
        }
        df <- df[,preserve_columns,drop=FALSE]
    }
    if(ncol(df)==0){
        df <- NULL
    }
    df
}

.is_consistent <- function(x1, x2){
  all(unlist(lapply(split(x1, x2), function(x) length(unique(x))==1)))
}
