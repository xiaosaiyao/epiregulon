.validate_clusters <- function(clusters, expMatrix){
  clusters <- tryCatch(as.vector(clusters), error = function(cond) {
    message("'clusters' argument should be coercible to a vector")
    stop(cond)
  })
  if (length(clusters) != ncol(expMatrix)) {
    stop("'clusters' length should be equal to the number of cells")
  }
  if (any(is.na(clusters))) {
    stop("'clusters' object contains NA")
  }
}

.balance_check <- function(peak_cutoff, exp_cutoff, peakMatrix, expMatrix){
  if (!is.null(peak_cutoff)) {
    prop <- sum(peakMatrix > peak_cutoff)/prod(dim(peakMatrix))
    if (prop < 1e-04 | prop > 0.9999) {
      warning(sprintf("Strong inbalance between groups after applying cutoff to peakMatrix. Consider %s value of the peak_cutoff",
                      c("increasing", "decreasing")[(prop < 1e-04) + 1]))
    }
  }
  if (!is.null(exp_cutoff)) {
    prop <- sum(expMatrix > exp_cutoff)/prod(dim(expMatrix))
    if (prop < 1e-04 | prop > 0.9999) {
      warning(sprintf("Strong inbalance between groups after applying cutoff to expMatrix. Consider %s value of the exp_cutoff",
                      c("increasing", "decreasing")[(prop < 1e-04) + 1]))
    }
  }
}


.validate_input_sce <- function(SCE, 
                                assay_name, 
                                row.ranges=FALSE,
                                accepted_classes = c("SingleCellExperiment", "RangedSummarizedExperiment")){
  checkmate::assert_multi_class(SCE, accepted_classes)
  stopifnot(assay_name %in% names(assays(SCE)))
  data_object_name <- as.character(substitute(SCE))
  if (any(dim(SCE) == 0)){
    stop(sprintf("%s with no data", data_object_name))
  }
  if(row.ranges){
    if (length(rowRanges(SCE)) == 0) {
      stop(sprintf("%s should contain non-empty rowRanges", data_object_name))
    }
    checkmate::assert_class(rowRanges(SCE), "GRanges")
    if (length(rowRanges(SCE)) == 0) {
      stop(sprintf("%s should contain non-empty rowRanges", data_object_name))
    }
  }
}

