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
    if (prop < 1e-04 | prop > 0.9999)
      warning(sprintf("Strong inbalance between groups after applying cutoff to peakMatrix. Consider %s value of the peak_cutoff",
                      c("increasing", "decreasing")[(prop < 1e-04) + 1]))
  }
  if (!is.null(exp_cutoff)) {
    prop <- sum(expMatrix > exp_cutoff)/prod(dim(expMatrix))
    if (prop < 1e-04 | prop > 0.9999)
      warning(sprintf("Strong inbalance between groups after applying cutoff to expMatrix. Consider %s value of the exp_cutoff",
                      c("increasing", "decreasing")[(prop < 1e-04) + 1]))
  }
}

.validate_input_sce <- function(expMatrix, exp_assay, peakMatrix=NULL, peak_assay=NULL, tf_re.merge=FALSE){
  checkmate::assert_class(expMatrix, "SingleCellExperiment")
  stopifnot(exp_assay %in% names(assays(expMatrix)))
  if (any(dim(expMatrix) == 0)) stop("SingleCellExperiment with no data")
  if(tf_re.merge | !is.null(peakMatrix)){
    checkmate::assert_class(peakMatrix, "SingleCellExperiment")
    stopifnot(peak_assay %in% names(assays(peakMatrix)))
    stopifnot(ncol(peakMatrix) == ncol(expMatrix))
    if(nrow(peakMatrix)==0) stop("peakMatrix with no data")
  }
}
