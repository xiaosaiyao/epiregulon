.validate_clusters <- function(clusters, expMatrix){
  clusters <- tryCatch(as.vector(clusters), error=function(cond) {
    message("'clusters' argument should be coercible to a vector")
    stop(cond)
  })
  if (!is.atomic(as.vector(clusters))) {
      stop("'clusters' argument should be coercible to an atomic vector")
  }
  if (length(as.vector(clusters)) != ncol(expMatrix)) {
    stop("'clusters' length should be equal to the number of cells")
  }
  if (any(is.na(as.vector(clusters)))) {
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
                                accepted_classes=c("SingleCellExperiment", "RangedSummarizedExperiment"),
                                unique_features=FALSE){
    checkmate::assert_multi_class(SCE, accepted_classes)
    stopifnot(assay_name %in% names(assays(SCE)))
    data_object_name <- as.character(substitute(SCE))
    if (any(dim(SCE) == 0)){
        stop(sprintf("%s with no data", data_object_name))
    }
    if (unique_features){
        if (any(duplicated(rownames(SCE)))){
            stop(sprintf("Feature names in %s should be unique", data_object_name))
        }
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

.validate_regulon <- function(regulon, required_columns=c("tf", "target", "idxATAC")){
    checkmate::assert_multi_class(regulon, c("DataFrame", "data.frame", "DFrame"))
    missing_cols <- setdiff(required_columns, colnames(regulon))
    if (length(missing_cols)>1) {
        stop("The following column(s) are missing in the regulon object: ", paste(missing_cols, collapse=", "))
    }
    if (nrow(regulon)==0){
        stop("regulon should not be empty")
    }
    columns_with_NA <- required_columns[unlist(lapply(regulon[,required_columns], function(x) any(is.na(x))))]
    if (length(columns_with_NA)>0){
        warning("The following regulon column(s) contain NA value(s): ", paste(columns_with_NA, collapse=", "))
    }
}
