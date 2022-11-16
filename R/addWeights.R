#' Calculate weights for the regulons by computing co-association between TF and target gene expression
#'
#' @param regulon A data frame consisting of tf (regulator) and target in the column names.
#' @param expMatrix A SingleCellExperiment object containing gene expression information
#' @param peakMatrix A SingleCellExperiment object or matrix containing peak accessibility with
#' peaks in the rows and cells in the columns
#' @param exp_assay String specifying the name of the assay to be retrieved from the SingleCellExperiment object
#' @param peak_assay String indicating the name of the assay in peakMatrix for chromatin accessibility
#' @param method String specifying the method of weights calculation. Four options are available: `corr`,`MI`, `wilcoxon` and `logFC`.
#' @param clusters A vector corresponding to the cluster labels of the cells
#' @param exp_cutoff A scalar indicating the minimum gene expression for transcription factor above which
#' cell is considered as having expressed transcription factor.
#' @param peak_cutoff A scalar indicating the minimum peak accessibility above which peak is
#' considered open.
#' @param block_factor String specifying the field in the colData of the SingleCellExperiment object to be used as blocking factor (such as batch)
#' @param aggregation_function Function being used for summarizing weights from the transcription factor-target gene pair with
#' many regulatory elements.
#' @param min_targets Integer specifying the minimum number of targets for each tf in the regulon with 10 targets as the default
#' @param tf_re.merge A logical to indicate whether to consider both TF expression and chromatin accessibility. See details.
#' @param BPPARAM A BiocParallelParam object specifying whether summation should be parallelized. Use BiocParallel::SerialParam() for
#' serial evaluation and use BiocParallel::MulticoreParam() for parallel evaluation
#'
#' @return A data frame with columns of corr and/or MI added to the regulon. TFs not found in the expression matrix and regulons not
#' meeting the minimal number of targets were filtered out.
#' @import SummarizedExperiment
#' @details

#'
#' @export
#'
#' @examples
#' # create a mock singleCellExperiment object for gene expression matrix
#' expMatrix <- scuttle::mockSCE()
#' expMatrix <- scuttle::logNormCounts(expMatrix)
#' expMatrix$cluster <- sample(LETTERS[1:5], ncol(expMatrix), replace = TRUE)
#'
#' # create a mock singleCellExperiment object for peak matrix
#' peakMatrix <- scuttle::mockSCE()
#' rownames(peakMatrix) <- 1:2000
#'
#' # create a mock regulon
#' regulon <- data.frame(tf = c(rep("Gene_0001",5), rep("Gene_0002",10)),
#'                       idxATAC = 1:15,
#'                       target = c(paste0("Gene_000",2:6), paste0("Gene_00",11:20)))
#' @author Xiaosai Yao, Shang-yang Chen, Tomasz Wlodarczyk


addWeights <- function(regulon,
                       expMatrix = NULL,
                       peakMatrix = NULL,
                       exp_assay = "logcounts",
                       peak_assay = "PeakMatrix",
                       exp_cutoff = 1,
                       peak_cutoff = 0,
                       method = c("corr", "MI", "lmfit","logFC"),
                       clusters,
                       block_factor = NULL,
                       aggregation_function = mean,
                       min_targets = 10,
                       tf_re.merge = FALSE,
                       BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)){



  # choose method
  method <- match.arg(method)
  message("adding weights using ", method)

  # extract matrices from SE
  if (checkmate::test_class(expMatrix,classes = "SummarizedExperiment")){
    expMatrix <- assay(expMatrix, exp_assay)
    expMatrix <- as(expMatrix, "dgCMatrix")
  }

  if (checkmate::test_class(peakMatrix, classes = "SummarizedExperiment")){
    peakMatrix <- assay(peakMatrix, peak_assay)
  }



  if (method %in% c("corr", "MI", "lmfit")){
    message("calculating average expression across clusters...")

    # define groupings
    groupings <- S4Vectors::DataFrame(cluster = clusters)
    if (!is.null(block_factor)) {
      groupings$block <- colData(expMatrix)[block_factor]
    }

    # compute average expression across clusters and batches
    averages.se.exp <- scuttle::sumCountsAcrossCells(
      expMatrix,
      ids = groupings,
      average = TRUE,
      BPPARAM = BPPARAM
    )

    # average expression across pseudobulk clusters
    expMatrix <- assays(averages.se.exp)$average

    # remove genes whose expressions are NA for all pseudobulks
    expMatrix <- expMatrix[!rowSums(is.na(expMatrix)) == ncol(expMatrix), ]


    averages.se.peak <- scuttle::sumCountsAcrossCells(
      peakMatrix,
      ids = groupings,
      average = TRUE,
      BPPARAM = BPPARAM
    )

    # average accessibility across pseudobulk clusters
    peakMatrix <- assays(averages.se.peak)$average

    # remove REs whose accessibility are NA for all pseudobulks
    #peakMatrix <- peakMatrix[!rowSums(is.na(peakMatrix)) == ncol(peakMatrix), ]

  } else if (method %in% c("logFC")){
    message("binarizing matrices...")
    peakMatrix <- binarize_matrix(peakMatrix, peak_cutoff)
    tfMatrix <- binarize_matrix(expMatrix, exp_cutoff)
  }

  # order regulon
  regulon <- regulon[order(regulon$tf, regulon$idxATAC, regulon$target),]

  # remove tfs not found in expression matrix
  regulon <- regulon[which(regulon$tf %in% rownames(expMatrix)),]

  # remove targets not found in expression matrix
  regulon <- regulon[which(regulon$target %in% rownames(expMatrix)),]

  # remove REs not found in peak matrix
  #regulon <- regulon[which(regulon$idxATAC %in% rownames(peakMatrix)),]

  # remove tfs with less than min_targets
  regulon <- regulon[regulon$tf %in% names(which(table(regulon$tf) >= min_targets)),]



  regulon.split <- split(regulon, regulon$tf)


  message("computing weights...")
  if (method == "corr") {

    output_df <- BiocParallel::bplapply(
      X = seq_len(length(regulon.split)),
      FUN = use_corr_method,
      regulon.split,
      expMatrix,
      peakMatrix,
      tf_re.merge,
      BPPARAM = BPPARAM)

  } else if (method == "MI") {

    output_df <- BiocParallel::bplapply(
      X = seq_len(length(regulon.split)),
      FUN = use_MI_method,
      regulon.split,
      expMatrix,
      peakMatrix,
      tf_re.merge,
      BPPARAM = BPPARAM)

  } else if (method == "lmfit") {

    output_df <- BiocParallel::bplapply(
      X = seq_len(length(regulon.split)),
      FUN = use_lmfit_method,
      regulon.split,
      expMatrix,
      peakMatrix,
      tf_re.merge,
      BPPARAM = BPPARAM)

  } else if (method == "logFC") {


    output_df <- BiocParallel::bplapply(X = seq_len(length(regulon.split)),
                                        FUN = compare_logFC_bp,
                                        regulon.split,
                                        expMatrix,
                                        tfMatrix,
                                        peakMatrix)

  }
  else {

    stop("method should be corr, MI or lmfit")

  }

  output_df <- do.call(rbind, output_df)

  ## Aggregate by REs
  regulon <- stats::aggregate(weight~tf+target, FUN = aggregation_function, na.rm = TRUE, data = output_df)
  regulon[order(regulon$tf),]

  regulon <- output_df

  return(regulon)
}

#' @keywords internal

use_corr_method <- function(n,
                            regulon.split,
                            expMatrix,
                            peakMatrix,
                            tf_re.merge,
                            BPPARAM = BPPARAM){

  if (tf_re.merge){
    tf_re <- expMatrix[regulon.split[[n]]$tf, , drop = FALSE] * peakMatrix[regulon.split[[n]]$idxATAC, , drop = FALSE]
  } else {
    tf_re <- expMatrix[regulon.split[[n]]$tf, , drop = FALSE]
  }
  tg <- expMatrix[regulon.split[[n]]$target, , drop = FALSE]
  regulon.split[[n]]$weight <- mapply(stats::cor, as.data.frame(t(tf_re)), as.data.frame(t(tg)), use = "everything")
  regulon.split[[n]]
}


use_lmfit_method <- function(n,
                            regulon.split,
                            expMatrix,
                            peakMatrix,
                            tf_re.merge,
                            BPPARAM = BPPARAM){

  if (tf_re.merge){
    tf_re <- expMatrix[regulon.split[[n]]$tf, , drop = FALSE] * peakMatrix[regulon.split[[n]]$idxATAC, , drop = FALSE]
  } else {
    tf_re <- expMatrix[regulon.split[[n]]$tf, , drop = FALSE]
  }
  tg <- expMatrix[regulon.split[[n]]$target, , drop = FALSE]
  regulon.split[[n]]$weight <- mapply(linearfit, as.data.frame(t(tf_re)), as.data.frame(t(tg)))
  regulon.split[[n]]
}

linearfit <- function(tf,tg){
  fit_res <- lm(tf~tg, data = data.frame(tf = tf, tg = tg))
  fit_res <- fit_res$coefficients[2]
}


compare_logFC_bp <- function(n,
                             regulon.split,
                             expMatrix,
                             tfMatrix,
                             peakMatrix){

  expMatrix <- expMatrix[regulon.split[[n]]$target,,drop = FALSE]
  tf_reMatrix <- tfMatrix[regulon.split[[n]]$tf,,drop = FALSE] *
    peakMatrix[regulon.split[[n]]$idxATAC,,drop = FALSE]

  group1 <- tf_reMatrix * expMatrix
  group0 <- (1-tf_reMatrix) * expMatrix

  regulon.split[[n]]$weight <- Matrix::rowSums(group1)/Matrix::rowSums(tf_reMatrix) -
    Matrix::rowSums(group0)/Matrix::rowSums((1-tf_reMatrix))

  return(regulon.split[[n]])

}


