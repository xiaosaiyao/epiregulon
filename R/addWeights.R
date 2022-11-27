#' Calculate weights for the regulons by computing co-association between TF and target gene expression
#'
#' @param regulon A data frame consisting of tf (regulator) and target in the column names.
#' @param expMatrix A SingleCellExperiment object containing gene expression information
#' @param peakMatrix A SingleCellExperiment object or matrix containing peak accessibility with
#' peaks in the rows and cells in the columns
#' @param exp_assay String specifying the name of the assay to be retrieved from the SingleCellExperiment object
#' @param peak_assay String indicating the name of the assay in peakMatrix for chromatin accessibility
#' @param method String specifying the method of weights calculation. Four options are available: `corr`,`MI`, `lmfit`,`wilcoxon` and `logFC`.
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
#' This function estimates the regulatory potential of transcription factor on its target genes, or in other words,
#' the magnitude of gene expression changes induced by transcription factor activity, using one of the four methods:
#' \itemize{
#' \item{`corr` - correlation between TF and target gene expression}
#' \item{`lmfit` - coefficients of target genes estimated from linear regression of TF ~ TG}
#' \item{`MI` - mutual information between the TF and target gene expression}
#' \item{`wilcoxon` - effect size of the Wilcoxon test between target gene expression in cells jointly expressing all 3 elements vs
#' cells that do not}
#' \item{`logFC` - log 2 fold difference of target gene expression in cells jointly expressing all 3 elements vs cells that do not}
#' }
#' Four measures (`corr`, `lmfit`, `wilcoxon` and `logFC`) give both the magnitude and directionality of changes whereas `MI` always outputs
#' positive weights. The correlation, linear fit and mutual information statistics are computed on the pseudobulked gene expression or accessibility
#' matrices, whereas the Wilcoxon and log fold change group cells based on the joint expression of TF, RE and TG in each single cell.
#'
#' When using the `corr` method, the default practice is to compute weights by correlating the pseudobulk target gene expression vs
#' the pseudobulk TF gene expression. However, often times, an inhibitor of TF does not alter the gene expression of the TF.
#' In rare cases, cells may even compensate by increasing the expression of the TF. In this case, the activity of the TF,
#' if computed by TF-TG correlation, may show a spurious increase in its activity. As an alternative to gene expression,
#' we may correlate the product of TF and RE against TG. When `tf_re.merge` is `TRUE`, we take the product of
#' the gene expression and chromatin accessibility.
#'

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
#'
#' # add weights to regulon
#' regulon.w <- addWeights(regulon = regulon, expMatrix = expMatrix, exp_assay = "logcounts",
#' peakMatrix = peakMatrix, peak_assay = "counts", clusters = expMatrix$cluster,
#' min_targets = 5, method = "logFC")
#

#' @author Xiaosai Yao, Shang-yang Chen, Tomasz Wlodarczyk


addWeights <- function(regulon,
                       expMatrix = NULL,
                       peakMatrix = NULL,
                       exp_assay = "logcounts",
                       peak_assay = "PeakMatrix",
                       method = c("corr", "MI", "lmfit","logFC","wilcoxon"),
                       clusters,
                       exp_cutoff = 1,
                       peak_cutoff = 0,
                       block_factor = NULL,
                       aggregation_function = mean,
                       min_targets = 10,
                       tf_re.merge = FALSE,
                       BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)){


  # choose method
  method <- match.arg(method)
  message("adding weights using ", method, "...")

  # extract matrices from SE
  if (checkmate::test_class(expMatrix,classes = "SummarizedExperiment")){
    expMatrix <- assay(expMatrix, exp_assay)
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
    expMatrix <- expMatrix[!Matrix::rowSums(is.na(expMatrix)) == ncol(expMatrix), ]


    if (tf_re.merge) {
      averages.se.peak <- scuttle::sumCountsAcrossCells(
        peakMatrix,
        ids = groupings,
        average = TRUE,
        BPPARAM = BPPARAM
      )

      # average accessibility across pseudobulk clusters
      peakMatrix <- assays(averages.se.peak)$average

    }


  } else if (method %in% c("logFC", "wilcoxon")){
    message("binarizing matrices...")
    peakMatrix <- binarize_matrix(peakMatrix, peak_cutoff)
    tfMatrix <- binarize_matrix(expMatrix, exp_cutoff)
  }

  # order regulon


  regulon <- regulon[order(regulon$tf),]

  # remove tfs not found in expression matrix
  regulon <- regulon[which(regulon$tf %in% rownames(expMatrix)),]

  # remove targets not found in expression matrix
  regulon <- regulon[which(regulon$target %in% rownames(expMatrix)),]

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

    n_pseudobulk <- length(unique(clusters))

    output_df <- BiocParallel::bplapply(
      X = seq_len(length(regulon.split)),
      FUN = use_MI_method,
      regulon.split,
      expMatrix,
      peakMatrix,
      tf_re.merge,
      n_pseudobulk,
      BPPARAM = BPPARAM
    )


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


    output_df <- #BiocParallel::bp
    lapply(X = seq_len(length(regulon.split)),
                                        FUN = compare_logFC_bp,
                                        regulon.split,
                                        expMatrix,
                                        tfMatrix,
                                        peakMatrix#,
                                        #BPPARAM = BPPARAM
           )

  } else if (method == "wilcoxon") {
    message("calculating rank...")
    tg_rank <- t(apply(expMatrix, 1, rank, ties.method = "average"))
    tie <- apply(tg_rank, 1, function(x){
      freq <- table(x)
      sum((freq^3 - freq)/12)
    })
    message("performing Mann-Whitney U-Test...")
    output_df <- BiocParallel::bplapply(X = seq_len(length(regulon.split)),
                                        FUN = compare_wilcox_bp,
                                        regulon.split,
                                        expMatrix,
                                        tfMatrix,
                                        peakMatrix,
                                        tg_rank,
                                        tie,
                                        BPPARAM = BPPARAM)



  } else {

    stop("method should be corr, MI, lmfit, logFC or wilcoxon")

  }

  output_df <- do.call(rbind, output_df)


  if (method == "wilcoxon") {
    # Calculate effect size for wilcoxon
    n_cells <- ncol(expMatrix)
    # if groups have the same ranks the result will be NaN
    output_df$weight[is.nan(output_df$weight)] <- 0
    # transform z-scores to effect size
    output_df$weight <- output_df$weight/sqrt(n_cells)

  }


  ## Aggregate by REs to have only TF-TG weights
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



use_MI_method <- function(n,
                          regulon.split,
                          expMatrix,
                          peakMatrix,
                          tf_re.merge,
                          n_pseudobulk){

  if (tf_re.merge){
    tf_re <- expMatrix[regulon.split[[n]]$tf, , drop = FALSE] * peakMatrix[regulon.split[[n]]$idxATAC, , drop = FALSE]
  } else {
    tf_re <- expMatrix[regulon.split[[n]]$tf, , drop = FALSE]
  }
  tg <- expMatrix[regulon.split[[n]]$target, , drop = FALSE]

  regulon.split[[n]]$weight <- mapply(MI_per_row, as.data.frame(t(tf_re)), as.data.frame(t(tg)))
  regulon.split[[n]]
}

MI_per_row <- function (tf_re, tg){
  y2d <- entropy::discretize2d(tf_re,
                               tg,
                               numBins1 =  max(10, unique(tf_re)),
                               numBins2 =  max(10, unique(tg)))
  MI <- entropy::mi.empirical(y2d)
}



use_lmfit_method <- function(n,
                            regulon.split,
                            expMatrix,
                            peakMatrix,
                            tf_re.merge){

  if (tf_re.merge){
    tf_re <- expMatrix[regulon.split[[n]]$tf, , drop = FALSE] * peakMatrix[regulon.split[[n]]$idxATAC, , drop = FALSE]
  } else {
    tf_re <- expMatrix[regulon.split[[n]]$tf, , drop = FALSE]
  }
  tg <- expMatrix[regulon.split[[n]]$target, , drop = FALSE]
  regulon.split[[n]]$weight <- mapply(function(y,x) {stats::cov(x,y)/stats::var(x)}, as.data.frame(t(tf_re)), as.data.frame(t(tg)))
  regulon.split[[n]]
}



compare_logFC_bp <- function(n,
                             regulon.split,
                             expMatrix,
                             tfMatrix,
                             peakMatrix){

  expMatrix <- expMatrix[regulon.split[[n]]$target,,drop = FALSE]
  tf_reMatrix <- tfMatrix[regulon.split[[n]]$tf,,drop = FALSE] *
    peakMatrix[regulon.split[[n]]$idxATAC,,drop = FALSE]

  group1 <- tf_reMatrix
  group0 <- (1-tf_reMatrix)

  regulon.split[[n]]$weight <- Matrix::rowSums(expMatrix * group1) / Matrix::rowSums(group1) -
    Matrix::rowSums(expMatrix * group0) / Matrix::rowSums(group0)

  return(regulon.split[[n]])

}


compare_wilcox_bp <- function(n,
                              regulon.split,
                              expMatrix,
                              tfMatrix,
                              peakMatrix,
                              tg_rank,
                              tie){
  expMatrix <- expMatrix[regulon.split[[n]]$target,,drop = FALSE]
  tf_reMatrix <- tfMatrix[regulon.split[[n]]$tf,,drop = FALSE] *
    peakMatrix[regulon.split[[n]]$idxATAC,,drop = FALSE]

  tg_rank <- tg_rank[regulon.split[[n]]$target,,drop = FALSE]
  tie <- tie[regulon.split[[n]]$target]
  regulon.split[[n]]$weight <- wilcoxTest(tf_reMatrix, tg_rank, tie)
  return(regulon.split[[n]])
}

# expMatrix <- scuttle::mockSCE()
# expMatrix <- scuttle::logNormCounts(expMatrix)
# expMatrix <- SingleCellExperiment::logcounts(expMatrix)
# expMatrix <- as(expMatrix, "dgCMatrix")
# tg_rank <- t(apply(expMatrix,1, rank, ties.method = "average"))
#
# regulon <- data.frame(tf = c(rep("Gene_0001",5), rep("Gene_0002",10)),
#                       idxATAC = 1:15,
#                       target = c(paste0("Gene_000",2:6), paste0("Gene_00",11:20)))
# tg_rank <- tg_rank[regulon$target,]
#
# tf_reMatrix <- as(matrix(rbinom(15*200,1,0.1), 15,200),"dgCMatrix")
# rownames(tf_reMatrix) <- rownames(tg_rank)
# colnames(tf_reMatrix) <- rownames(colnames)
#
#
# wilcox_stats <- wilcox(tf_reMatrix, tg_rank)

wilcoxTest <- function(tf_reMatrix, tg_rank, tie) {

  n <- ncol(tf_reMatrix)
  n1 <- Matrix::rowSums(tf_reMatrix)
  n2 <- n - n1

  T1 <- Matrix::rowSums(tg_rank * tf_reMatrix)
  T2 <- Matrix::rowSums(tg_rank * (1-tf_reMatrix))

  U1 <- n1*n2 + n1*(n1+1)/2 - T1
  U2 <- n1*n2 + n2*(n2+1)/2 - T2


  U <- cbind(U1, U2)
  U <- apply(U, 1, min)

  mu <- n1*n2/2


  sigma <- sqrt(n1*n2/n/(n-1)) * sqrt((n^3-n)/12 - tie)
  stats <- (U-mu)/sigma * sign(U1-U2)

}



