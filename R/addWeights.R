#' Calculate weights for the regulons by computing co-association between TF and target gene expression
#'
#' @param regulon A DataFrame object consisting of tf (regulator) and target in the column names.
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
#' @return A DataFrame with columns of corr and/or MI added to the regulon. TFs not found in the expression matrix and regulons not
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
#' expMatrix$cluster <- sample(LETTERS[1:5], ncol(expMatrix), replace=TRUE)
#'
#' # create a mock singleCellExperiment object for peak matrix
#' peakMatrix <- scuttle::mockSCE()
#' rownames(peakMatrix) <- 1:2000
#'
#' # create a mock regulon
#' regulon <- DataFrame(tf=c(rep("Gene_0001",5), rep("Gene_0002",10)),
#'                       idxATAC=1:15,
#'                       target=c(paste0("Gene_000",2:6), paste0("Gene_00",11:20)))
#'
#' # add weights to regulon
#' regulon.w <- addWeights(regulon=regulon, expMatrix=expMatrix, exp_assay="logcounts",
#' peakMatrix=peakMatrix, peak_assay="counts", clusters=expMatrix$cluster,
#' min_targets=5, method="logFC")
#

#' @author Xiaosai Yao, Shang-yang Chen, Tomasz Wlodarczyk


addWeights <- function(regulon,
                       expMatrix=NULL,
                       peakMatrix=NULL,
                       exp_assay="logcounts",
                       peak_assay="PeakMatrix",
                       method=c("corr", "MI", "lmfit","logFC","wilcoxon"),
                       clusters=NULL,
                       exp_cutoff=1,
                       peak_cutoff=0,
                       block_factor=NULL,
                       aggregation_function=mean,
                       min_targets=10,
                       tf_re.merge=FALSE,
                       BPPARAM=BiocParallel::SerialParam(progressbar=TRUE)){


  # choose method
  method <- match.arg(method)
  message("adding weights using ", method, "...")

  # extract matrices from SE
  if (checkmate::test_class(expMatrix,classes="SummarizedExperiment")){
    expMatrix <- assay(expMatrix, exp_assay)
  }

  if(any(dim(expMatrix)==0)) stop("expMatrix with no data")

  if (checkmate::test_class(peakMatrix, classes="SummarizedExperiment")){
    peakMatrix <- assay(peakMatrix, peak_assay)
  }

  if(!is.null(peakMatrix)){
    checkmate::testMultiClass(peakMatrix, c("matrix", "dgeMatrix",
                                            "lgCMatrix", "dgCMatrix"))
  }

  checkmate::testMultiClass(regulon, c("data.frame", "DFrame"))

  if(nrow(regulon) == 0) stop("Regulon with zero rows")

  checkmate::assert_logical(tf_re.merge, len = 1)

  if(!is.null(clusters)){
    if(!is.character(clusters) | !is.vector(clusters))
    tryCatch(clusters <- as.character(as.vector(clusters)), error = function(e) stop("'clusters' agrument should be coercible to a character vector"))
  }

  if(method %in% c("logFC", "wilcoxon") | tf_re.merge){
    if(is.null(peakMatrix)) stop("Peak matrix should be provided")
    if(any(dim(peakMatrix) == 0)) stop("Peak matrix is empty")
  }


  expMatrix <- as(expMatrix, "dgCMatrix")

  regulon <- S4Vectors::DataFrame(regulon)

  # order regulon
  regulon <- regulon[order(regulon$tf),]

  # remove genes not found in regulon
  expMatrix <- expMatrix[which(rownames(expMatrix) %in% unique(c(regulon$tf, regulon$target))),]

  keep <- regulon$tf %in% rownames(expMatrix) &
    regulon$target %in% rownames(expMatrix)

  regulon <- regulon[keep,]

  if(nrow(regulon) == 0) stop("Gene names in the regulon should match those in the expMatrix")

  # remove tfs with less than min_targets
  regulon <- regulon[regulon$tf %in% names(which(table(regulon$tf) >= min_targets)),]

  # if a cluster is named "all", replace it to distinguish from all cells
  if (!is.null(clusters)){
    clusters[clusters == "all"] <- "clusters_all"

  }
  unique_clusters <- sort(unique(clusters))

  # define weight matrix
  if (method %in%  c("corr", "MI", "lmfit")) {
    regulon$weight <- NA
    regulon.split <- split(regulon, regulon$tf)

  }  else if (method %in% c("logFC", "wilcoxon")){
    regulon$weight <- initiateMatCluster(clusters, nrow = nrow(regulon))
    regulon.split <- split(regulon, regulon$tf)
  }

  if (!is.null(peakMatrix)) {
    peakMatrix <- as(peakMatrix, "dgCMatrix")
    # name peakMatrix
    rownames(peakMatrix) <- seq_len(nrow(peakMatrix))
    # remove peaks not found in regulon
    peakMatrix <- peakMatrix[which(rownames(peakMatrix) %in% unique(regulon$idxATAC)),]
  }

  if (method == "wilcoxon") {
    keep <- regulon$idxATAC >= 1 & regulon$idxATAC <= nrow(peakMatrix)
    regulon <- regulon[keep,]
    peakMatrix <- binarize_matrix(peakMatrix, cutoff = peak_cutoff)
    copy <- regulon
    all.targets <- sort(unique(regulon$target))
    all.tfs <- sort(unique(regulon$tf))
    copy$tf <- match(copy$tf, all.tfs)
    copy$target <- match(copy$target, all.targets)
    if (!is.null(clusters)){
      # binarize expression matrix for each cluster separately
      expMatrix_tfs_clusters <- expMatrix[all.tfs,,drop = FALSE]
      for(cluster in unique(clusters)){
        cluster_ind <- which(clusters == cluster)
        expMatrix_tfs_clusters[,cluster_ind, drop = FALSE] <- binarize_matrix(expMatrix_tfs_clusters[,cluster_ind, drop = FALSE],
                                                                 cutoff = exp_cutoff)
      }
    }
    expMatrix_tfs <- binarize_matrix(expMatrix[all.tfs,,drop = FALSE], cutoff = exp_cutoff)
    exprs_trans_target <- Matrix::t(expMatrix[all.targets,,drop=FALSE])
    exprs_trans_tf <- Matrix::t(expMatrix_tfs)
    all.peaks <- sort(unique(copy$idxATAC))
    peak_trans <- Matrix::t(peakMatrix[all.peaks,,drop=FALSE])
    if (!is(peak_trans, "dgCMatrix")) {
        peak_trans <- as(peak_trans, "dgCMatrix")
    }
    copy$idxATAC <- match(copy$idxATAC, all.peaks)
    reg.order <- order(copy$target, copy$tf, copy$idxATAC)
    copy <- copy[reg.order,,drop=FALSE]
    # calculate stats for all clusters
    output <- fast_wilcox(
      exprs_x = exprs_trans_target@x,
      exprs_i = exprs_trans_target@i,
      exprs_p = exprs_trans_target@p,
      exprs_tf_x = as.logical(exprs_trans_tf@x),
      exprs_tf_i = exprs_trans_tf@i,
      exprs_tf_p = exprs_trans_tf@p,
      peak_x = peak_trans@x,
      peak_i = peak_trans@i,
      peak_p = peak_trans@p,
      target_id = copy$target - 1L,
      tf_id = copy$tf - 1L,
      peak_id = copy$idxATAC - 1L,
      clusters = integer(0),
      cell_numb = nrow(exprs_trans_target)
    )

    if(!is.null(clusters)){
      # calculate stats for each cluster separately
      exprs_trans_tf_clusters <- Matrix::t(expMatrix_tfs_clusters)
      fclusters <- factor(clusters)
      fclusters_order <- order(levels(fclusters))
      iclusters <- as.integer(fclusters)
      output_clusters <- fast_wilcox(
        exprs_x = exprs_trans_target@x,
        exprs_i = exprs_trans_target@i,
        exprs_p = exprs_trans_target@p,
        exprs_tf_x = as.logical(exprs_trans_tf_clusters@x),
        exprs_tf_i = exprs_trans_tf_clusters@i,
        exprs_tf_p = exprs_trans_tf_clusters@p,
        peak_x = peak_trans@x,
        peak_i = peak_trans@i,
        peak_p = peak_trans@p,
        target_id = copy$target - 1L,
        tf_id = copy$tf - 1L,
        peak_id = copy$idxATAC - 1L,
        clusters = iclusters - 1L,
        cell_numb = nrow(exprs_trans_target)
      )
      output <- mapply(function(x,y) rbind(x,y), output, output_clusters, SIMPLIFY = FALSE)
    }

    AUC <- output$auc
    ties <- output$ties
    n1 <- output$total0
    n2 <- output$total1

    prod <- n1 * n2
    n <- n1 + n2 # technically same across all rows, but we'll just do this for simplicity.
    sigma <- sqrt(prod / 12 * (n + 1 - ties / n / (n-1)))

    mu <- prod/2
    stats <- (AUC - mu)/sigma
    # set z-score to zero if the size of the of the groups is equal to 0
    stats[n1==0 | n2==0] <- 0
    stats[,reg.order] <- stats
    regulon$weight <- t(stats)
    # Calculate effect size
    n_cells <- ncol(expMatrix)
    n_cells <- c(n_cells, table(clusters))
    # transform z-scores to effect size
    regulon$weight <- t(t(regulon$weight)/sqrt(n_cells))
    return(regulon)
  }

  if (method %in% c("corr", "MI", "lmfit")){
    message("calculating average expression across clusters...")

    # define groupings
    groupings <- S4Vectors::DataFrame(cluster=clusters)
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


  }



  message("computing weights...")
  if (method == "corr") {

    output_df <- BiocParallel::bplapply(
      X=seq_len(length(regulon.split)),
      FUN=use_corr_method,
      regulon.split,
      expMatrix,
      peakMatrix,
      tf_re.merge,
      BPPARAM=BPPARAM)

  } else if (method == "MI") {

    n_pseudobulk <- length(unique(clusters))

    output_df <- BiocParallel::bplapply(
      X=seq_len(length(regulon.split)),
      FUN=use_MI_method,
      regulon.split,
      expMatrix,
      peakMatrix,
      tf_re.merge,
      n_pseudobulk,
      BPPARAM=BPPARAM
    )


  } else if (method == "lmfit") {

    output_df <- BiocParallel::bplapply(
      X=seq_len(length(regulon.split)),
      FUN=use_lmfit_method,
      regulon.split,
      expMatrix,
      peakMatrix,
      tf_re.merge,
      BPPARAM = BPPARAM)

  } else if (method == "logFC") {


    for (cluster in c("all", unique_clusters)) {
      if (cluster == "all" & !is.null(clusters)){
        cluster.current <- unique_clusters
      } else if (cluster == "all" & is.null(clusters)) {
        cluster.current <- "all"
        clusters <- rep("all", ncol(expMatrix))
      } else if (cluster != "all"){
        cluster.current <- cluster
        regulon.split <- output_df
      }

      peakMatrix.bi <- binarize_matrix(peakMatrix[,clusters %in% cluster.current], peak_cutoff)
      tfMatrix.bi <- binarize_matrix(expMatrix[,clusters %in% cluster.current], exp_cutoff)
      output_df <- BiocParallel::bplapply(X=seq_len(length(regulon.split)),
                            FUN=compare_logFC_bp,
                            regulon.split,
                            expMatrix[,clusters %in% cluster.current],
                            tfMatrix.bi,
                            peakMatrix.bi,
                            cluster,
                            BPPARAM=BPPARAM
      )
    }


  } else {

    stop("method should be corr, MI, lmfit, logFC or wilcoxon")

  }

  output_df <- do.call(rbind, output_df)



  ## Aggregate by REs to have only TF-TG weights
  #regulon <- stats::aggregate(weight~tf+target, FUN = aggregation_function, na.rm = TRUE, data = output_df)
  #regulon[order(regulon$tf),]

  #regulon <- output_df

  return(output_df)

}

#' @keywords internal

use_corr_method <- function(n,
                            regulon.split,
                            expMatrix,
                            peakMatrix,
                            tf_re.merge,
                            BPPARAM = BPPARAM){

  if (tf_re.merge){
    tf_re <- expMatrix[regulon.split[[n]]$tf, , drop=FALSE] * peakMatrix[as.character(regulon.split[[n]]$idxATAC), , drop=FALSE]
  } else {
    tf_re <- expMatrix[regulon.split[[n]]$tf, , drop=FALSE]
  }
  tg <- expMatrix[regulon.split[[n]]$target, , drop=FALSE]
  regulon.split[[n]]$weight <- mapply(stats::cor, as.data.frame(t(tf_re)), as.data.frame(t(tg)), use="everything")
  regulon.split[[n]]
}



use_MI_method <- function(n,
                          regulon.split,
                          expMatrix,
                          peakMatrix,
                          tf_re.merge,
                          n_pseudobulk){

  if (tf_re.merge){
    tf_re <- expMatrix[regulon.split[[n]]$tf, , drop=FALSE] * peakMatrix[as.character(regulon.split[[n]]$idxATAC), , drop=FALSE]
  } else {
    tf_re <- expMatrix[regulon.split[[n]]$tf, , drop=FALSE]
  }
  tg <- expMatrix[regulon.split[[n]]$target, , drop=FALSE]

  regulon.split[[n]]$weight <- mapply(MI_per_row, as.data.frame(t(tf_re)), as.data.frame(t(tg)))
  regulon.split[[n]]
}

MI_per_row <- function (tf_re, tg){
  y2d <- entropy::discretize2d(tf_re,
                               tg,
                               numBins1=max(10, unique(tf_re)),
                               numBins2=max(10, unique(tg)))
  MI <- entropy::mi.empirical(y2d)
}



use_lmfit_method <- function(n,
                             regulon.split,
                             expMatrix,
                             peakMatrix,
                             tf_re.merge){

  if (tf_re.merge){
    tf_re <- expMatrix[regulon.split[[n]]$tf, , drop=FALSE] * peakMatrix[as.character(regulon.split[[n]]$idxATAC), , drop = FALSE]
  } else {
    tf_re <- expMatrix[regulon.split[[n]]$tf, , drop=FALSE]
  }
  tg <- expMatrix[regulon.split[[n]]$target, , drop=FALSE]
  regulon.split[[n]]$weight <- mapply(function(y,x) {stats::cov(x,y)/stats::var(x)}, as.data.frame(t(tf_re)), as.data.frame(t(tg)))
  regulon.split[[n]]
}





compare_logFC_bp <- function(n,
                             regulon.split,
                             expMatrix,
                             tfMatrix,
                             peakMatrix,
                             cluster){

  expMatrix <- expMatrix[regulon.split[[n]]$target,,drop=FALSE]
  tf_reMatrix <- tfMatrix[regulon.split[[n]]$tf,,drop=FALSE] *
    peakMatrix[as.character(regulon.split[[n]]$idxATAC),,drop=FALSE]

  group1 <- tf_reMatrix
  group0 <- (1-tf_reMatrix)

  regulon.split[[n]]$weight[,cluster] <- Matrix::rowSums(expMatrix * group1) / Matrix::rowSums(group1) -
    Matrix::rowSums(expMatrix * group0) / Matrix::rowSums(group0)

  return(regulon.split[[n]])

}
