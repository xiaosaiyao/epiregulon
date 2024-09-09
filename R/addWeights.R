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
#' @param min_targets Integer specifying the minimum number of targets for each tf in the regulon with 10 targets as the default
#' @param tf_re.merge A logical to indicate whether to consider both TF expression and chromatin accessibility. See details.
#' @param aggregateCells A logical to indicate whether to aggregate cells into groups determined by cellNum. This option can be used to
#' overcome data sparsity when using `wilcoxon`.
#' @param useDim String indicating the name of the dimensionality reduction matrix in expMatrix used for cell aggregation
#' @param cellNum An integer specifying the number of cells per cluster for cell aggregation. Default is 10.
#' @param BPPARAM A BiocParallelParam object specifying whether summation should be parallelized. Use BiocParallel::SerialParam() for
#' serial evaluation and use BiocParallel::MulticoreParam() for parallel evaluation
#'
#' @return A DataFrame with columns of corr and/or MI added to the regulon. TFs not found in the expression matrix and regulons not
#' meeting the minimal number of targets were filtered out.
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assay assays colData
#' @importFrom S4Vectors split
#' @details
#' This function estimates the regulatory potential of transcription factor on its target genes, or in other words,
#' the magnitude of gene expression changes induced by transcription factor activity, using one of the four methods:
#' \itemize{
#' \item{`corr` - correlation between TF and target gene expression}
#' \item{`MI` - mutual information between the TF and target gene expression}
#' \item{`wilcoxon` - effect size of the Wilcoxon test between target gene expression in cells jointly expressing all 3 elements vs
#' cells that do not}}
#' Two measures (`corr` and `wilcoxon`) give both the magnitude and directionality of changes whereas `MI` always outputs
#' positive weights. The correlation and mutual information statistics are computed on the pseudobulked gene expression or accessibility
#' matrices, whereas the Wilcoxon method groups cells based on the joint expression of TF, RE and TG in each single cell.
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
#' regulon <- S4Vectors::DataFrame(tf=c(rep('Gene_0001',5), rep('Gene_0002',10)),
#'                       idxATAC=1:15,
#'                       target=c(paste0('Gene_000',2:6), paste0('Gene_00',11:20)))
#'
#' # add weights to regulon
#' regulon.w <- addWeights(regulon=regulon, expMatrix=expMatrix, exp_assay='logcounts',
#' peakMatrix=peakMatrix, peak_assay='counts', clusters=expMatrix$cluster,
#' min_targets=5, method='wilcox')
#'
#' # add weights with cell aggregation
#' expMatrix <- scater::runPCA(expMatrix)
#' regulon.w <- addWeights(regulon=regulon, expMatrix=expMatrix, exp_assay='logcounts',
#' peakMatrix=peakMatrix, peak_assay='counts', clusters=expMatrix$cluster,
#' min_targets=5, method='wilcox', aggregateCells=TRUE, cellNum=3, useDim = 'PCA')
#'
#' @author Xiaosai Yao, Shang-yang Chen, Tomasz Wlodarczyk


addWeights <- function(regulon,
                       expMatrix = NULL,
                       peakMatrix = NULL,
                       exp_assay = "logcounts",
                       peak_assay = "PeakMatrix",
                       method = c("wilcoxon", "corr", "MI"),
                       clusters = NULL,
                       exp_cutoff = 1,
                       peak_cutoff = 0,
                       block_factor = NULL,
                       min_targets = 10,
                       tf_re.merge = FALSE,
                       aggregateCells = FALSE,
                       useDim = "IterativeLSI_ATAC",
                       cellNum = 10,
                       BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)) {
    # choose method
    method <- match.arg(method)
    message("adding weights using ", method, "...")
    checkmate::assert_logical(tf_re.merge, len = 1)
    .validate_input_sce(expMatrix, exp_assay, peakMatrix, peak_assay, tf_re.merge)

    if(!is.null(clusters)) .validate_clusters(clusters, expMatrix)

    checkmate::testMultiClass(regulon, c("data.frame", "DFrame"))

    if (nrow(regulon) == 0)
      stop("Regulon with zero rows")

    # pseudobulk
    if (aggregateCells && method != "wilcoxon")
        message("Cell aggregation is possible only with 'wilcoxon' method.")
    else if (aggregateCells && method == "wilcoxon") .aggregateCells(cellNum, expMatrix, peakMatrix, environment(),
                                                                     useDim, exp_assay, peak_assay, BPPARAM, clusters)

    # extract matrices from SCE
    expMatrix <- assay(expMatrix, exp_assay)

    if(!is.null(peakMatrix)){
      peakMatrix <- assay(peakMatrix, peak_assay)
      checkmate::testMultiClass(peakMatrix, c("matrix", "dgeMatrix", "lgCMatrix",
                                              "dgCMatrix", "CsparseMatrix"))
    }

    if (method=="wilcoxon" & is.null(peakMatrix)) {
        stop("Peak matrix should be provided")
    }

    expMatrix <- as(expMatrix, "CsparseMatrix")

    regulon <- S4Vectors::DataFrame(regulon)

    # order regulon
    regulon <- regulon[order(regulon$tf), , drop=FALSE]

    # remove genes not found in regulon
    expMatrix <- expMatrix[which(rownames(expMatrix) %in% unique(c(regulon$tf, regulon$target))),
        ]

    keep <- regulon$tf %in% rownames(expMatrix) & regulon$target %in% rownames(expMatrix)

    regulon <- regulon[keep, , drop=FALSE]

    if (nrow(regulon) == 0)
        stop("Gene names in the regulon should match those in the expMatrix")

    # remove tfs with less than min_targets
    regulon <- regulon[regulon$tf %in% names(which(table(regulon$tf) >= min_targets)),
        ]

    if (nrow(regulon) == 0) {
        warning("No transcription factor has the required number of target genes")
        return(NULL)
    }

    # if a cluster is named 'all', replace it to distinguish from all cells
    clusters <- renameCluster(clusters)

    unique_clusters <- sort(unique(clusters))

    # define weight matrix
    if (method %in% c("corr", "MI")) {
        regulon$weight <- NA
        regulon.split <- split(regulon, regulon$tf)

    } else if (method %in% c("logFC", "wilcoxon")) {
        regulon$weight <- initiateMatCluster(clusters, nrow = nrow(regulon))
        regulon.split <- split(regulon, regulon$tf)
    }

    if (!is.null(peakMatrix)) {
        peakMatrix <- as(peakMatrix, "CsparseMatrix")
        # name peakMatrix
        rownames(peakMatrix) <- seq_len(nrow(peakMatrix))
    }

    if (method == "wilcoxon") {
        .balance_check(peak_cutoff, exp_cutoff, peakMatrix, expMatrix)
        keep <- regulon$idxATAC >= 1 & regulon$idxATAC <= nrow(peakMatrix)
        regulon <- regulon[keep, , drop=FALSE]
        peakMatrix <- binarize_matrix(peakMatrix, cutoff = peak_cutoff)
        copy <- regulon
        all.targets <- sort(unique(regulon$target))
        all.tfs <- sort(unique(regulon$tf))
        copy$tf <- match(copy$tf, all.tfs)
        copy$target <- match(copy$target, all.targets)
        if (!is.null(clusters)) {
            # binarize expression matrix for each cluster separately
            expMatrix_tfs_clusters <- expMatrix[all.tfs, , drop = FALSE]
            for (cluster in unique(clusters)) {
                cluster_ind <- which(clusters == cluster)
                expMatrix_tfs_clusters[, cluster_ind, drop = FALSE] <- binarize_matrix(expMatrix_tfs_clusters[,
                  cluster_ind, drop = FALSE], cutoff = exp_cutoff)
            }
        }
        expMatrix_tfs <- binarize_matrix(expMatrix[all.tfs, , drop = FALSE], cutoff = exp_cutoff)
        exprs_trans_target <- Matrix::t(expMatrix[all.targets, , drop = FALSE])
        exprs_trans_tf <- Matrix::t(expMatrix_tfs)
        all.peaks <- sort(unique(copy$idxATAC))
        peak_trans <- Matrix::t(peakMatrix[all.peaks, , drop = FALSE])
        if (!is(peak_trans, "CsparseMatrix")) {
            peak_trans <- as(peak_trans, "CsparseMatrix")
        }
        copy$idxATAC <- match(copy$idxATAC, all.peaks)
        reg.order <- order(copy$target, copy$tf, copy$idxATAC)
        copy <- copy[reg.order, , drop = FALSE]
        # calculate stats for all clusters
        output <- fast_wilcox(exprs_x = exprs_trans_target@x, exprs_i = exprs_trans_target@i,
            exprs_p = exprs_trans_target@p, exprs_tf_x = as.logical(exprs_trans_tf@x),
            exprs_tf_i = exprs_trans_tf@i, exprs_tf_p = exprs_trans_tf@p, peak_x = peak_trans@x,
            peak_i = peak_trans@i, peak_p = peak_trans@p, target_id = copy$target -
                1L, tf_id = copy$tf - 1L, peak_id = copy$idxATAC - 1L, clusters = integer(0),
            cell_numb = nrow(exprs_trans_target))
        if (!is.null(clusters)) {
            # calculate stats for each cluster separately
            exprs_trans_tf_clusters <- Matrix::t(expMatrix_tfs_clusters)
            fclusters <- factor(clusters)
            iclusters <- as.integer(fclusters)
            output_clusters <- fast_wilcox(exprs_x = exprs_trans_target@x, exprs_i = exprs_trans_target@i,
                exprs_p = exprs_trans_target@p, exprs_tf_x = as.logical(exprs_trans_tf_clusters@x),
                exprs_tf_i = exprs_trans_tf_clusters@i, exprs_tf_p = exprs_trans_tf_clusters@p,
                peak_x = peak_trans@x, peak_i = peak_trans@i, peak_p = peak_trans@p,
                target_id = copy$target - 1L, tf_id = copy$tf - 1L, peak_id = copy$idxATAC -
                  1L, clusters = iclusters - 1L, cell_numb = nrow(exprs_trans_target))
            output <- mapply(function(x, y) rbind(x, y), output, output_clusters,
                SIMPLIFY = FALSE)
            # find cluster column indices in the weight matrix
            cluster_col_ind <- iclusters[match(colnames(regulon$weight)[2:ncol(regulon$weight)],
                clusters)]
        }

        AUC <- output$auc
        ties <- output$ties
        n1 <- output$total0
        n2 <- output$total1

        prod <- n1 * n2
        n <- n1 + n2  # technically same across all rows, but we'll just do this for simplicity.
        sigma <- sqrt(prod/12 * (n + 1 - ties/n/(n - 1)))
        mu <- prod/2
        stats <- (AUC - mu)/sigma
        # set z-score to zero if the size of the of the groups is equal to 0
        stats[n1 == 0 | n2 == 0 | sigma == 0] <- 0
        stats[, reg.order] <- stats
        if (is.null(clusters)) {
            regulon$weight[, 1] <- t(stats)
        } else {
            regulon$weight[, c(1, order(cluster_col_ind) + 1)] <- t(stats)
        }
        # Calculate effect size
        n_cells <- ncol(expMatrix)
        if (!is.null(clusters)) {
            # calculate cluster sizes
            # first column in weight matrix is reserved for all cells
            n_cells <- c(n_cells, table(clusters)[colnames(regulon$weight)[2:ncol(regulon$weight)]])
        }
        # transform z-scores to effect size
        regulon$weight <- t(t(regulon$weight)/sqrt(n_cells))
        return(regulon)
    }

    if (method %in% c("corr", "MI")) {
        if (is.null(clusters)) {
            stop("'clusters' argument should be provided for ", method, " method")
        }
        message("calculating average expression across clusters...")

        # define groupings
        groupings <- S4Vectors::DataFrame(cluster = clusters)
        if (!is.null(block_factor)) {
            groupings$block <- colData(expMatrix)[block_factor]
        }

        # compute average expression across clusters and batches
        averages.se.exp <- scuttle::sumCountsAcrossCells(expMatrix, ids = groupings,
            average = TRUE, BPPARAM = BPPARAM)

        # average expression across pseudobulk clusters
        expMatrix <- assays(averages.se.exp)$average

        # remove genes whose expressions are NA for all pseudobulks
        expMatrix <- expMatrix[!Matrix::rowSums(is.na(expMatrix)) == ncol(expMatrix),
            ]


        if (tf_re.merge) {
            averages.se.peak <- scuttle::sumCountsAcrossCells(peakMatrix, ids = groupings,
                average = TRUE, BPPARAM = BPPARAM)

            # average accessibility across pseudobulk clusters
            peakMatrix <- assays(averages.se.peak)$average

        }
        message("computing weights...")
        output_df <- BiocParallel::bplapply(X = seq_len(length(regulon.split)), FUN = use_specific_method,
                                            regulon.split, expMatrix, peakMatrix, tf_re.merge, BPPARAM = BPPARAM,
                                            method = method)
    }

    output_df <- do.call(rbind, output_df)
    output_df[["weight"]][is.na(output_df[["weight"]])] <- 0

    return(output_df)

}

#' @keywords internal

use_specific_method <- function(n, regulon.split, expMatrix, peakMatrix,
                                tf_re.merge, BPPARAM = BPPARAM, method) {
  association_fun <- switch(method, "corr"=stats::cor, "MI" = MI_per_row)
  if (tf_re.merge) {
    tf_re <- expMatrix[regulon.split[[n]]$tf, , drop = FALSE] *
      peakMatrix[as.character(regulon.split[[n]]$idxATAC), , drop = FALSE]
  } else {
    tf_re <- expMatrix[regulon.split[[n]]$tf, , drop = FALSE]
  }
  tg <- expMatrix[regulon.split[[n]]$target, , drop = FALSE]
  regulon.split[[n]]$weight <- mapply(association_fun, as.data.frame(t(tf_re)),
                                      as.data.frame(t(tg)), use = "everything")
  regulon.split[[n]]
}

MI_per_row <- function(tf_re, tg, ...) {
    if (length(unique(tf_re)) == 1 | length(unique(tg)) == 1)
        return(0)
    y2d <- entropy::discretize2d(tf_re, tg, numBins1 = min(10,
        length(unique(tf_re))), numBins2 = min(10, length(unique(tg))))
    entropy::mi.empirical(y2d)
}

