#' Calculate the per cell activity of master regulators based on a regulon
#'
#' @param expMatrix A SingleCellExperiment object containing gene expression information with rows representing genes and columns represent cells.
#' Rownames (either gene symbols or geneID) must be consistent with the naming convention in the regulon.
#' @param exp_assay String specifying the name of the assay to be retrieved from the SingleCellExperiment object. Set to
#' 'logcounts' as the default
#' @param regulon  A DataFrame object consisting of tf (regulator) and target in the column names, with additional columns
#' indicating degree of association between tf and target such as 'mor' or 'corr' obtained from `addWeights`.
#' @param normalize Logical indicating whether row means should be subtracted from expression matrix. default is FALSE
#' @param mode String indicating the name of column to be used as the weights
#' @param method String indicating the method for calculating activity. Available methods are `weightedMean` or `aucell`
#' @param genesets A list of genesets. Each list element can be a dataframe with the first column indicating the genes and second column indicating the weights.
#' Alternatively, each list element is a character vector corresponding to the genes in the geneset. A feature set collection in the form of CompressedSplitDataFrameList that
#' contains genes in the first column and weights in the second column. See details
#' @param clusters A vector indicating cluster assignment
#' @param FUN function to aggregate the weights
#' @param ncore Integer specifying the number of cores to be used in AUCell
#' @param BPPARAM A BiocParallelParam object specifying whether summation should be parallelized. Use BiocParallel::SerialParam() for
#' serial evaluation and use BiocParallel::MulticoreParam() for parallel evaluation
#' @return A matrix of inferred transcription factor (row) activities in single cells (columns)
#' @export
#' @import methods utils
#' @details
#' This function calculates activity score from a regulon that is a DataFrame consisting of a tf column,
#' a target column and a weight column. Alternatively, instead of a regulon, this function also accepts weighted
#' signature sets where each gene set or signature is a data frame or unweighted signature sets where each gene set is a character vector.
#' The user has the option of computing signature score by weighted mean of target gene expression or
#' the relative ranking of the target genes computed by AUCell.
#'
#' @examples
#' # create a mock singleCellExperiment object for gene expMatrixession matrix
#' set.seed(1000)
#' gene_sce <- scuttle::mockSCE()
#' gene_sce <- scuttle::logNormCounts(gene_sce)
#' rownames(gene_sce) <- paste0('Gene_',1:2000)
#'
#' # create a mock singleCellExperiment object for peak matrix
#' peak_gr <- GRanges(seqnames = 'chr1',
#'                    ranges = IRanges(start = seq(from = 1, to = 10000, by = 100), width = 100))
#' peak_counts <- matrix(sample(x = 0:4, size = ncol(gene_sce)*length(peak_gr), replace = TRUE),
#'                       nrow = length(peak_gr), ncol=ncol(gene_sce))
#' peak_sce <- SingleCellExperiment(list(counts = peak_counts), colData = colData(gene_sce))
#' rownames(peak_sce) <- paste0('Peak_',1:100)
#'
#' # create a mock regulon
#' regulon <- data.frame(tf = c(rep('Gene_1',10), rep('Gene_2',10)),
#'                       idxATAC = sample(1:100, 20),
#'                       target = c(paste0('Gene_', sample(3:2000,10)),
#'                                  paste0('Gene_',sample(3:2000,10))))
#'
#' #  # prune regulon
#' pruned.regulon <- pruneRegulon(expMatrix = gene_sce,
#'                                exp_assay = 'logcounts',
#'                                peakMatrix = peak_sce,
#'                                peak_assay = 'counts',
#'                                regulon = regulon,
#'                                clusters = gene_sce$Treatment,
#'                                regulon_cutoff = 0.5,
#'                                p_adj = TRUE)
#'
#' regulon.w <- addWeights(regulon = regulon,
#'                         expMatrix = gene_sce,
#'                         clusters = gene_sce$Treatment,
#'                         exp_assay = 'logcounts',
#'                         min_targets = 5,
#'                         method = 'corr')
#'
#' # calculate activity
#' activity <- calculateActivity(expMatrix = gene_sce,
#'                               regulon = regulon.w,
#'                               exp_assay = 'logcounts')
#'
#' # calculate cluster-specific activity if cluster-specific weights are supplied
#' regulon.w$weight <- matrix(runif(nrow(regulon.w)*2, -1,1), nrow(regulon.w),2)
#' colnames(regulon.w$weight) <- c('treat1','treat2')
#'
#' activity.cluster <- calculateActivity(gene_sce,
#' regulon = regulon.w, clusters = gene_sce$Treatment,
#' exp_assay = 'logcounts', FUN = 'mean')
#'
#' # compute signature scores from weighted genesets
#' weighted_genesets <- list(set1 = data.frame(genes = c('Gene_1', 'Gene_2', 'Gene_3'),
#' weights = c(1,2,3)), set2 = data.frame(genes = c('Gene_4', 'Gene_5', 'Gene_6'), weights = c(4,5,6)))
#'
#' activity <- calculateActivity(gene_sce, genesets = weighted_genesets)
#'
#' # compute signature scores from unweighted genesets
#' unweighted_genesets <- list(set1 = c('Gene_1', 'Gene_2', 'Gene_3'),
#'                             set2 = c('Gene_4', 'Gene_5', 'Gene_6'))
#' activity <- calculateActivity(gene_sce, genesets = unweighted_genesets)
#'

#' @author Xiaosai Yao, Shang-yang Chen

calculateActivity <- function(expMatrix = NULL, exp_assay = "logcounts", regulon = NULL,
    normalize = FALSE, mode = "weight", method = c("weightedmean", "aucell"), genesets = NULL,
    clusters = NULL, FUN = c("mean", "sum"), ncore = 1, BPPARAM = BiocParallel::SerialParam()) {
    if (!is.null(regulon))
        checkmate::assertMultiClass(regulon, c("data.frame", "DFrame"))
    if (!is.null(regulon) && !mode %in% colnames(regulon))
        stop("No such column in the regulon: ", mode)
    .validate_input_sce(expMatrix, exp_assay, env=environment())
    if(!is.null(clusters)) .validate_clusters(clusters, expMatrix)
    method <- tolower(method)
    method <- match.arg(method)
    FUN <- match.arg(FUN)

    regulon <- S4Vectors::DataFrame(regulon)
    # convert expMatrix to CsparseMatrix
    expMatrix <- assay(expMatrix, exp_assay)
    expMatrix <- as(expMatrix, "CsparseMatrix")

    # convert genesets to regulon
    if (!is.null(genesets)) {
        if (is.list(genesets) | is(genesets, "CompressedSplitDFrameList")) {
            regulon <- genesets2regulon(genesets)
        } else {
            stop("genesets should be a list of data frames or character vectors")
        }
    }

    # remove rows with zero weight
    if (is.matrix(regulon[[mode]]) & !is.null(clusters)) {
        #  check all cluster columns
        regulon <- regulon[which(apply(regulon[[mode]], 1, function(x) any(x[2:length(x)] !=
            0))), ,drop=FALSE]
    } else if (is.matrix(regulon[[mode]]) & is.null(clusters)) {
        # check all cells column
        regulon <- regulon[which(apply(regulon[[mode]], 1, function(x) x[1] != 0)), ,
                           drop=FALSE]
    } else {
        regulon <- regulon[which(regulon[, mode] != 0), , drop=FALSE]
    }
    if (nrow(regulon) == 0) {
        warning("No non-zero weight in the regulon")
        return(NULL)
    }

    # check that rownames match regulon
    fraction_genes <- length(which(regulon$target %in% rownames(expMatrix)))/length(regulon$target)
    if (fraction_genes < 0.01) {
        warning("Less than 1% of target genes in the regulon are found in expression matrix. Check rownames of gene expression matrix ")
    }


    # remove genes in regulons not found in expMatrix
    regulon <- regulon[which(regulon$target %in% rownames(expMatrix)), , drop = FALSE]

    # calculate activity
    if (method == "weightedmean") {
        message("calculating TF activity from regulon using ", method)


        # aggregate weights across the same tf-target pairs
        if (is.null(clusters) & length(regulon[1, mode]) > 1) {
            warning("The ", mode, " column contains multiple subcolumns but no cluster information was provided. Using first column to compute activity...")
            regulon[, mode] <- regulon[, mode][, 1]
        }

        # aggregate weights across the same tf-target pairs
        if (!is.null(clusters)) {
            regulon[, mode] <- I(as.matrix(regulon[, mode]))
        }

        message("aggregating regulons...")
        aggregated.regulon <- aggregateMatrix(regulon, mode, FUN)


        # create tf x target matrix of weights
        message("creating weight matrix...")
        tf_target_mat <- createTfTgMat(aggregated.regulon, mode, clusters = clusters)


        # if cluster information is provided and if there are cluster-specific weights provided,
        # compute total activity by summation of cluster-specific activity
        if (is.null(clusters)) {
            # if no cluster information is provided, calculate activity for all cells

            message("calculating activity scores...")
            # cross product of expMatrix and tf_target matrix
            score.combine <- calculateScore(expMatrix, tf_target_mat)

            # need to normalize
            if (normalize) {
                message("normalize by mean...")
                meanExpr <- Matrix::rowMeans(expMatrix[rownames(tf_target_mat), ,drop=FALSE])
                mean_activity <- meanExpr %*% tf_target_mat
                score.combine <- sweep(score.combine, 2, mean_activity, "-")
            }
            message("normalize by the number of targets...")
            #normalize by number of targets
            freq <- calculateFrequency(regulon = aggregated.regulon, mode = mode)
            score.combine <- normalizeByFrequency(score.combine, freq)

        } else if (!is.null(clusters)) {

            # Calculate the number of targets per cluster
            # freq is a table of tf x clusters and the elements represent the number of targets per tf
            freq <- initiateMatCluster(clusters, nrow = length(unique(regulon$tf)),
                value = 1)
            rownames(freq) <- unique(regulon$tf)

            message("calculating frequency...")
            freq <- calculateFrequency(freq, aggregated.regulon, mode = mode)

            # Calculating scores
            score.combine <- matrix(0, nrow = ncol(expMatrix), ncol = length(unique(regulon$tf)))
            rownames(score.combine) <- colnames(expMatrix)
            colnames(score.combine) <- colnames(tf_target_mat[[1]])
            score.combine <- as(score.combine, "CsparseMatrix")

            message("calculating activity scores...")
            score.combine <- calculateScore(expMatrix, tf_target_mat, clusters = clusters,
                score.combine)


            # if normalize gene expression (taking the mean across all cells)
            if (normalize) {
                message("normalize by mean...")
                meanExpr <- Matrix::rowMeans(expMatrix[rownames(tf_target_mat[[1]]),
                  ])
                for (cluster in sort(unique(clusters))) {
                  # calculate cluster-specific mean
                  mean_activity <- meanExpr %*% tf_target_mat[[cluster]]
                  score.combine[clusters == cluster, ] <- sweep(score.combine[clusters ==
                    cluster, ,drop=FALSE], 2, mean_activity, "-")
                }
            }

            message("normalize by number of targets...")
            # normalize by the number of target genes
            score.combine <- normalizeByFrequency(score.combine, freq, clusters = clusters)

        }
        score.combine <- Matrix::t(score.combine)

    } else if (method == "aucell") {
        message("calculating TF activity from regulon using ", method)
        geneSets <- split(regulon$target, regulon$tf)
        message("ranking cells...")
        cells_rankings <- AUCell::AUCell_buildRankings(expMatrix, splitByBlocks = TRUE,
            plotStats = FALSE, BPPARAM = BPPARAM)
        message("calculating AUC...")
        cells_AUC <- AUCell::AUCell_calcAUC(geneSets, rankings = cells_rankings,
            nCores = ncore)

        score.combine <- AUCell::getAUC(cells_AUC)
    }
    return(score.combine)
}

genesets2regulon <- function(genesets) {
    regulon <- list()
    for (i in seq_len(length(genesets))) {
        if (is(genesets[[i]], "DFrame") | is(genesets[[i]], "data.frame")) {
            regulon[[i]] <- S4Vectors::DataFrame(tf = names(genesets)[i],
                target = genesets[[i]][, 1], weight = genesets[[i]][, 2])
        } else if (is.vector(genesets[[i]])) {
            regulon[[i]] <- S4Vectors::DataFrame(tf = names(genesets)[i],
                target = genesets[[i]], weight = 1)
        }
    }

    regulon <- do.call(rbind, as.list(regulon))
}


genesets2regulon <- function(genesets) {
  regulon <- lapply(seq_along(genesets), function(i){
    if(is(genesets[[i]], "DFrame") | is(genesets[[i]], "data.frame")){
      S4Vectors::DataFrame(tf = names(genesets)[i],
                           target = genesets[[i]][, 1], weight = genesets[[i]][, 2])
    }
    else if (is.vector(genesets[[i]])) {
      S4Vectors::DataFrame(tf = names(genesets)[i], target = genesets[[i]], weight = 1)
    }
  })
  regulon <- do.call(rbind, as.list(regulon))
}


createTfTgMat <- function(regulon, mode, clusters = NULL) {

    regulon$tfidx <- as.numeric(as.factor(regulon$tf))
    regulon$targetidx <- as.numeric(as.factor(regulon$target))

    n_target <- length(unique(regulon$target))
    n_tf <- length(unique(regulon$tf))

    if (is.null(clusters)) {
        tf_target_mat <- Matrix::sparseMatrix(x = as.vector(regulon[, mode]), i = regulon$targetidx,
            j = regulon$tfidx, dims = c(n_target, n_tf))

        colnames(tf_target_mat) <- levels(as.factor(regulon$tf))
        rownames(tf_target_mat) <- levels(as.factor(regulon$target))
        tf_target_mat[is.na(tf_target_mat)] <- 0

    } else if (!is.null(clusters)) {
        tf_target_mat <- list()
        for (cluster in unique(clusters)) {
            tf_target_mat[[cluster]] <- Matrix::sparseMatrix(x = as.vector(regulon[,
                mode][, cluster]), i = regulon$targetidx, j = regulon$tfidx, dims = c(n_target,
                n_tf))

            colnames(tf_target_mat[[cluster]]) <- levels(as.factor(regulon$tf))
            rownames(tf_target_mat[[cluster]]) <- levels(as.factor(regulon$target))
            tf_target_mat[[cluster]][is.na(tf_target_mat[[cluster]])] <- 0
        }
    }

    tf_target_mat
}



calculateScore <- function(expMatrix, tf_target_mat, clusters = NULL,
    score.combine = NULL) {
    if (is.null(clusters)) {
        score.combine <- Matrix::t(expMatrix[rownames(tf_target_mat),
            , drop = FALSE]) %*% tf_target_mat
        rownames(score.combine) <- colnames(expMatrix)
        colnames(score.combine) <- colnames(tf_target_mat)

    } else {

        for (cluster in sort(unique(clusters))) {
            expr_data <- expMatrix[rownames(tf_target_mat[[cluster]]),
                clusters == cluster, drop = FALSE]
            score.combine[clusters == cluster, ] <- Matrix::t(expr_data) %*%
                tf_target_mat[[cluster]]
        }
    }
    score.combine
}

calculateFrequency <- function(freq = NULL, regulon, mode) {
    if (any(is.null(ncol(regulon[, mode])), ncol(regulon[, mode]) == 1)) {
        freq <- table(regulon$tf[as.vector(regulon[, mode]) != 0])
        freq[freq == 0 | is.na(freq)] <- 1
    } else {
        for (cluster in colnames(regulon[, mode])) {
            freq_counts <- table(regulon$tf[which(regulon[, mode][, cluster] !=
                0)])
            freq_counts[freq_counts == 0] <- 1
            freq[names(freq_counts), cluster] <- freq_counts
        }
    }

    freq
}

normalizeByFrequency <- function(score.combine,
    freq, clusters = NULL) {
    if (is.null(clusters)) {
        score.combine[, names(freq)] <- sweep(score.combine[,
            names(freq), drop = FALSE], 2, freq,
            "/")
    } else {
        for (cluster in unique(clusters)) {
            score.combine[clusters == cluster,
                rownames(freq)] <- sweep(score.combine[clusters ==
                cluster, rownames(freq), drop = FALSE],
                2, freq[, cluster], "/")
        }
    }
    score.combine
}
