#' Prune regulons for true transcription factor - regulatory elements - target genes relationships
#'
#' @param regulon A dataframe informing the gene regulatory relationship with the ```tf``` column
#' representing transcription factors, ```idxATAC``` corresponding to the index in the peakMatrix and
#'  ```target``` column representing target genes
#' @param expMatrix A SingleCellExperiment object or matrix containing gene expression with
#' genes in the rows and cells in the columns
#' @param peakMatrix A SingleCellExperiment object or matrix containing peak accessibility with
#' peaks in the rows and cells in the columns
#' @param exp_assay String indicating the name of the assay in expMatrix for gene expression
#' @param peak_assay String indicating the name of the assay in peakMatrix for chromatin accessibility
#' @param test String indicating whether `binom` or `chi.sq` test should be performed
#' @param clusters A vector corresponding to the cluster labels of the cells if
#' cluster-specific joint probabilities are also required. If left ```NULL```, joint probabilities
#' are calculated for all cells
#' @param exp_cutoff A scalar indicating the minimum gene expression above which
#' gene is considered active. Default value is 1. Applied to both transcription
#' factors and target genes.
#' @param peak_cutoff A scalar indicating the minimum peak accessibility above which peak is
#' considered open. Default value is 0
#' @param regulon_cutoff A scalar indicating the maximal value for p-value for a tf-idxATAC-target trio
#' to be retained in the pruned regulon.
#' @param p_adj A logical indicating whether p adjustment should be performed
#' @param prune_value String indicating whether to filter regulon based on `pval` or `padj`.
#' @param aggregateCells A logical to indicate whether to aggregate cells into groups determined by cellNum. Use
#' option to overcome data sparsity if needed
#' @param useDim String indicating the name of the dimensionality reduction matrix in expMatrix used for cell aggregation
#' @param cellNum An integer specifying the number of cells per cluster for cell aggregation. Default is 10.
#' @param BPPARAM A BiocParallelParam object specifying whether calculation should be parallelized.
#' Default is set to BiocParallel::MulticoreParam()
#
#' @return A DataFrame of pruned regulons with p-values indicating the probability of independence
#' either for all cells or for individual clusters, z-score statistics for binomial tests or chi-square statistics
#' for chi-square test and q-adjusted values.
#'
#' @details
#' The function prunes the network by performing tests of independence on the observed number of cells
#' jointly expressing transcription factor (TF), regulatory element (RE) and target gene (TG) vs
#' the expected number of cells if TF/RE and TG are independently expressed.
#'
#' In other words, if no regulatory relationship exists, the expected probability of cells expressing all
#' three elements is P(TF, RE) * P(TG), that is, the product of (1) proportion of cells both expressing transcription factor
#' and having accessible corresponding regulatory element, and (2) proportion of cells expressing
#' target gene. The expected number of cells expressing all three elements is therefore n*P(TF, RE)*P(TG),
#' where n is the total number of cells. However, if a TF-RE-TG relationship exists,
#' we expect the observed number of cells jointly having all three elements (TF, RE, TG) to deviate from
#' the expected number of cells predicted from an independent relationship.
#'
#' If the user provides cluster assignment, the tests of independence are performed on a per-cluster basis
#' in addition to providing all cells statistics. This enables pruning by cluster, and thus yields cluster-specific
#' gene regulatory relationships.
#'
#' We implement two tests, the binomial test and the chi-square test.
#'
#' In the binomial test, the expected probability is P(TF, RE) * P(TG), and the number of trials is the number of cells,
#' and the observed successes is the number of cells jointly expressing all three elements.
#'
#' In the chi-square test, the expected probability for having all 3 elements active is also P(TF, RE) * P(TG) and
#' the probability otherwise is 1- P(TF, RE) * P(TG). The observed cell count for the active category is the number of cells
#' jointly expressing all three elements, and the cell count for the inactive category is n - n_triple.
#'
#'
#'
#' @export
#' @import utils SingleCellExperiment
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
#'                   ranges = IRanges(start = seq(from = 1, to = 10000, by = 100), width = 100))
#' peak_counts <- matrix(sample(x = 0:4, size = ncol(gene_sce)*length(peak_gr), replace = TRUE),
#'                      nrow = length(peak_gr), ncol=ncol(gene_sce))
#' peak_sce <- SingleCellExperiment(list(counts = peak_counts), colData = colData(gene_sce))
#' rownames(peak_sce) <- paste0('Peak_',1:100)
#'
#' # create a mock regulon
#' regulon <- data.frame(tf = c(rep('Gene_1',10), rep('Gene_2',10)),
#'                      idxATAC = sample(1:100, 20),
#'                      target = c(paste0('Gene_', sample(3:2000,10)),
#'                                 paste0('Gene_',sample(3:2000,10))))
#'
#' # prune regulon
#' pruned.regulon <- pruneRegulon(expMatrix = gene_sce,
#' exp_assay = 'logcounts', peakMatrix = peak_sce, peak_assay = 'counts',
#' regulon = regulon, clusters = gene_sce$Treatment, regulon_cutoff = 0.5)
#'
#' # add weights with cell aggregation
#' gene_sce <- scater::runPCA(gene_sce)
#' pruned.regulon <- pruneRegulon(expMatrix = gene_sce, exp_assay = 'logcounts',
#' peakMatrix = peak_sce, peak_assay = 'counts', regulon = regulon,
#' clusters = gene_sce$Treatment, regulon_cutoff = 0.5,
#' aggregateCells=TRUE, cellNum=3, useDim = 'PCA')
#'
#' @author Xiaosai Yao, Tomasz Wlodarczyk

pruneRegulon <- function(regulon, expMatrix = NULL, peakMatrix = NULL, exp_assay = "logcounts",
    peak_assay = "PeakMatrix", test = c("chi.sq", "binom"), clusters = NULL, exp_cutoff = 1,
    peak_cutoff = 0, regulon_cutoff = 0.05, p_adj = TRUE, prune_value = "pval", aggregateCells = FALSE,
    useDim = "IterativeLSI_ATAC", cellNum = 10, BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)) {
    tryCatch(ncol(peakMatrix), error = function(cond) stop("'ncol' method should be defined for the 'peakMatrix' object"))
    tryCatch(ncol(expMatrix), error = function(cond) stop("'ncol' method should be defined for the 'expMatrix' object"))
    stopifnot(ncol(peakMatrix) == ncol(expMatrix))
    if (!is.null(clusters)) {
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
    # choose test method
    test <- match.arg(test)
    if(test=="chi.sq" && any(duplicated(rownames(expMatrix)))) stop("Gene names provided in 'expMatrix' are not unique.")
    message("pruning network with ", test, " tests using a regulon cutoff of ", prune_value,
        "<", regulon_cutoff)

    # pseudobulk
    if (aggregateCells) {
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
    }


    # extracting assays from SE
    if (checkmate::test_class(expMatrix, classes = "SummarizedExperiment")) {
        expMatrix <- assay(expMatrix, exp_assay)
    }

    if (checkmate::test_class(peakMatrix, classes = "SummarizedExperiment")) {
        peakMatrix <- assay(peakMatrix, peak_assay)
    }

    expMatrix <- as(expMatrix, "CsparseMatrix")
    peakMatrix <- as(peakMatrix, "CsparseMatrix")



    if (!is.null(peak_cutoff)) {
        prop = sum(peakMatrix > peak_cutoff)/prod(dim(peakMatrix))
        if (prop < 1e-04 | prop > 0.9999)
            warning(sprintf("Strong inbalance between groups after applying cutoff to peakMatrix. Consider %s value of the peak_cutoff",
                c("increasing", "decreasing")[(prop < 1e-04) + 1]))
    }
    if (!is.null(exp_cutoff)) {
        prop = sum(expMatrix > exp_cutoff)/prod(dim(expMatrix))
        if (prop < 1e-04 | prop > 0.9999)
            warning(sprintf("Strong inbalance between groups after applying cutoff to expMatrix. Consider %s value of the exp_cutoff",
                c("increasing", "decreasing")[(prop < 1e-04) + 1]))
    }

    unique_clusters <- c("all", as.character(sort(unique(clusters))))

    # clean up regulons by removing tf and targets not found in regulons
    regulon <- regulon[regulon$tf %in% rownames(expMatrix), ]
    regulon <- regulon[regulon$target %in% rownames(expMatrix), ]
    regulon <- regulon[order(regulon$tf), ]

    message("pruning regulons")
    if (test == "binom") {


        # remove genes not found in regulon
        expMatrix <- expMatrix[which(rownames(expMatrix) %in% unique(c(regulon$tf,
            regulon$target))), ]

        # name peakMatrix
        rownames(peakMatrix) <- seq_len(nrow(peakMatrix))

        # remove peaks not found in regulon
        peakMatrix <- peakMatrix[which(rownames(peakMatrix) %in% unique(regulon$idxATAC)),
            ]

        # binarize peak and expression matrices according to cutoff
        message("binarizing matrices")
        peakMatrix.bi <- binarize_matrix(peakMatrix, peak_cutoff)
        expMatrix.bi <- tfMatrix.bi <- binarize_matrix(expMatrix, exp_cutoff)

        res <- list()
        regulon.split <- split(regulon, regulon$tf)

        # Perform binomial test
        res <- BiocParallel::bplapply(X = seq_len(length(regulon.split)), FUN = binom_bp,
            regulon.split, expMatrix.bi, peakMatrix.bi, tfMatrix.bi, clusters, unique_clusters,
            BPPARAM = BPPARAM)

        res <- do.call("rbind", res)

    } else if (test == "chi.sq") {
        # Perform chi-square test

        if (is.null(clusters)) {
            cluster_id <- factor(integer(ncol(peakMatrix)))
        } else {
            cluster_id <- factor(clusters, levels = as.character(sort(unique(clusters))))
        }
        stats <- countCells(regulon, expMatrix, peakMatrix, cluster_id, peak_cutoff,
            exp_cutoff, clusters)

        if (is.null(clusters)) {
            cluster_freq <- length(cluster_id)
        } else {
            cluster_freq <- c(length(cluster_id), as.numeric(table(cluster_id)))
        }

        cluster_freq <- matrix(nrow = nrow(regulon), ncol = length(cluster_freq),
            cluster_freq, byrow = TRUE)
        peak.prop <- stats$peak/cluster_freq
        target.prop <- stats$target/cluster_freq
        null_probability <- peak.prop * target.prop
        # if p=0 or 1 chi square test would produce NaN values
        test_unavailable_ind <- null_probability%%1 == 0
        res <- chisqTest(k = stats$triple, size = cluster_freq, p = null_probability)
        colnames(res$p) <- sprintf("pval_%s", unique_clusters)
        colnames(res$stat) <- sprintf("stats_%s", unique_clusters)
        res$p[test_unavailable_ind] <- 1
        res$stat[test_unavailable_ind] <- 0
        res <- cbind(res$p, res$stat)

    } else {

        stop("test must be either binom or chi.sq")
    }


    # append test stats to regulon

    pvalue <- res[, grep("^pval_", colnames(res)), drop = FALSE]
    stats <- res[, grep("^stats_", colnames(res)), drop = FALSE]

    colnames(pvalue) <- unique_clusters
    colnames(stats) <- unique_clusters

    regulon.combined <- S4Vectors::DataFrame(regulon, pval = I(pvalue), stats = I(stats))


    # add p-value adjustment

    if (p_adj) {
        "performing multiple testing correction..."

        qvalue <- apply(regulon.combined$pval, 2, function(x) {
            stats::p.adjust(x, method = "holm", n = nrow(regulon.combined))
        })
        colnames(qvalue) <- unique_clusters
        regulon.combined <- S4Vectors::DataFrame(regulon.combined, qval = I(qvalue))
    }

    # prune by p-value
    regulon.prune_value <- regulon.combined[, prune_value, drop = FALSE]
    prune_value_min <- apply(regulon.prune_value, 1, function(x) {
        if (sum(is.na(x)) == length(x))
            Inf else min(x, na.rm = TRUE)
    })
    regulon.combined <- regulon.combined[which(prune_value_min < regulon_cutoff),
        ]

    return(regulon.combined)

}


binarize_matrix <- function(matrix_obj, cutoff = NULL) {
    if (is.null(cutoff))
        cutoff <- Matrix::rowMeans(matrix_obj) else cutoff <- rep(cutoff, nrow(matrix_obj))
    if (is(matrix_obj, "CsparseMatrix")) {
        matrix_obj@x <- as.double(matrix_obj@x > cutoff[matrix_obj@i +
            1])
        matrix_obj
    } else {
        cutoff <- rep(cutoff, ncol(matrix_obj))
        matrix_obj.bi.index <- Matrix::which(matrix_obj >
            cutoff, arr.ind = TRUE)
        matrix_obj <- Matrix::sparseMatrix(x = rep(1,
            nrow(matrix_obj.bi.index)), i = matrix_obj.bi.index[,
            1], j = matrix_obj.bi.index[, 2], dims = dim(matrix_obj),
            dimnames = dimnames(matrix_obj))
    }
}


binom_bp <- function(n, regulon.split, expMatrix.bi, peakMatrix.bi,
    tfMatrix.bi, clusters, unique_clusters, BPPARAM = BPPARAM) {

    full_ncells <- ncol(peakMatrix.bi)
    has_tf <- tfMatrix.bi[regulon.split[[n]]$tf[1], ] == 1
    expMatrix.bi <- expMatrix.bi[regulon.split[[n]]$target, , drop = FALSE]
    expMatrix.tf.bi <- expMatrix.bi[, has_tf, drop = FALSE]
    peakMatrix.bi <- peakMatrix.bi[as.character(regulon.split[[n]]$idxATAC),
        has_tf, drop = FALSE]

    triple.bi <- peakMatrix.bi * expMatrix.tf.bi
    tf_re.bi <- peakMatrix.bi

    res <- list()

    for (selected_cluster in unique_clusters) {
        if (selected_cluster != "all") {
            is_current_cluster <- as.logical(clusters == selected_cluster)
            expCurrent <- as.vector(expMatrix.bi %*% is_current_cluster)
            n_cells <- sum(is_current_cluster)

            #subset is_current_cluster to only cells with tf greater than cutoff
            is_current_cluster <- is_current_cluster[has_tf]
            tf_reCurrent <- as.vector(tf_re.bi %*% is_current_cluster)
            n_triple <- as.vector(triple.bi %*% is_current_cluster)

        } else {
            expCurrent <- Matrix::rowSums(expMatrix.bi)
            n_triple <- Matrix::rowSums(triple.bi)
            tf_reCurrent <- Matrix::rowSums(tf_re.bi)
            n_cells <- full_ncells
        }

        null_probability <- tf_reCurrent * expCurrent/n_cells^2
        p.value <- binomTest(n_triple, n_cells, p = null_probability)
        z_score <- stats::qnorm(p.value/2) * sign(null_probability -
            n_triple/n_cells)

        res[[selected_cluster]] <- cbind(p.value, z_score)
        colnames(res[[selected_cluster]]) <- c(paste0("pval_",
            selected_cluster), paste0("stats_", selected_cluster))
    }

    res <- do.call("cbind", res)


}

binomTest <- function(k, size, p) {
    if (size >= 10000) {
        e1 <- p * size
        e2 <- size - e1
        chi <- (k - e1)^2/e1 + (size - k - e2)^2/e2
        return(stats::pchisq(chi, df = 1, lower.tail = FALSE))
    }

    p.value <- rep_len(1, length(k))
    for (ip in unique(p)) {
        current <- p == ip
        d <- binom_distribution(n = size, p = ip)
        o <- order(d)
        cumsump <- cumsum(d[o])[order(o)]
        p.value[current] <- cumsump[k[current] + 1]
    }
    p.value
}

binom_distribution <- function(n, p) {
    # calculate argmax(binom_probability)
    start_k <- as.integer(round(p * n))
    # calculate the highest probability
    start_dbinom <- stats::dbinom(start_k, n, p)
    # fill distribution going to the left end
    to_start_res <- to_start(start_dbinom, n, start_k, p)
    # fill distribution going to the right end
    to_end_res <- to_end(start_dbinom, n, start_k, p)
    c(rev(to_start_res), to_end_res[2:length(to_end_res)])
}

to_start <- function(start_dbinom, n, start_k, p) .Call("to_start",
    start_dbinom, n, start_k, p)

to_end <- function(start_dbinom, n, start_k, p) .Call("to_end", start_dbinom, n,
    start_k, p)


countCells <- function(regulon, expMatrix, peakMatrix, cluster_id,
    peak_cutoff, exp_cutoff, clusters) {
    if (!is.null(peak_cutoff))
        peak_cutoff <- as.numeric(peak_cutoff)
    if (!is.null(exp_cutoff))
        exp_cutoff <- as.numeric(exp_cutoff)
    peak_id <- regulon$idxATAC
    target_id <- factor(regulon$target, levels = rownames(expMatrix))
    tf_id <- factor(regulon$tf, levels = rownames(expMatrix))

    p_o <- order(peak_id, tf_id, target_id)
    p_peak_id <- as.integer(peak_id[p_o]) - 1L
    p_target_id <- as.integer(target_id[p_o]) - 1L
    p_tf_id <- as.integer(tf_id[p_o]) - 1L

    t_o <- order(target_id)
    t_target_id <- as.integer(target_id[t_o]) - 1L

    cluster_id2 <- as.integer(cluster_id) - 1L
    if (!is.null(exp_cutoff) & !is.null(peak_cutoff)) {
        exp_cutoff_mat <- matrix(exp_cutoff, nrow = nrow(expMatrix),
            ncol = nlevels(cluster_id))
        peak_cutoff_mat <- matrix(peak_cutoff, nrow = nrow(peakMatrix),
            ncol = nlevels(cluster_id))
        stats <- fast_chisq(peak_ordered = p_peak_id, tf_by_peak = p_tf_id,
            target_by_peak = p_target_id, target_ordered = t_target_id,
            npeaks = nrow(peakMatrix), peakmat_x = peakMatrix@x,
            peakmat_i = peakMatrix@i, peakmat_p = peakMatrix@p,
            peak_cutoff = peak_cutoff_mat, ngenes = nrow(expMatrix),
            expmat_x = expMatrix@x, expmat_i = expMatrix@i,
            expmat_p = expMatrix@p, exp_cutoff = exp_cutoff_mat,
            nclusters = nlevels(cluster_id), clusters = cluster_id2)
        if (!is.null(clusters)) {
            stats$peak <- cbind(rowSums(stats$peak), stats$peak)
            stats$triple <- cbind(rowSums(stats$triple),
                stats$triple)
            stats$target <- cbind(rowSums(stats$target),
                stats$target)
        }

    } else {
        if (is.null(peak_cutoff)) {
            peak_cutoff_mat <- matrix(Matrix::rowMeans(peakMatrix),
                nrow = nrow(peakMatrix), ncol = nlevels(cluster_id))
        } else {
            peak_cutoff_mat <- matrix(peak_cutoff, nrow = nrow(peakMatrix),
                ncol = nlevels(cluster_id))
        }
        if (is.null(exp_cutoff)) {
            exp_cutoff_mat <- matrix(Matrix::rowMeans(expMatrix),
                nrow = nrow(expMatrix), ncol = nlevels(cluster_id))
        } else {
            exp_cutoff_mat <- matrix(exp_cutoff, nrow = nrow(expMatrix),
                ncol = nlevels(cluster_id))
        }
        stats <- fast_chisq(peak_ordered = p_peak_id, tf_by_peak = p_tf_id,
            target_by_peak = p_target_id, target_ordered = t_target_id,
            npeaks = nrow(peakMatrix), peakmat_x = peakMatrix@x,
            peakmat_i = peakMatrix@i, peakmat_p = peakMatrix@p,
            peak_cutoff = peak_cutoff_mat, ngenes = nrow(expMatrix),
            expmat_x = expMatrix@x, expmat_i = expMatrix@i,
            expmat_p = expMatrix@p, exp_cutoff = exp_cutoff_mat,
            nclusters = 1L, clusters = rep(0L, ncol(expMatrix)))

        if (!is.null(clusters)) {
            if (is.null(exp_cutoff)) {
                for (cluster in unique(cluster_id2)) {
                  cluster_ind <- which(cluster_id2 == cluster)
                  exp_cutoff_mat[, as.numeric(cluster) +
                    1] <- Matrix::rowMeans(expMatrix[, cluster_ind,
                    drop = FALSE])
                }
            }
            if (is.null(peak_cutoff)) {
                for (cluster in unique(cluster_id2)) {
                  cluster_ind <- which(cluster_id2 == cluster)
                  peak_cutoff_mat[, as.numeric(cluster) +
                    1] <- Matrix::rowMeans(peakMatrix[, cluster_ind,
                    drop = FALSE])
                }
            }
            stats_clusters <- fast_chisq(peak_ordered = p_peak_id,
                tf_by_peak = p_tf_id, target_by_peak = p_target_id,
                target_ordered = t_target_id, npeaks = nrow(peakMatrix),
                peakmat_x = peakMatrix@x, peakmat_i = peakMatrix@i,
                peakmat_p = peakMatrix@p, peak_cutoff = peak_cutoff_mat,
                ngenes = nrow(expMatrix), expmat_x = expMatrix@x,
                expmat_i = expMatrix@i, expmat_p = expMatrix@p,
                exp_cutoff = exp_cutoff_mat, nclusters = nlevels(cluster_id),
                clusters = cluster_id2)

            stats$peak <- cbind(stats$peak, stats_clusters$peak)
            stats$triple <- cbind(stats$triple, stats_clusters$triple)
            stats$target <- cbind(stats$target, stats_clusters$target)

        }
    }
    stats$peak[p_o, ] <- stats$peak
    stats$triple[p_o, ] <- stats$triple
    stats$target[t_o, ] <- stats$target
    stats
}

chisqTest <- function(k, size, p) {
    e1 <- p * size
    e2 <- size - e1
    chi <- (k - e1)^2/e1 + (size - k - e2)^2/e2
    list(p = stats::pchisq(chi, df = 1, lower.tail = FALSE), stats = chi)
}

