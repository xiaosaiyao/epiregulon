#' Establish peak to gene links based on correlations between ATAC-seq peaks and RNA-seq genes
#'
#' @param peakMatrix A SingleCellExperiment object containing counts of chromatin accessibility at each peak region or genomic bin from scATAC-seq.
#' `rowRanges` should contain genomic positions of the peaks in the form of `GRanges`.
#' @param expMatrix A SingleCellExperiment object containing gene expression counts from scRNA-seq. `rowRanges` should contain genomic positions of
#' the genes in the form of `GRanges`. `rowData` should contain a column of gene symbols with column name matching the `gene_symbol` argument.
#' @param reducedDim A matrix of dimension reduced values
#' @param cor_cutoff A numeric scalar to specify the correlation cutoff between ATAC-seq peaks and RNA-seq genes to assign peak to gene links.
#'  Default correlation cutoff is 0.5.
#' @param useDim String specifying the name of the reduced dimension matrix supplied by the user
#' @param cellNum An integer to specify the number of cells to include in each K-means cluster
#' @param maxDist An integer to specify the base pair extension from transcription start start for overlap with peak regions
#' @param exp_assay String indicating the name of the assay in expMatrix for gene expression
#' @param peak_assay String indicating the name of the assay in peakMatrix for chromatin accessibility
#' @param gene_symbol String indicating the column name in the rowData of expMatrix that corresponds to gene symbol
#' @param clusters A vector corresponding to the cluster labels for calculation of correlations within each cluster. If left NULL, correlation is calculated across
#' all clusters. See details for the use of clusters
#' @param cor_method String indicating which correlation coefficient is to be computed. One of 'pearson' (default), 'kendall', or 'spearman'.
#' @param BPPARAM A BiocParallelParam object specifying whether summation should be parallelized. Use BiocParallel::SerialParam() for
#' serial evaluation and use BiocParallel::MulticoreParam() for parallel evaluation
#'
#' @return A DataFrame of Peak to Gene correlation
#' @details Cluster information is sometimes helpful to avoid the [Simpsons's paradox](https://en.wikipedia.org/wiki/Simpson%27s_paradox) in which baseline differences
#' between cell lines or cell types can create artificial or even inverse correlations between peak accessibility and gene expression. If Cluster information is provided,
#' correlation is performed within cell aggregates of each cluster.
#' @import SummarizedExperiment SingleCellExperiment GenomicRanges
#' @export
#'
#' @examples
#' # create a mock singleCellExperiment object for gene expression matrix
#' set.seed(1000)
#' gene_sce <- scuttle::mockSCE()
#' gene_sce <- scuttle::logNormCounts(gene_sce)
#' gene_gr <- GRanges(seqnames = Rle(c('chr1', 'chr2', 'chr3','chr4'), nrow(gene_sce)/4),
#'                    ranges = IRanges(start = seq(from = 1, length.out=nrow(gene_sce), by = 1000),
#'                    width = 100))
#' rownames(gene_sce) <- rownames(gene_sce)
#' gene_gr$name <- rownames(gene_sce)
#' rowRanges(gene_sce) <- gene_gr
#'
#' # create a mock singleCellExperiment object for peak matrix
#' peak_gr <- GRanges(seqnames = 'chr1',
#'                    ranges = IRanges(start = seq(from = 1, to = 10000, by = 1000), width = 100))
#' peak_counts <- matrix(sample(x = 0:4, size = ncol(gene_sce)*length(peak_gr), replace = TRUE),
#'                       nrow = length(peak_gr), ncol=ncol(gene_sce))
#' peak_sce <- SingleCellExperiment(list(counts = peak_counts), colData = colData(gene_sce))
#' rowRanges(peak_sce) <- peak_gr
#' rownames(peak_sce) <- paste0('peak',1:10)

#' # create a mock reducedDim matrix
#' reducedDim_mat <- matrix(runif(ncol(gene_sce)*50, min = 0, max = 1), nrow = ncol(gene_sce), 50)
#' p2g <- calculateP2G(peakMatrix = peak_sce, expMatrix = gene_sce, reducedDim = reducedDim_mat,
#'                     cellNum = 20, clusters = gene_sce$Treatment)
#' @author Xiaosai Yao, Shang-yang Chen

calculateP2G <- function(peakMatrix = NULL, expMatrix = NULL, reducedDim = NULL,
    useDim = "IterativeLSI", maxDist = 250000,
    cor_cutoff = 0.5, cellNum = 100, exp_assay = "logcounts", peak_assay = "counts",
    gene_symbol = "name", clusters = NULL, cor_method = c("pearson", "kendall", "spearman"),
    BPPARAM = BiocParallel::SerialParam()) {


    if (!is.null(peakMatrix) & !is.null(expMatrix) & !is.null(reducedDim)) {

        writeLines("Using epiregulon to compute peak to gene links...")

        cor_method <- match.arg(cor_method)

        checkmate::assert_class(rowRanges(peakMatrix), "GRanges")

        if (length(rowRanges(peakMatrix)) == 0) {
            stop("peakMatrix should contain non-empty rowRanges")
        }

        checkmate::assert_class(rowRanges(expMatrix), "GRanges")

        if (length(rowRanges(expMatrix)) == 0) {
            stop("expMatrix should contain non-empty rowRanges")
        }

        if (!gene_symbol %in% colnames(rowData(expMatrix))) {
            stop("colData of expMatrix does not contain ", gene_symbol)
        }
        if (cellNum > ncol(expMatrix)) stop("The value of 'cellNum' parameter cannot be greater than the total number of cells")

        # Package expression matrix and peak matrix into a single sce
        sce <- combineSCE(expMatrix, exp_assay, peakMatrix, peak_assay, reducedDim,
            useDim)

        message("performing k means clustering to form metacells")

        kNum <- trunc(ncol(sce)/cellNum)
        kclusters <- scran::clusterCells(sce, use.dimred = useDim, BLUSPARAM = bluster::KmeansParam(centers = kNum,
            iter.max = 5000))
        kclusters <- as.character(kclusters)
        cluster_numb_warning <- length(unique(kclusters)) < 5

        # aggregate sce by k-means clusters
        sce_grouped <- applySCE(sce, scuttle::aggregateAcrossCells, ids = kclusters,
            statistics = "mean", BPPARAM = BPPARAM)

        # some sces have strand information in metadata that conflicts with genomic ranges
        mcols(expMatrix)$strand <- NULL

        # keep track of original ATAC and expression indices
        rowData(sce_grouped)$old.idxRNA <- seq_len(nrow(sce_grouped))
        rowData(altExp(sce_grouped))$old.idxATAC <- seq_len(nrow(altExp(sce_grouped)))

        # remove genes and peaks that are equal to 0
        sce_grouped <- sce_grouped[which(rowSums(assay(sce_grouped)) != 0), ]
        altExp(sce_grouped) <- altExp(sce_grouped)[which(rowSums(assay(altExp(sce_grouped),
            "counts")) != 0), ]

        # extract gene expression and peak matrix
        expGroupMatrix <- assay(sce_grouped, "counts")
        peakGroupMatrix <- assay(altExp(sce_grouped), "counts")


        # get gene information
        geneSet <- rowRanges(sce_grouped)
        geneStart <- promoters(geneSet)

        # get peak range information
        peakSet <- rowRanges(altExp(sce_grouped))

        # find overlap after resizing
        o <- S4Vectors::DataFrame(findOverlaps(resize(geneStart, maxDist, "center"),
            peakSet, ignore.strand = TRUE))



        #Get Distance from Fixed point A B
        o$distance <- distance(geneStart[o[, 1]], peakSet[o[, 2]])
        colnames(o) <- c("RNA", "ATAC", "distance")


        # add old idxRNA and idxATAC
        o$old.idxRNA <- rowData(sce_grouped)[o[, 1], "old.idxRNA"]
        o$old.idxATAC <- rowData(altExp(sce_grouped))[o[, 2], "old.idxATAC"]

        #add metadata to o
        o$Gene <- rowData(sce_grouped)[o[, 1], gene_symbol]
        o$chr <- as.character(seqnames(rowRanges(altExp(sce_grouped))[o[, 2]]))
        o$start <- GenomicRanges::start(rowRanges(altExp(sce_grouped))[o[, 2], ])
        o$end <- GenomicRanges::end(rowRanges(altExp(sce_grouped))[o[, 2], ])

        # Calculate correlation
        expCorMatrix <- expGroupMatrix[as.integer(o$RNA), ]
        peakCorMatrix <- peakGroupMatrix[as.integer(o$ATAC), ]

        writeLines("Computing correlation")

        # if a cluster is named 'all', replace it to distinguish from all cells
        clusters <- renameCluster(clusters)

        unique_clusters <- sort(unique(clusters))
        if(any(unique_clusters=="")) stop("Some of the culster lables are empty strings.")

        o$Correlation <- initiateMatCluster(clusters, nrow = nrow(expCorMatrix))
        o$Correlation[, "all"] <- mapply(stats::cor, as.data.frame(t(expCorMatrix)),
            as.data.frame(t(peakCorMatrix)), MoreArgs = list(method = cor_method))

        # compute correlation within each cluster
        if (!is.null(clusters)) {
            # composition of kcluster
            cluster_composition <- table(clusters, kclusters)
            cluster_composition <- sweep(cluster_composition, 2, STATS = colSums(cluster_composition),
                FUN = "/")
            for (cluster in unique_clusters) {
                clusters_idx <- colnames(cluster_composition)[cluster_composition[cluster,
                  ] >= 1/length(unique_clusters)]
                if(length(clusters_idx)<5) cluster_numb_warning <- TRUE
                if(length(clusters_idx)<3) o$Correlation[, cluster] <- NA
                else{
                  o$Correlation[, cluster] <- mapply(stats::cor,
                                                     as.data.frame(t(expCorMatrix[,clusters_idx])),
                                                     as.data.frame(t(peakCorMatrix[, clusters_idx])))
                }
            }
        }
        if(cluster_numb_warning) {
          suggested_numb <- sqrt(ncol(sce))
          if(!is.null(clusters)) suggested_numb <- suggested_numb/length(unique_clusters)
          if(round(suggested_numb)<10)
            warning("The number of aggregated cells in user-specified cluster is low. Consider providing lesser number of clusters")
          else
            warning(sprintf("The number of aggregated cells in user-specified cluster is low. Consider dropping cells from small clusters or changing cellNum parameter to %d", round(suggested_numb)))
          }

        p2g_merged <- o[, c("old.idxATAC", "chr", "start", "end", "old.idxRNA", "Gene",
            "Correlation", "distance")]
        colnames(p2g_merged) <- c("idxATAC", "chr", "start", "end", "idxRNA", "target",
            "Correlation", "distance")

        correlation_max <- apply(p2g_merged$Correlation, 1, max, na.rm = TRUE)
        p2g_merged <- p2g_merged[correlation_max > cor_cutoff, , drop = FALSE]



    } else {
        stop("Input obj must be either an 'ArchR' path or all 3 matrices: gene expression, chromatin accessibility and dimensionality reduction")

    }
    p2g_merged <- p2g_merged[order(p2g_merged$idxATAC, p2g_merged$idxRNA), , drop = FALSE]
    return(p2g_merged)

}


combineSCE <- function(expMatrix, exp_assay, peakMatrix, peak_assay, reducedDim,
    useDim) {

    # convert expMatrix and peakMatrix in case they weren't already so
    expMatrix <- as(expMatrix, "SingleCellExperiment")
    peakMatrix <- as(peakMatrix, "SingleCellExperiment")


    sce <- SingleCellExperiment(list(counts = assay(expMatrix, exp_assay)), altExps = list(peakMatrix = SingleCellExperiment(list(counts = assay(peakMatrix,
        peak_assay)))))

    rowRanges(sce) <- rowRanges(expMatrix)
    rowRanges(altExp(sce)) <- rowRanges(peakMatrix)

    # add reduced dimension information to sce object
    reducedDim(sce, useDim) <- reducedDim

    sce

}
