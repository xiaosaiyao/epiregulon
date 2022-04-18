#' A function to establish peak to gene links based on correlations between ATAC-seq peaks and RNA-seq genes
#'
#' @param peakMatrix A SingleCellExperiment object containing counts of chromatin accessibility at each peak region or genomic bin from scATAC-seq
#' @param expMatrix A SingleCellExperiment object containing gene expression counts from scRNA-seq
#' @param reducedDim A matrix of dimension reduced values, for example derived from IterativeLSI algorithm of ArchR
#' @param ArchR_path String specifying the path to an ArchR project if ArchR's implementation of addPeak2GeneLinks is desired
#' @param reducedDimName String specifying the name of the reduced dimension matrix
#' @param cor_cutoff A numeric scalar to specify the correlation cutoff between ATAC-seq peaks and RNA-seq genes to assign peak to gene links.
#'  Default correlation cutoff is 0.5.
#' @param useDim String specifying which dimensional reduction representation in the ArchR project to use
#' @param useMatrix String specifying which data matrix in the ArchR project to use
#' @param cellNum An integer to specify the number of cells to include in each K-means cluster. Default is 200 cells.
#' @param seed An integer scalar to specify the seed for K-means clustering
#' @param ... other parameters to pass to addPeak2GeneLinks from ArchR package
#'
#' @return A Peak2Gene correlation datafrane
#' @export
#'
#' @examples
#' # create a mock singleCellExperiment object for gene expression matrix
#' set.seed(1000)
#' gene_sce <- scuttle::mockSCE()
#' gene_sce <- scuttle::logNormCounts(gene_sce)
#' gene_gr <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr3","chr4"), nrow(gene_sce)/4),
#'                    ranges = IRanges(start = seq(from = 1, length.out=nrow(gene_sce), by = 1000), width = 100))
#' rownames(gene_sce) <- rownames(gene_sce)
#' gene_gr$name <- rownames(gene_sce)
#' rowRanges(gene_sce) <- gene_gr
#'
#' # create a mock singleCellExperiment object for peak matrix
#' peak_gr <- GRanges(seqnames = "chr1",
#'                    ranges = IRanges(start = seq(from = 1, to = 10000, by = 1000), width = 100))
#' peak_counts <- matrix(sample(x = 0:4, size = ncol(gene_sce)*length(peak_gr), replace = TRUE),
#'                       nrow = length(peak_gr), ncol=ncol(gene_sce))
#' peak_sce <- SingleCellExperiment(list(counts = peak_counts), colData = colData(gene_sce))
#' rowRanges(peak_sce) <- peak_gr
#' rownames(peak_sce) <- paste0("peak",1:10)

#' # create a mock reducedDim matrix
#' reducedDim_mat <- matrix(runif(ncol(gene_sce)*50, min = 0, max = 1), nrow = ncol(gene_sce), 50)
#' p2g <- calculateP2G(peakMatrix = peak_sce, expMatrix = gene_sce, reducedDim = reducedDim_mat,
#'                     cellNum = 20)

calculateP2G = function(peakMatrix = NULL,
                        expMatrix = NULL,
                        reducedDim = NULL,
                        ArchR_path = NULL,
                        reducedDimName = "LSI",
                        useDim = "IterativeLSI",
                        useMatrix = "GeneIntegrationMatrix",
                        cor_cutoff = 0.5,
                        cellNum = 200,
                        seed = 1,
                        ...) {
  set.seed(seed)

  if (!is.null(ArchR_path)) {
    ArchR::addArchRLogging(useLogs = FALSE)

    suppressMessages(obj <- ArchR::loadArchRProject(ArchR_path))

    obj <- ArchR::addPeak2GeneLinks(
      ArchRProj = obj,
      reducedDims = useDim,
      useMatrix = useMatrix,
      logFile = "x",
      ...
    )

    p2g <- ArchR::getPeak2GeneLinks(
      ArchRProj = obj,
      corCutOff = cor_cutoff,
      resolution = 1000,
      returnLoops = FALSE
    )

  } else if (!is.null(peakMatrix) &
             !is.null(expMatrix) & !is.null(reducedDim)) {
    # retrieve peak matrix
    #peakMatrix = obj[["PeakMatrix"]]

    # retrieve expression matrix
    #expMatrix = obj[[useMatrix]]
    #rownames(expMatrix) <- rowData(expMatrix)$name

    # retrieve dimensionality reduction
    #reducedDim <- SingleCellExperiment::reducedDims(obj[['TileMatrix500']])[[reducedDims]]

    # create sce object from expression matrix
    sce <-
      SingleCellExperiment::SingleCellExperiment(list(counts = assay(expMatrix)), altExps = list(peakMatrix =
                                                                                                   peakMatrix))
    # add reduced dimension information to sce object
    SingleCellExperiment::reducedDim(sce, reducedDimName) <-
      reducedDim

    # K-means clustering
    kNum <- trunc(ncol(sce) / cellNum)
    kclusters <-
      scran::clusterCells(
        sce,
        use.dimred = reducedDimName,
        BLUSPARAM = bluster::KmeansParam(centers = kNum, iter.max = 5000)
      )
    kclusters <- as.character(kclusters)

    # aggregate matrix by k-means clusters
    sce_grouped <-
      SingleCellExperiment::applySCE(sce,
                                     scuttle::aggregateAcrossCells,
                                     ids = kclusters,
                                     statistics = "sum")
    expGroupMatrix <- SummarizedExperiment::assay(sce_grouped)
    peakGroupMatrix <-
      SummarizedExperiment::assay(SingleCellExperiment::altExp(sce_grouped))

    # normalize group matrix
    expGroupMatrix <-
      t(t(expGroupMatrix) / colSums(expGroupMatrix)) * 10 ^ 4
    peakGroupMatrix <-
      t(t(peakGroupMatrix) / colSums(peakGroupMatrix)) * 10 ^ 4

    # get gene information
    geneSet <- rowRanges(expMatrix)
    geneStart <- GenomicRanges::resize(geneSet, 1, "start")

    # get peak range information
    peakSet <- rowRanges(peakMatrix)

    # find overlap after resizeing
    o <- DataFrame(
      GenomicRanges::findOverlaps(
        GenomicRanges::resize(geneStart, 2 * 100000 + 1, "center"),
        GenomicRanges::resize(peakSet, 1, "center"),
        ignore.strand = TRUE
      )
    )

    #Get Distance from Fixed point A B
    o$distance <-
      GenomicRanges::distance(geneStart[o[, 1]] , peakSet[o[, 2]])
    colnames(o) <- c("RNA", "ATAC", "distance")

    # Calculate correlation
    expCorMatrix <- expGroupMatrix[as.integer(o$RNA), ]
    peakCorMatrix <- peakGroupMatrix[as.integer(o$ATAC), ]
    o$Correlation <-
      mapply(cor, as.data.frame(t(expCorMatrix)), as.data.frame(t(peakCorMatrix)))
    o$VarATAC <- matrixStats::rowVars(peakCorMatrix)
    o$VarRNA <- matrixStats::rowVars(expCorMatrix)
    o$TStat <-
      (o$Correlation / sqrt((
        pmax(1 - o$Correlation ^ 2, 0.00000000000000001, na.rm = TRUE)
      ) / (ncol(peakCorMatrix) - 2))) #T-statistic P-value
    o$Pval <- 2 * pt(-abs(o$TStat), ncol(peakCorMatrix) - 2)
    o$FDR <- p.adjust(o$Pval, method = "fdr")
    p2g <-
      o[, c("RNA", "ATAC", "Correlation", "FDR", "VarATAC", "VarRNA")]
    colnames(p2g) <-
      c("idxRNA",
        "idxATAC",
        "Correlation",
        "FDR",
        "VarATAC",
        "VarRNA")
    metadata(p2g)$peakSet <- peakSet
    metadata(p2g)$geneSet <- geneStart

  } else {
    stop(
      "Input obj must be either 'ArchR' or all of the 3 matrices: gene expression, chromatin accessibility and dimensionality reduction"
    )

  }

  # Get metadata from p2g object and turn into df with peak indexes
  peak_metadata = as.data.frame(S4Vectors::metadata(p2g)[[1]]) # shows  chromosome, start, and end coordinates for each peak
  peak_metadata$idxATAC <- seq_along(rownames(peak_metadata))

  gene_metadata = as.data.frame(S4Vectors::metadata(p2g)[[2]]) # shows gene name and RNA index of genomic ranges
  gene_metadata$idxRNA = seq_along(rownames(gene_metadata))

  # Add gene names, chromosome num, chrom start, and chrom end to dataframe
  p2g <- as.data.frame(p2g)
  p2g_merged <-
    merge(p2g, gene_metadata, by = "idxRNA") # merge by gene ID
  p2g_merged <-
    merge(p2g_merged, peak_metadata, by = "idxATAC") # merge by peak ID

  p2g_merged <-
    p2g_merged[, c("idxATAC", "seqnames.x", "idxRNA", "name", "Correlation")]
  colnames(p2g_merged) <-
    c("idxATAC", "Chrom", "idxRNA", "Gene", "Correlation")
  p2g_merged <-
    p2g_merged[order(p2g_merged$idxATAC, p2g_merged$idxRNA), ]
  p2g_merged <- subset(p2g_merged, Correlation > 0.5)

  return(p2g_merged)
}
