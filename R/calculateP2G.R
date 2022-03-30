#' A function to compute correlations between ATAC-seq peaks and RNA-seq genes
#'
#' @param peakMatrix SingleCellExperiment object containing counts of scATAC-seq peaks
#' @param expMatrix SingleCellExperiment object containing counts of scRNA-seq expressions
#' @param reducedDim matrix of dimensional reduced values, can be derived from IterativeLSI algorithm of ArchR
#' @param ArchR_path string specifying optional path to a ArchR project if ArchR's implementation of addPeak2GeneLinks is desired
#' @param cor_cutoff cutoff for correlations between ATAC-seq peaks and RNA-seq genes
#' @param cellNum number of cells to include in each K-means cluster, to control number of k
#' @param seed number of random seed to use for K-means clustering
#' @param ... other parameters to pass to addPeak2GeneLinks from ArchR package
#'
#' @return A Peak2Gene correlation datafrane
#' @export
#'
#' @examples
#' p2g <- calculateP2G(peakmatrix, expmatrix, reducedDim)
#' head(p2g)
#'
calculateP2G = function(peakMatrix = NULL, expMatrix = NULL, reducedDim = NULL, ArchR_path = NULL, cor_cutoff = 0.5, cellNum = 200, seed = 1, ...){

  set.seed(seed)

  if (!is.null(ArchR_path)){

    ArchR::addArchRLogging(useLogs = FALSE)

    suppressMessages(obj <- ArchR::loadArchRProject(Archr_path))

    obj <- ArchR::addPeak2GeneLinks(
      ArchRProj = obj,
      reducedDims = reducedDims,
      useMatrix = useMatrix,
      logFile = "x", ...
    )

    p2g <- ArchR::getPeak2GeneLinks(
      ArchRProj = obj,
      corCutOff = cor_cutoff,
      resolution = 1000,
      returnLoops = FALSE
    )

  } else if (!is.null(peakMatrix) & !is.null(expMatrix) & !is.null(reducedDim)) {

    # retrieve peak matrix
    #peakMatrix = obj[["PeakMatrix"]]

    # retrieve expression matrix
    #expMatrix = obj[[useMatrix]]
    #rownames(expMatrix) <- rowData(expMatrix)$name

    # retrieve dimensionality reduction
    #reducedDim <- SingleCellExperiment::reducedDims(obj[['TileMatrix500']])[[reducedDims]]

    # create sce object from expression matrix
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts=assay(expMatrix)), altExps = list(peakMatrix=peakMatrix))
    # add reduced dimension information to sce object
    SingleCellExperiment::reducedDim(sce, "LSI") <- reducedDim

    # K-means clustering
    kNum <- trunc(ncol(sce)/cellNum)
    kclusters <- scran::clusterCells(sce, use.dimred="LSI", BLUSPARAM=bluster::KmeansParam(centers=kNum, iter.max = 5000))
    kclusters <- as.character(kclusters)

    # aggregate matrix by k-means clusters
    sce_grouped <- SingleCellExperiment::applySCE(sce, scuttle::aggregateAcrossCells, ids=kclusters, statistics = "sum")
    expGroupMatrix <- SummarizedExperiment::assay(sce_grouped)
    peakGroupMatrix <- SummarizedExperiment::assay(SingleCellExperiment::altExp(sce_grouped))

    # normalize group matrix
    expGroupMatrix <- t(t(expGroupMatrix) / colSums(expGroupMatrix)) * 10^4
    peakGroupMatrix <- t(t(peakGroupMatrix) / colSums(peakGroupMatrix)) * 10^4

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
      ))

    #Get Distance from Fixed point A B
    o$distance <- GenomicRanges::distance(geneStart[o[,1]] , peakSet[o[,2]])
    colnames(o) <- c("RNA", "ATAC", "distance")

    # Calculate correlation
    expCorMatrix <- expGroupMatrix[as.integer(o$RNA),]
    peakCorMatrix <- peakGroupMatrix[as.integer(o$ATAC),]
    o$Correlation <- mapply(cor, as.data.frame(t(expCorMatrix)), as.data.frame(t(peakCorMatrix)))
    o$VarATAC <- matrixStats::rowVars(peakCorMatrix)
    o$VarRNA <- matrixStats::rowVars(expCorMatrix)
    o$TStat <- (o$Correlation / sqrt((pmax(1-o$Correlation^2, 0.00000000000000001, na.rm = TRUE))/(ncol(peakCorMatrix)-2))) #T-statistic P-value
    o$Pval <- 2*pt(-abs(o$TStat), ncol(peakCorMatrix) - 2)
    o$FDR <- p.adjust(o$Pval, method = "fdr")
    p2g <- o[, c("RNA", "ATAC","Correlation", "FDR", "VarATAC", "VarRNA")]
    colnames(p2g) <- c("idxRNA", "idxATAC", "Correlation", "FDR", "VarATAC", "VarRNA")
    metadata(p2g)$peakSet <- peakSet
    metadata(p2g)$geneSet <- geneStart

  } else {

    stop("Input obj must be either 'ArchR' or individual assays from 'MultiAssayExperiment'")

  }

  # Get metadata from p2g object and turn into df with peak indexes
  peak_metadata = as.data.frame(S4Vectors::metadata(p2g)[[1]]) # shows  chromosome, start, and end coordinates for each peak
  peak_metadata$idxATAC <- seq_along(rownames(peak_metadata))

  gene_metadata = as.data.frame(S4Vectors::metadata(p2g)[[2]]) # shows gene name and RNA index of genomic ranges
  gene_metadata$idxRNA = seq_along(rownames(gene_metadata))

  # Add gene names, chromosome num, chrom start, and chrom end to dataframe
  p2g <- as.data.frame(p2g)
  p2g_merged <- merge(p2g, gene_metadata, by = "idxRNA") # merge by gene ID
  p2g_merged <- merge(p2g_merged, peak_metadata, by = "idxATAC") # merge by peak ID

  p2g_merged <- p2g_merged[, c("idxATAC","seqnames.x","idxRNA","name","Correlation")]
  colnames(p2g_merged) <- c("idxATAC","Chrom","idxRNA", "Gene","Correlation")
  p2g_merged <- p2g_merged[order(p2g_merged$idxATAC,p2g_merged$idxRNA),]
  p2g_merged <- subset(p2g_merged, Correlation > 0.5)

  return(p2g_merged)
}
