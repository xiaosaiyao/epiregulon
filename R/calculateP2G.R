#' Establish peak to gene links based on correlations between ATAC-seq peaks and RNA-seq genes
#'
#' @param peakMatrix A SingleCellExperiment object containing counts of chromatin accessibility at each peak region or genomic bin from scATAC-seq
#' @param expMatrix A SingleCellExperiment object containing gene expression counts from scRNA-seq
#' @param reducedDim A matrix of dimension reduced values, for example derived from IterativeLSI algorithm of ArchR
#' @param ArchR_path String specifying the path to an ArchR project if ArchR's implementation of addPeak2GeneLinks is desired
#' @param cor_cutoff A numeric scalar to specify the correlation cutoff between ATAC-seq peaks and RNA-seq genes to assign peak to gene links.
#'  Default correlation cutoff is 0.5.
#' @param useDim String specifying the dimensional reduction representation in the ArchR project to use or the name of the reduced dimension matrix supplied by the user
#' @param useMatrix String specifying which the name of the gene expression matrix in the ArchR project to use.
#' It is often called the "GeneExpressionMatrix" for multiome and "GeneIntegrationMatrix" for unpaired data in ArchR project.
#' @param cellNum An integer to specify the number of cells to include in each K-means cluster. Default is 200 cells.
#' @param seed An integer scalar to specify the seed for K-means clustering
#' @param maxDist An integer to specify the base pair extension from transcription start start for overlap with peak regions
#' @param ... other parameters to pass to addPeak2GeneLinks from ArchR package
#'
#' @return A Peak2Gene correlation dataframe
#' @import SummarizedExperiment stats SingleCellExperiment GenomicRanges
#' @export
#'
#' @examples
#' # create a mock singleCellExperiment object for gene expression matrix
#' set.seed(1000)
#' gene_sce <- scuttle::mockSCE()
#' gene_sce <- scuttle::logNormCounts(gene_sce)
#' gene_gr <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr3","chr4"), nrow(gene_sce)/4),
#'                    ranges = IRanges(start = seq(from = 1, length.out=nrow(gene_sce), by = 1000),
#'                    width = 100))
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
#' @author Xiaosai Yao, Shang-yang Chen

calculateP2G <- function(peakMatrix = NULL,
                        expMatrix = NULL,
                        reducedDim = NULL,
                        ArchR_path = NULL,
                        useDim = "IterativeLSI",
                        useMatrix = "GeneIntegrationMatrix",
                        maxDist = 250000,
                        cor_cutoff = 0.5,
                        cellNum = 200,
                        seed = 1,
                        ...) {
  set.seed(seed)

  if (!is.null(ArchR_path)) {
    ArchR::addArchRLogging(useLogs = FALSE)

    writeLines("Using ArchR to compute peak to gene links...")
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

    # Get metadata from p2g object and turn into df with peak indexes
    peak_metadata <- as.data.frame(S4Vectors::metadata(p2g)[[1]]) # shows  chromosome, start, and end coordinates for each peak
    peak_metadata$idxATAC <- seq_along(rownames(peak_metadata))

    gene_metadata <- as.data.frame(S4Vectors::metadata(p2g)[[2]]) # shows gene name and RNA index of genomic ranges
    gene_metadata$idxRNA <- seq_along(rownames(gene_metadata))

    # Add gene names and peak positions to dataframe
    p2g <- as.data.frame(p2g)
    p2g_merged <- merge(p2g, gene_metadata, by = "idxRNA") # merge by gene ID
    p2g_merged <- merge(p2g_merged, peak_metadata, by = "idxATAC") # merge by peak ID

    # Extract relevant columns
    p2g_merged <- p2g_merged[, c("idxATAC", "seqnames.x", "start.x","start.y", "idxRNA", "name", "Correlation","FDR")]
    colnames(p2g_merged) <- c("idxATAC", "chr", "start","end", "idxRNA", "target", "Correlation","FDR")




  } else if (!is.null(peakMatrix) &
             !is.null(expMatrix) & !is.null(reducedDim)) {

    writeLines("Using epiregulon to compute peak to gene links...")

    # Package expression matrix and peak matrix into a single sce
    sce <- SingleCellExperiment(list(counts = assay(expMatrix)),
                                altExps = list(peakMatrix = peakMatrix))

    # add reduced dimension information to sce object
    reducedDim(sce,useDim) <- reducedDim


    writeLines("performing k means clustering to form metacells")
    # K-means clustering
    kNum <- trunc(ncol(sce) / cellNum)
    kclusters <- scran::clusterCells(sce,
                                     use.dimred =useDim,
                                     BLUSPARAM = bluster::KmeansParam
                                     (centers = kNum, iter.max = 5000))
    kclusters <- as.character(kclusters)

    # aggregate sce by k-means clusters
    sce_grouped <- applySCE(sce,
                            scuttle::aggregateAcrossCells,
                            ids = kclusters,
                            statistics = "mean")

    # some sces has strand information in metadata that conflicts with genomic ranges
    mcols(expMatrix)$strand <- NULL

    # transfer rowRanges(expMatrix) to rowranges(sce_grouped)
    rowRanges(sce_grouped) <- rowRanges(expMatrix)

    # rowRanges(altExp(sce_grouped)) already preserved

    # keep track of original ATAC and expression indices
    rowData(sce_grouped)$old.idxRNA <- 1:nrow(sce_grouped)
    rowData(altExp(sce_grouped))$old.idxATAC <- 1:nrow(altExp(sce_grouped))

    # remove genes and peaks that are equal to 0
    sce_grouped <- sce_grouped[which(rowSums(assay(sce_grouped)) != 0),]
    altExp(sce_grouped) <- altExp(sce_grouped)[which(rowSums(assay(altExp(sce_grouped))) != 0),]

    # extract gene expression and peak matrix
    expGroupMatrix <- assay(sce_grouped)
    peakGroupMatrix <- assay(altExp(sce_grouped))


    # get gene information
    geneSet <- rowRanges(sce_grouped)
    geneStart <- promoters(geneSet)

    # get peak range information
    peakSet <- rowRanges(altExp(sce_grouped))

    # find overlap after resizing
    o <- DataFrame(findOverlaps(resize(geneStart, maxDist, "center"),
                                peakSet,
                                ignore.strand = TRUE))


    #Get Distance from Fixed point A B
    o$distance <- distance(geneStart[o[, 1]] , peakSet[o[, 2]])
    colnames(o) <- c("RNA", "ATAC", "distance")


    # add old idxRNA and idxATAC
    o$old.idxRNA <- rowData(sce_grouped)[o[,1],"old.idxRNA"]
    o$old.idxATAC <- rowData(altExp(sce_grouped))[o[,2],"old.idxATAC"]

    # Calculate correlation
    expCorMatrix <- expGroupMatrix[as.integer(o$RNA), ]
    peakCorMatrix <- peakGroupMatrix[as.integer(o$ATAC), ]

    writeLines("Computing correlation")
    o$Correlation <- mapply(stats::cor,
                            as.data.frame(t(expCorMatrix)),
                            as.data.frame(t(peakCorMatrix)))

    o$VarATAC <- matrixStats::rowVars(peakCorMatrix)
    o$VarRNA <- matrixStats::rowVars(expCorMatrix)
    o$TStat <- (o$Correlation /
                  sqrt((pmax(1 - o$Correlation ^ 2, 0.00000000000000001, na.rm = TRUE))
                       / (ncol(peakCorMatrix) - 2))) #T-statistic P-value
    o$Pval <- 2 * pt(-abs(o$TStat), ncol(peakCorMatrix) - 2)
    o$FDR <- p.adjust(o$Pval, method = "fdr")

    #add metadata to o
    o$Gene <-  rowData(sce_grouped)[o[,1],"name"]
    o$chr <- as.character(seqnames(rowRanges(altExp(sce_grouped))[o[,2]]))
    o$start <- GenomicRanges::start(rowRanges(altExp(sce_grouped))[o[,2],])
    o$end <- GenomicRanges::end(rowRanges(altExp(sce_grouped))[o[,2],])

    o <- data.frame(o)

    p2g_merged <- o[, c("old.idxATAC", "chr","start","end", "old.idxRNA", "Gene", "Correlation", "FDR")]
    colnames(p2g_merged) <- c("idxATAC", "chr", "start","end", "idxRNA", "target", "Correlation","FDR")
    p2g_merged <- p2g_merged[p2g_merged$Correlation > cor_cutoff,,drop=FALSE]



  } else {
    stop(
      "Input obj must be either 'ArchR' or all of the 3 matrices: gene expression, chromatin accessibility and dimensionality reduction"
    )

  }
  p2g_merged <- p2g_merged[order(p2g_merged$idxATAC, p2g_merged$idxRNA),,drop=FALSE]
  return(p2g_merged)

}
