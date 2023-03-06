#' Compute correlations between ATAC-seq peaks and RNA-seq genes
#'
#' @param archr_path String indicating the path to a ArchR project that has already performed LSI dimensionality reduction and scRNA-seq integration
#' @param cor_cutoff A numeric scalar to indicate the cutoff for correlations between ATAC-seq peaks and RNA-seq genes
#' @param reducedDims String specifying which dimensional reduction representation in the ArchR project to use
#' @param useMatrix String specifying which data matrix in the ArchR project to use
#' @param ... other parameters to pass to addPeak2GeneLinks from ArchR package
#'
#' @return A Peak2Gene correlation data frame
#' @import ArchR utils
#' @export
#' @author Shang-yang Chen


getP2Glinks <- function(archr_path,
                        cor_cutoff = 0.5,
                        reducedDims = "IterativeLSI",
                        useMatrix = "GeneIntegrationMatrix",
                        ...){

 .Deprecated("calculateP2G")

 ArchR::addArchRLogging(useLogs = FALSE)

  suppressMessages(proj <- ArchR::loadArchRProject(archr_path))

  proj <- ArchR::addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = reducedDims,
    useMatrix = useMatrix,
    logFile = "x",
    ...
  )

  p2g <- ArchR::getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = cor_cutoff,
    resolution = 1000,
    returnLoops = FALSE
  )

  # Get metadata from p2g object and turn into df with peak indexes
  peak_metadata <- as.data.frame(S4Vectors::metadata(p2g)[[1]]) # shows  chromosome, start, and end coordinates for each peak
  peak_metadata$idxATAC <- as.numeric(rownames(peak_metadata))

  gene_metadata <- as.data.frame(S4Vectors::metadata(p2g)[[2]]) # shows gene name and RNA index of genomic ranges
  gene_metadata$idxRNA <- as.numeric(rownames(gene_metadata))

  # Add gene names, chromosome num, chrom start, and chrom end to dataframe
  p2g <- as.data.frame(p2g)
  p2g_merged <- merge(p2g, gene_metadata, by = "idxRNA") # merge by gene ID
  p2g_merged <- merge(p2g_merged, peak_metadata, by = "idxATAC") # merge by peak ID

  p2g_merged <- p2g_merged[, c("idxATAC","seqnames.x","idxRNA","name","Correlation")]
  colnames(p2g_merged) <- c("idxATAC","Chrom","idxRNA", "Gene","Correlation")
  p2g_merged <- p2g_merged[order(p2g_merged$idxATAC,p2g_merged$idxRNA),]

  return(p2g_merged)
}


#' An accessor function to retrieve TF motif info from [scMultiome](https://github.com/xiaosaiyao/scMultiome/blob/master/R/tfBinding.R)
#'
#' Combined transcription factor ChIP-seq data from ChIP-Atlas and ENCODE or
#' from CistromeDB and ENCODE.
#' @inherit scMultiome::tfBinding params return references
#' @examples
#' \dontrun{
#' getTFMotifInfo("mm10", "atlas")
#' }
#'
#' @export
#'
getTFMotifInfo <- scMultiome::tfBinding

#' Add TF binding motif occupancy information to the peak2gene object
#'
#' @param p2g A Peak2Gene dataframe created by ArchR or getP2Glinks() function
#' @param grl GRangeList object containing reference TF binding information. We recommend retrieving `grl` from `getTFMotifInfo` which
#' contains TF occupancy data derived from public and ENCODE ChIP-seq peaks. Alternatively, if the users would like to provide a GRangeList
#' of motif annotations. This can be derived using `motifmatchr::matchMotifs`. See details
#' @param peakMatrix A matrix of scATAC-seq peak regions with peak ids as rows
#' @param archR_project_path Path to an ArchR project that have performed LSI dimensionality reduction and scRNA-seq integration
#'
#' @return A data frame containing overlapping ids of scATAC-seq peak regions and reference TF binding regions
#' @details This function annotates each regulatory element with possible transcription factors. We can either provide a GRangeList of known ChIP-seq
#' binding sites (TF occupancy) or positions of TF motifs (TF motifs). While public ChIP-seq data may not fully align with the ground truth TF occupancy in users' data
#' (due to technical challenges of ChIP-seq or cell type nature of TF occupancy), it does offer a few important advantages over TF motif information:
#' \enumerate{
#' \item{TF occupancy allows co-activators to be included. Co-activators are chromatin modifiers that do not directly bind to DNA but nonetheless play an important role
#' in gene regulation}
#' \item{TF occupancy can distinguish between members of the same class that may share similar motifs but that may have drastically different binding sites}
#' }
#' If multiple ChIP-seq are available for the same TF, we merge the ChIP-seq data to represent an universal set of possible binding sites. The predicted TF
#' occupancy is further refined by \code{\link{pruneRegulon}}.
#'
#' If the users prefer to use TF motifs instead of TF occupancy, the users can create a GRangeList of motif annotation using `motifmatchr::matchMotifs`.
#' Here, we demonstrate how to annotate peaks with cisbp motif database
#' ```
#' library(motifmatchr)
#' library(chromVARmotifs)
#' data("human_pwms_v1")
#' peaks <- GenomicRanges::GRanges(seqnames = c("chr1","chr2","chr2"),
#'                                ranges = IRanges::IRanges(start = c(76585873,42772928, 100183786),
#'                                                          width = 500))
#' grl <- matchMotifs(human_pwms_v1, peaks, genome = "hg38", out = "positions")
#' # retain only TF symbols. TF symbols need to be consistent with gene names in regulon
#' names(grl) <- sapply(strsplit(names(grl), "_"), "[",3)
#' ```
#'
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' # create a mock peak-to-gene matrix
#' p2g <- data.frame(idxATAC = c(rep(1,5), rep(2,5)), Chrom = "chr1", idxRNA = 1:10,
#' Gene = paste0("Gene_",1:10),Correlation = runif(10, 0,1))
#'
#' # create mock a GRanges list of TF binding sites
#' grl <- GRangesList("TF1" = GRanges(seqnames = "chr1",
#' ranges = IRanges(start = c(50,1050), width = 100)),
#' "TF2" = GRanges(seqnames = "chr1",
#' ranges = IRanges(start = c(1050), width = 100))
#' )
#'
#' # create a mock singleCellExperiment object for peak matrix
#' peak_gr <- GRanges(seqnames = "chr1",
#'              ranges = IRanges(start = seq(from = 1, to = 10000, by = 1000),
#'              width = 100))
#' peak_counts <- matrix(sample(x = 0:4, size = 100*length(peak_gr), replace = TRUE),
#' nrow = length(peak_gr), ncol = 100)
#' peak_sce <- SingleCellExperiment(list(counts = peak_counts))
#' rowRanges(peak_sce) <- peak_gr
#' rownames(peak_sce) <- paste0("peak",1:10)
#'
#' # create overlaps between p2g matrix, TF binding sites and peak matrix
#' overlap <- addTFMotifInfo(p2g, grl, peakMatrix = peak_sce)
#' utils::head(overlap)
#' @author Xiaosai Yao, Shang-yang Chen

addTFMotifInfo <- function(p2g,
                           grl,
                           peakMatrix = NULL,
                           archR_project_path = NULL){

  if (!is.null(archR_project_path)) {
    proj <- loadArchRProject(path = archR_project_path, showLogo = FALSE)
    peakSet <- getPeakSet(ArchRProj = proj)
  } else {
    peakSet <- rowRanges(peakMatrix)
  }

  message("Computing overlap...")
  overlap <- GenomicRanges::findOverlaps(peakSet, grl)
  overlap <- data.frame(overlap)
  colnames(overlap) <- c("idxATAC", "idxTF")
  overlap <- overlap[which(overlap$idxATAC %in% p2g$idxATAC), ]
  overlap$tf <- names(grl)[overlap$idxTF]
  message("Success!")

  return(overlap)

}


#' Combine the TF binding info and peak to gene correlations to generate regulons
#'
#' @param p2g A Peak2Gene data frame created by ArchR or getP2Glinks() function
#' @param overlap A data frame storing overlaps between the regions of the peak matrix with the bulk TF ChIP-seq binding sites computed from addTFMotifInfo
#' @param aggregate logical to specify whether peak and gene ids are kept in regulon output or not
#' @param FUN function to aggregate the weights
#'
#' @return A DataFrame consisting of tf(regulator), target and a column indicating degree of association between TF and target such as "mor" or "corr".
#'
#' @export
#'
#'
#' @examples
#' set.seed(1)
#' # create a mock peak-to-gene matrix
#' p2g <- data.frame(idxATAC = c(rep(1,5), rep(2,5)), Chrom = "chr1", idxRNA = 1:10,
#' target = paste0("Gene_", 1:10), Correlation = runif(10, 0, 1))
#'
#' # create a Granges list of TF binding sites
#' grl <- GRangesList("TF1" = GRanges(seqnames = "chr1",
#' ranges = IRanges(start = c(50,1050), width = 100)),
#' "TF2" = GRanges(seqnames = "chr1",
#' ranges = IRanges(start = c(1050), width = 100))
#' )
#'
#' # Create a mock peak matrix
#' peak_gr <- GRanges(seqnames = "chr1",
#'                    ranges = IRanges(start = seq(from = 1, to = 10000, by = 1000), width = 100))
#'
#' peak_counts <- matrix(sample(x = 0:4, size = 100*length(peak_gr), replace = TRUE),
#' nrow = length(peak_gr),ncol = 100)
#' peak_sce <- SingleCellExperiment(list(counts = peak_counts))
#' rowRanges(peak_sce) <- peak_gr
#' rownames(peak_sce) <- paste0("peak", 1:10)
#'
#' # create overlaps between p2g matrix, TF binding sites and peak matrix
#' overlap <- addTFMotifInfo(p2g, grl, peakMatrix = peak_sce)
#' utils::head(overlap)
#'
#' # aggregate gene expression if the gene is bound by the same TF at regulatory elements
#' regulon <- getRegulon(p2g, overlap, aggregate = FALSE)
#' @author Xiaosai Yao, Shang-yang Chen

getRegulon <- function(p2g,
                       overlap,
                       aggregate=FALSE,
                       FUN=colMeans){

  p2g <- S4Vectors::DataFrame(p2g)
  regulon_df <- S4Vectors::merge(p2g, overlap, by="idxATAC")

  Correlation.rownames <- colnames(regulon_df)[grep("Correlation.|Correlation", colnames(regulon_df))]
  corr_matrix <- regulon_df[,Correlation.rownames, drop=FALSE]

  if (any(grepl("Correlation.", Correlation.rownames))){
    colnames(corr_matrix) <- gsub("Correlation.","", Correlation.rownames)
  }

  regulon_df[,grep("Correlation.",colnames(regulon_df))] <- NULL
  regulon_df$Correlation <- as.matrix(corr_matrix)


  if (aggregate) {
    "aggregating regulon ..."
    regulon_df <- aggregateMatrix.DF(regulon_df[,c("tf","target","Correlation")], "Correlation", colMeans)
  }
  colnames(regulon_df)[colnames(regulon_df) == "Correlation"] <- "corr"
  return(regulon_df)

}


aggregateMatrix.DF <- function(regulon, mode, FUN){
  regulon$tf <- as.factor(regulon$tf)
  regulon$target <- as.factor(regulon$target)
  groupings <- interaction(regulon$tf,regulon$target, sep = '_')
  index <- order(groupings)
  regulon <- regulon[index,]
  breaks <- which(!duplicated(groupings[index]))
  aggregated <- lapply(seq_len(length(breaks)-1), function(i){
    FUN(as.matrix(regulon[breaks[i]:(breaks[i+1]-1), mode, drop=FALSE]))})

  aggregated[[length(breaks)]] <-
    FUN(as.matrix(regulon[breaks[length(breaks)]:nrow(regulon), mode, drop=FALSE]))
  aggregated <- do.call(rbind, aggregated)
  aggregated <- S4Vectors::DataFrame(tf=regulon$tf[breaks],
                                     target=regulon$target[breaks],
                                     Correlation=I(aggregated))
}
