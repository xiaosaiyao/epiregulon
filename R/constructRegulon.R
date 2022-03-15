#' A function to compute correlations between ATAC-seq peaks and RNA-seq genes
#'
#' @param archr_path Path to a ArchR project that have performed LSI dimensionality reduction and scRNA-seq integration
#' @param cor_cutoff Cutoff for correlations between ATAC-seq peaks and RNA-seq genes
#' @param reducedDims String specifying which dimensional reduction representation in the ArchR project to use
#' @param useMatrix String specifying which data matrix in the ArchR project to use
#' @param ... other parameters to pass to addPeak2GeneLinks from ArchR package
#'
#' @return A Peak2Gene correlation datafrane
#' @import ArchR
#' @import utils
#' @export
#'
#' @examples
#' p2g <- getP2Glinks("/gstore/project/lineage/sam/heme_GRN/OUTPUT")
#' head(p2g)
#'
getP2Glinks <- function(archr_path, cor_cutoff = 0.5, reducedDims = "IterativeLSI", useMatrix = "GeneIntegrationMatrix", ...){

 ArchR::addArchRLogging(useLogs = FALSE)

  suppressMessages(proj <- ArchR::loadArchRProject(archr_path))

  proj <- ArchR::addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = reducedDims,
    useMatrix = useMatrix,
    logFile = "x", ...
  )

  p2g <- ArchR::getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = cor_cutoff,
    resolution = 1000,
    returnLoops = FALSE
  )

  # Get metadata from p2g object and turn into df with peak indexes
  peak_metadata = as.data.frame(S4Vectors::metadata(p2g)[[1]]) # shows  chromosome, start, and end coordinates for each peak
  peak_metadata$idxATAC <- as.numeric(rownames(peak_metadata))

  gene_metadata = as.data.frame(S4Vectors::metadata(p2g)[[2]]) # shows gene name and RNA index of genomic ranges
  gene_metadata$idxRNA = as.numeric(rownames(gene_metadata))

  # Add gene names, chromosome num, chrom start, and chrom end to dataframe
  p2g <- as.data.frame(p2g)
  p2g_merged <- merge(p2g, gene_metadata, by = "idxRNA") # merge by gene ID
  p2g_merged <- merge(p2g_merged, peak_metadata, by = "idxATAC") # merge by peak ID

  p2g_merged <- p2g_merged[, c("idxATAC","seqnames.x","idxRNA","name","Correlation")]
  colnames(p2g_merged) <- c("idxATAC","Chrom","idxRNA", "Gene","Correlation")
  p2g_merged <- p2g_merged[order(p2g_merged$idxATAC,p2g_merged$idxRNA),]

  return(p2g_merged)
}


#' An accessor function to retrieve TF motif info from genomitory repository
#'
#' @param genome a string specifying the genome for bed files
#'
#' @return A GRangeList object containing binding site information of 1274 TFs
#' @export
#'
#' @examples
#' grl <- getTFMotifInfo(genome="mm10")
#' head(grl)
#'
getTFMotifInfo <- function(genome = "hg19"){

  if (genome == 'hg19'){

    id <- genomitory::packID("GMTY156", # project name
                             "hg19_motif_bed_granges.rds", # path within the project
                             "REVISION-2") # version
    grl <- genomitory::getFeatures(id)

  } else if (genome == "hg38") {

    id <- genomitory::packID("GMTY162", # project name
                             "hg38_motif_bed_granges.rds", # path within the project
                             "REVISION-1") # version
    grl <- genomitory::getFeatures(id)

  } else if (genome == "mm10"){
    id <- genomitory::packID("GMTY181", # project name
                             "mm10_motif_bed_granges.rds", # path within the project
                             "REVISION-2") # version
    grl <- genomitory::getFeatures(id)
  }

  return(grl)
}

#' A function to add TF binding motif occupancy information to the peak2gene object
#'
#' @param p2g A Peak2Gene dataframe created by ArchR or getP2Glinks() function
#' @param grl GRangeList object containing reference TF binding information
#' @param peakMatrix A matrix of scATAC-seq peak regions with peak ids as rows
#' @param archR_project_path Path to a ArchR project that have performed LSI dimensionality reduction and scRNA-seq integration
#'
#' @return A dataframe containing overlapping ids of scATAC-seq peak regions and reference TF binding regions
#' @export
#'
#' @examples
#' p2g <- getP2Glinks("/gstore/project/lineage/sam/heme_GRN/OUTPUT")
#' grl <- getTFMotifInfo(genome="mm10")
#' overlap <- addTFMotifInfo(p2g, grl, archR_project_path = "/gstore/project/lineage/sam/heme_GRN/OUTPUT")
#' head(overlap)
addTFMotifInfo <- function(p2g, grl, peakMatrix=NULL, archR_project_path=NULL){

  if (!is.null(archR_project_path)) {
    proj <- loadArchRProject(path = archR_project_path, showLogo = F)
    peakSet= getPeakSet(ArchRProj = proj)
  } else {
    peakSet = rowRanges(peakMatrix)
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


#'  A function to combine the TF binding motif info and peak to gene correlations to generate regulons
#'
#' @param p2g A Peak2Gene dataframe created by ArchR or getP2Glinks() function
#' @param overlap dataframe storing overlaps between the regions of the peak matrix with the bulk TF ChIP-seq binding sites
#' @param aggregate boolean value that specify whether peak and gene ids are kept in regulon output or not
#'
#' @return a tall format dataframe consisting of tf(regulator), target and a column indicating degree of association between TF and target such as "mor" or "corr".
#'           example regulon:
#'           tf      target  corr
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' p2g <- getP2Glinks("/gstore/project/lineage/sam/heme_GRN/OUTPUT")
#' grl <- getTFMotifInfo(genome="mm10")
#' overlap <- addTFMotifInfo(p2g, grl, archR_project_path = "/gstore/project/lineage/sam/heme_GRN/OUTPUT")
#' regulon <- getRegulon(p2g, overlap, aggregate = T)
#' head(regulon)
getRegulon <- function(p2g, overlap, aggregate = TRUE){

  regulon_df <- merge(overlap, p2g, by="idxATAC")

  if (aggregate) {
    regulon_df <- regulon_df[, c("tf", "Gene", "Correlation")]
    colnames(regulon_df) <- c("tf", "target", "corr")
    regulon_df <- aggregate(corr ~ tf + target, data = regulon_df,
                            FUN = mean, na.rm = TRUE)
  } else {
    colnames(regulon_df) <- c("idxATAC","idxTF", "tf", "Chrom", "idxRNA", "target","corr")
  }
  return(regulon_df)

}
