#' A function to compute correlations between ATAC-seq peaks and RNA-seq genes
#'
#' @param archr_path Path to a ArchR project that have performed LSI dimensionality reduction and scRNA-seq integration
#' @param cor_cutoff Cutoff for correlations between ATAC-seq peaks and RNA-seq genes
#'
#' @return A Peak2Gene object (DFrame)
#' @import ArchR
#' @import utils
#' @export
#'
#' @examples 1+1
#'
getP2Glinks <- function(archr_path, cor_cutoff = 0.5, reducedDims = "IterativeLSI", useMatrix = "GeneIntegrationMatrix", ...){

 ArchR::addArchRLogging(useLogs = FALSE)

  proj <- ArchR::loadArchRProject(archr_path)

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

  return(p2g)
}


#' An accessor function to retrieve TF motif info from genomitory repository
#'
#' @param genome a string specifying the genome for bed files
#'
#' @return A GRangeList object containing binding site information of 1274 TFs
#' @export
#'
#' @examples 1+1
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

  }

  return(grl)
}

#' A function to add TF binding motif occupancy information to the peak2gene object
#'
#' @param archr_path Path to a ArchR project that have performed LSI dimensionality reduction and scRNA-seq integration
#' @param organism human or mouse are currently supported
#' @param grl GRangeList object containing TF binding information
#'
#' @return None, a motif binary matrix file will be created in the Annotations folder of supplied ArchR project path as ChIP-Matches-In-Peaks.rds
#' @export
#'
#' @examples 1+1
addTFMotifInfo <- function(archr_path, grl, organism = "human"){

  proj <- ArchR::loadArchRProject(archr_path)

  if (organism=="human"){
    proj <- ArchR::addPeakAnnotations(ArchRProj = proj, regions = grl, name = "ChIP", force = TRUE,
                                      logFile = "x")
  } else{
    proj <- ArchR::addPeakAnnotations(ArchRProj = proj, regions = grl, name = "ChIP", force = TRUE,
                                      logFile = "x")
  }

  #### check if binary matrix is created
  if (file.exists(paste0(archr_path,"/Annotations/ChIP-Matches-In-Peaks.rds"))) {

    return ("Add TF motif info sucessful!")

  } else {

    stop("TF binary matrix file was not created.")

  }
}


#'  A function to combine the TF binding motif info and peak to gene correlations to generate regulons
#'
#' @param archr_path Path to a ArchR project that have performed LSI dimensionality reduction and scRNA-seq integration
#' @param p2g_object A Peak2Gene object (DFrame) created by ArchR or getP2Glinks() function
#' @param keep_id boolean value that specify whether peak and gene ids are kept in regulon output or not
#'
#' @return a tall format dataframe consisting of tf(regulator), target and a column indicating degree of association between TF and target such as "mor" or "corr".
#'           example regulon:
#'           tf      target  corr
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples 1+1
getRegulon <- function(archr_path, p2g_object, keep_id = FALSE){

  # Get metadata from p2g object and turn into df with peak indexes
  peak_metadata = as.data.frame(S4Vectors::metadata(p2g_object)[[1]]) # shows  chromosome, start, and end coordinates for each peak
  peak_metadata$idxATAC <- as.numeric(rownames(peak_metadata))

  gene_metadata = as.data.frame(S4Vectors::metadata(p2g_object)[[2]]) # shows gene name and RNA index of genomic ranges
  gene_metadata$idxRNA = as.numeric(rownames(gene_metadata))

  # Add gene names, chromosome num, chrom start, and chrom end to dataframe
  p2g <- as.data.frame(p2g_object)
  p2g_merged <- merge(p2g, gene_metadata, by = "idxRNA") # merge by gene ID
  p2g_merged <- merge(p2g_merged, peak_metadata, by = "idxATAC") # merge by peak ID

  p2g_merged <- p2g_merged[, c("idxATAC","seqnames.x","idxRNA","name","Correlation")]
  colnames(p2g_merged) <- c("idxATAC","Chrom","idxRNA", "Gene","Correlation")
  p2g_merged <- p2g_merged[order(p2g_merged$idxATAC,p2g_merged$idxRNA),]

  ## Append TF occupancy information to create regulon
  tf.binary <- readRDS(paste0(archr_path,'/Annotations/ChIP-Matches-In-Peaks.rds'))
  tf.binary.matrix <- as.matrix(tf.binary@assays@data$matches)
  tf.binary.matrix <- as.data.frame(tf.binary.matrix)
  tf.binary.matrix$idxATAC <- as.numeric(rownames(tf.binary.matrix))
  regulon_wide <- merge(p2g_merged, tf.binary.matrix, by="idxATAC")

  ### convert into long, readable matrix format
  #regulon_df <- within(regulon_wide, rm("idxATAC","Chrom","idxRNA"))
  regulon_df <- tidyr::pivot_longer(regulon_wide, -c(.data$idxATAC, .data$idxRNA, .data$Chrom, .data$Gene, .data$Correlation), names_to = "TF")
  regulon_df <- regulon_df[regulon_df$value==TRUE, ]

  if (keep_id){
    regulon_df <- regulon_df[,c("idxATAC","idxRNA","TF","Gene","Correlation")]
    colnames(regulon_df) <- c("idxATAC","idxRNA", "tf", "target","corr")

  } else {
    ### aggregate multiple tf-target rows by mean
    regulon_df <- regulon_df[,c("TF","Gene","Correlation")]
    colnames(regulon_df) <- c("tf", "target","corr")
    regulon_df <- aggregate(corr ~ tf + target, data = regulon_df, FUN = mean, na.rm = TRUE)
  }

  return (regulon_df)
}
