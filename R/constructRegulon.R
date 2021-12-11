#' A function to compute correlations between ATAC-seq peaks and RNA-seq genes
#'
#' @param archr_path Path to a ArchR project that have performed LSI dimensionality reduction and scRNA-seq integration
#' @param cor_cutoff Cutoff for correlations between ATAC-seq peaks and RNA-seq genes
#'
#' @return A Peak2Gene object (DFrame)
#' @export
#'
#' @examples 1+1
#'
getP2Glinks <- function(archr_path, cor_cutoff = 0.5){

  proj <- ArchR::loadArchRProject(archr_path)

  proj <- ArchR::addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    logFile = createLogFile("addPeak2GeneLinks",logDir = paste0(archr_path,"/ArchRLogs"))
  )

  p2g <- ArchR::getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = cor_cutoff,
    resolution = 1000,
    returnLoops = FALSE
  )

  return(p2g)
}

addTFMotifInfo <- function(archr_path, organism = "human"){

  proj <- ArchR::loadArchRProject(archr_path)

  if (organism=="human"){
    proj <- ArchR::addPeakAnnotations(ArchRProj = proj, regions = human_bed, name = "ChIP", force = TRUE,
                                      logFile = createLogFile("addPeakAnnotations",logDir = paste0(archr_path,"/ArchRLogs")))
  } else{
    proj <- ArchR::addPeakAnnotations(ArchRProj = proj, regions = human_bed, name = "ChIP", force = TRUE,
                                      logFile = createLogFile("addPeakAnnotations",logDir = paste0(archr_path,"/ArchRLogs")))
  }

  #### check if binary matrix is created
  if (file.exists(paste0(archr_path,"/Annotations/ChIP-Matches-In-Peaks.rds"))) {

    return ("Add TF motif info sucessful!")

  } else {

    stop("TF binary matrix file was not created.")

  }
}

getRegulon <- function(archr_path, p2g_object){

  # Get metadata from p2g object and turn into df with peak indexes
  peak_metadata = as.data.frame(metadata(p2g)[[1]]) %>% # shows  chromosome, start, and end coordinates for each peak
    rownames_to_column("idxATAC") %>% mutate(idxATAC = as.numeric(idxATAC))

  gene_metadata = as.data.frame(metadata(p2g)[[2]]) %>% # shows gene name and RNA index of genomic ranges
    rownames_to_column("idxRNA") %>% mutate(idxRNA =  as.numeric(idxRNA))

  # Add gene names, chromosome num, chrom start, and chrom end to dataframe
  p2g_merged <- as.data.frame(p2g) %>% left_join(peak_metadata %>% select(idxATAC:end)) %>% # add peak info
    left_join(gene_metadata %>% select(idxRNA,name)) %>% # add gene info
    select(idxATAC,seqnames:end,idxRNA,name,Correlation) %>% arrange(idxATAC, idxRNA) %>%
    dplyr::rename(Chrom = seqnames, Gene = name)

  ## Append TF occupancy information to create regulon
  tf.binary <- readRDS(paste0(archr_path,'/Annotations/ChIP-Matches-In-Peaks.rds'))
  tf.binary.matrix <- as.matrix(tf.binary@assays@data$matches)
  tf.binary.matrix <- as.data.frame(tf.binary.matrix) %>% rownames_to_column(var="idxATAC") %>%
    mutate(idxATAC = as.numeric(idxATAC))
  regulon_wide <- p2g_merged %>% left_join(tf.binary.matrix)

  ### convert into long, readable matrix format
  regulon_df <- regulon_wide %>% select(Gene:last_col())  %>%
    pivot_longer(-c(Gene,Correlation), names_to = "TF") %>%
    filter(value==TRUE)  %>% select(TF, Gene, Correlation)

  return (regulon_df)
}



