#' Add Motif Scores
#'
#' @param regulon A DataFrame consisting of tf (regulator) and target in the column names.
#' @param field_name Character string	indicating the column name of the regulon to add the motif information to
#' @param peaks A GRanges object indicating the peaks to perform motif annotation on if ArchR project is not provided.
#' The peak indices should match the `re` column in the regulon
#' @param pwms A PWMatrixList for annotation of motifs using 'motifmatchr::matchMotifs'
#' @param species Character string indicating species. Currently supported species is human or mouse
#' @param genome Character string indicating the genomic build
#' @param ... Additional arguments to pass into motifmatchr::matchMotifs
#'
#' @return A DataFrame with motif matches added with 1s indicating the presence of motifs and
#' 0s indicating the absence of motifs
#' @export
#'
#' @examples
#' regulon <- S4Vectors::DataFrame(tf = c("AR","AR","AR","ESR1","ESR1","NKX2-1"),
#' idxATAC = 1:6)
#' peaks <- GRanges(seqnames = c("chr12","chr19","chr19","chr11","chr6","chr1"),
#' ranges = IRanges(start = c(124914563,50850845, 50850844, 101034172, 151616327, 1000),
#' end = c(124914662,50850929, 50850929, 101034277, 151616394,2000)))
#' regulon <- addMotifScore(regulon, peaks=peaks)



addMotifScore <- function(regulon,
                          field_name="motif",
                          peaks=NULL,
                          pwms=NULL,
                          species=c("human","mouse"),
                          genome=c("hg38", "hg19", "mm10"),
                          ...) {

  species <- match.arg(species)
  genome <- match.arg(genome)

  if (!is.null(peaks)&&class(peaks)=="GRanges") {
    if(length(peaks)==0) stop("No peaks provided.")
    message ("annotating peaks with motifs")
    BS.genome <- switch(genome,
                        hg38 = "BSgenome.Hsapiens.UCSC.hg38",
                        hg19 = "BSgenome.Hsapiens.UCSC.hg19",
                        mm10 = "BSgenome.Mmusculus.UCSC.mm10")

    peaks.pruned <- GenomeInfoDb::keepStandardChromosomes(peaks, pruning.mode = "coarse")
    if(length(peaks.pruned)==0) {
      warning("No peaks in standard chromosomes. NAs returned.")
      regulon[,field_name] <- NA
      return(regulon)
    }

    # store original peak indices to match them to regulon idxATAC
    peaks.idx <- GenomicRanges::match(peaks.pruned, peaks)
    peaks.pruned <- peaks.pruned[peaks.idx %in% regulon$idxATAC]
    peaks.idx <- peaks.idx[peaks.idx %in% regulon$idxATAC]
    motifs <- annotateMotif(species, peaks.pruned, BS.genome, pwms, ...)
    motifs <- assay(motifs,"motifMatches")

    # Convert motifs to gene names
    motif_names <- unlist(lapply(strsplit(colnames(motifs), split="_|\\."), "[", 3))

    colnames(motifs) <- matchNames(motif_names, regulon)

  } else {

    stop("supply a GenomicRanges object for peaks")
  }



  # Remove motifs not found in regulon
  motifs <- motifs[, colnames(motifs) %in% unique(regulon$tf), drop=FALSE]

  # Add motif information
  regulon[,field_name] <- NA

  tfs_with_motif <- intersect(colnames(motifs), unique(regulon$tf))

  for (tf in tfs_with_motif){
    regulon[which(regulon$tf ==tf), field_name] <- motifs[match(regulon$idxATAC[which(regulon$tf ==tf)], peaks.idx),tf]
  }

  regulon[,field_name] <- as.numeric(regulon[,field_name])

  regulon


}


annotateMotif <- function(species, peaks, genome, pwms = NULL, ...) {

  if (is.null(pwms)){
    pwms <- species_motif(species)
  }
  motifs <- motifmatchr::matchMotifs(pwms=pwms,
                                     subject=peaks,
                                     genome=genome,
                                     ...)

}

species_motif <- function(species) {
  if(species == "human"){
    bfc <- BiocFileCache::BiocFileCache()
    tryCatch(
      path_to_database <- BiocFileCache::bfcrpath(bfc, "https://github.com/GreenleafLab/chromVARmotifs/raw/master/data/human_pwms_v2.rda"),
      error = function(cond) {
        message(cond)
        stop("Downloading of motif database failed. Please report a bug on github.")
      }
    )
    load(path_to_database)
    return(human_pwms_v2)
  }
  else if (species == "mouse"){
    bfc <- BiocFileCache::BiocFileCache()
    tryCatch(
      path_to_database <- BiocFileCache::bfcrpath(bfc, "https://github.com/GreenleafLab/chromVARmotifs/raw/master/data/mouse_pwms_v2.rda"),
      error = function(cond) {
        message(cond)
        stop("Downloading of motif database failed. Please report a bug on github.")
      }
    )
    load(path_to_database)
    return(mouse_pwms_v2)
  }
  else stop("Species should be 'human' or 'mouse'. Please provide your own pwms for other species")
}

matchNames <- function(motif_names, regulon){
  unique_tf <- data.frame(original = unique(regulon$tf),
                          after = unique(regulon$tf))
  # remove underscore
  unique_tf$after <- gsub("\\_","", unique_tf$original)

  # remove dash
  unique_tf$after <- gsub("\\-","", unique_tf$after)

  # remove dot
  unique_tf$after <- gsub("\\.","", unique_tf$after)


  official <- unique_tf$original[match(motif_names, unique_tf$after)]
  official[which(is.na(official))] <- motif_names[which(is.na(official))]
  official
}
