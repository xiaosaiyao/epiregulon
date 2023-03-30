#' Add Motif Scores
#'
#' @param regulon A DataFrame consisting of tf (regulator) and target in the column names.
#' @param archr_path Character string indicating the path of the ArchR project to retrieve motif information if
#' motif enrichment was already performed. If no motif enrichment has been performed, first annotate the ArchR using
#' `addMotifAnnotations`. If no ArchR project is provided, the user can also provide peaks in the form of GRanges and
#' this function will annotate the peaks with Cisbp
#' @param motif_name Character string	indicating name of the peakAnnotation object (i.e. Motifs) to retrieve from the designated ArchRProject.
#' @param peaks A GRanges object indicating the peaks to perform motif annotation on if ArchR project is not provided.
#' The peak indices should match the `re` column in the regulon
#' @param pwms A PWMatrixList for annotation of motifs using 'motifmatchr::matchMotifs'
#' @param species Character string indicating species. Currently supported species is human or mouse
#' @param genome Character string indicating the genomic build
#'
#' @return A DataFrame with motif matches added with 1s indicating the presence of motifs and
#' 0s indicating the absence of motifs
#' @export
#'
#' @examples
#' regulon <- S4Vectors::DataFrame(tf = c("AR","AR","AR","ESR1","ESR1"),
#' idxATAC = 1:5)
#' peaks <- GRanges(seqnames = c("chr12","chr19","chr19","chr11","chr6"),
#' ranges = IRanges(start = c(124914563,50850844, 50850844, 101034172, 151616327),
#' end = c(124914662,50850929, 50850929, 101034277, 151616394)))
#' regulon <- addMotifScore(regulon, peaks=peaks)



addMotifScore <- function(regulon,
                          archr_path=NULL,
                          field_name="motif",
                          motif_name="Motif",
                          peaks=NULL,
                          pwms=NULL,
                          species=c("human","mouse"),
                          genome=c("hg38","hg19","mm10")) {

  species <- match.arg(species)
  genome <- match.arg(genome)

  if (!is.null(archr_path) & is.null(peaks)) {
    message("retrieving motif information from ArchR project")
    ArchProj <-
      ArchR::loadArchRProject(path = archr_path, showLogo = FALSE)
    matches <- ArchR::getMatches(ArchProj, name = motif_name)
    motifs <- assay(matches, "matches")

    # Convert motifs to gene names
    colnames(motifs) <- unlist(lapply(strsplit(colnames(motifs), split="_"), "[", 1))

  } else if (is.null(archr_path) & !is.null(peaks)) {
    message ("annotating peaks with motifs")
    opts <- list()

    require(chromVARmotifs)
    data("human_pwms_v1")
    data("mouse_pwms_v1")

    species_motif <- function(species) {
      switch(species,
             human=human_pwms_v1,
             mouse=mouse_pwms_v1)
    }

    if (is.null(pwms)){
      pwms <- species_motif(species)
    }
    motifs <- motifmatchr::matchMotifs(pwms=pwms,
                                       subject=peaks,
                                       genome=genome)
    motifs <- assay(motifs,"motifMatches")
    # Convert motifs to gene names
    colnames(motifs) <- lapply(strsplit(colnames(motifs), split="_|\\."), "[", 3)

  } else {

    stop("specify either an ArchR project path OR supply a GenomicRanges object for peaks")
  }



  # Remove motifs not found in regulon
  motifs <- motifs[, colnames(motifs) %in% unique(regulon$tf), drop=FALSE]

  # Add motif information
  regulon[,field_name] <- NA

  tfs_with_motif <- intersect(colnames(motifs), unique(regulon$tf))

  for (tf in tfs_with_motif){
    regulon[which(regulon$tf ==tf), field_name] <- motifs[regulon$idxATAC[which(regulon$tf ==tf)],tf]
  }

  regulon[,field_name] <- as.numeric(regulon[,field_name])

  regulon


}
