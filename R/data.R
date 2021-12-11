#' Datasets containing the genomic coordinates of human TF binding sites
#'
#' A list of GRanges objects containing the genomic coordinates of human TF binding sites
#' curated from Cistrome and ENCODE databases
#'
#' @format A list of GRanges objects
#' \describe{
#'   \item{seqnames}{chromosome}
#'   \item{ranges}{genomic coordinates}
#'   \item{strand}{strand information}
#'   ...
#' }
#' @source \url{http://http://cistrome.org/} \url{https://www.encodeproject.org/}
"human_bed"
