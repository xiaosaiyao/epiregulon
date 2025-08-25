#' custom classes
#'
#' Additional class definitions.
#'
#' @section Classes:
#' \itemize{
#' \item{\code{SingleCellAccessibilityExperiment}, contains \code{SingleCellExperiment}}
#' }
#'
#' @keywords internal
#'
#' @name customClasses
#'
#' @rdname customClasses
setClass("SingleCellAccessibilityExperiment", contains = "SingleCellExperiment")
