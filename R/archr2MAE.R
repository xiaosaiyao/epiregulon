#' Create a MultiAssayExperiment from an ArchR project
#'
#' Import the matrices from a ArchR project to into a MultiAssayExperiment,
#' where each ArchR matrix is a SingleCellExperiment.
#' Metadata is not added to the MultiAssayExperiment or Experiments.
#'
#' @param archrDir character string specifying the ArchR project directory,
#'                 which must contain a file called Save-ArchR-Project.rds
#' @param defaultMatrix character string specifying the name of the matrix that
#' unmapped embeddings and reduced dimension matrices should be associated with
#'
#' @return A \linkS4class{MultiAssayExperiment}
#'
#' @author Natalie Fox
#'
#' @import MultiAssayExperiment
#'
#' @export
#'
archr2MAE <- function(archrDir, defaultMatrix = "TileMatrix") {
  checkmate::assertDirectoryExists(archrDir, access = "r")
  checkmate::assertFileExists(file.path(archrDir, "Save-ArchR-Project.rds"), access = "r")
  checkmate::assertString(defaultMatrix)

  # Load ArchR project
  archr.logging <- ArchR::getArchRLogging()
  ArchR::addArchRLogging(useLogs = archr.logging)

  # Check that there is a saved ArchR project in the directory
  if (! "Save-ArchR-Project.rds" %in% list.files(archrDir)) {
    stop("No Save-ArchR-Project.rds file found in the ArchR project directory.",
         "Are you sure this is an ArchR project directory?")
  }

  archrProj <- suppressMessages(ArchR::loadArchRProject(archrDir))

  # Check that defaultMatrix is found in ArchR project
  if (!defaultMatrix %in% ArchR::getAvailableMatrices(archrProj)) {
    stop(defaultMatrix, " not found in ArchR project")
  }


  # Get list of Single Cell Experiments
  all.exp <- create.exp.list.from.archr(archrProj  = archrProj, defaultMatrix = defaultMatrix)

  # create the sample Map for the MultiAssayExperiment
  el <- ExperimentList(all.exp)
  maplist <- lapply(all.exp, function(se) {
    data.frame(primary = se$Sample, colname = colnames(se), stringsAsFactors = FALSE)})

  sampMap <- listToMap(maplist)

  # create the MultiAssayExperiment
  mae <- MultiAssayExperiment(el,
                              sampleMap = sampMap,
                              colData = S4Vectors::DataFrame(row.names=unique(sampMap$primary)));

  return(mae)
}



#' @importFrom rhdf5 h5closeAll h5ls h5read
#' @keywords internal
create.exp.list.from.archr <- function(archrProj, defaultMatrix = "TileMatrix") {

  # map which Embeddings should be saved to which SingleCellExperiment
  embedding.map <- create.embedding.map(archrProj, defaultMatrix = defaultMatrix)

  # import the existing matrices into Experiment objects
  all.exp <- list()
  available.matrices <- ArchR::getAvailableMatrices(archrProj)

  for (matrix.type in available.matrices) {
    all.exp[[matrix.type]] <- load.matrix(matrix.type, archrProj, embedding.map, defaultMatrix = defaultMatrix)
  }


  return(all.exp)
}



#' @keywords internal
create.embedding.map <- function(archrProj, defaultMatrix = "TileMatrix") {

  # map which matrix was used to create each reduced dimension set
  reduced.dim.matrix.type <- lapply(archrProj@reducedDims, function(x) x$useMatrix)

  # map which reduced dim set was used to create each embedding
  embedding.reduced.dim.mapping <- lapply(
    names(archrProj@embeddings), function(x)
      unique(sub("#.*", "", colnames(ArchR::getEmbedding(archrProj, x))))
  )
  names(embedding.reduced.dim.mapping) <- names(archrProj@embeddings)

  # map which matrix each embedding was created from
  embedding.map <- list()
  for (embedding.name in names(embedding.reduced.dim.mapping)) {

    # match embeddings without a specified matrix to TileMatrix
    matrix.name <- defaultMatrix
    reduced.dim.name <- embedding.reduced.dim.mapping[[embedding.name]]
    if (reduced.dim.name %in% names(reduced.dim.matrix.type)) {
      if (!is.null(reduced.dim.matrix.type[[reduced.dim.name]])){
        matrix.name <- reduced.dim.matrix.type[[reduced.dim.name]]
      }

    }

    embedding.map[[matrix.name]] <- c(embedding.map[[matrix.name]], embedding.name)
  }

  return(embedding.map)
}





#' @import SummarizedExperiment
#' @import SingleCellExperiment
#' @keywords internal
load.matrix <- function(matrix.type, archrProj, embedding.map = NULL, defaultMatrix) {
  message(matrix.type)

  # check if the matrix was created as binary or not
  matrixClass <- rhdf5::h5read(ArchR::getArrowFiles(archrProj)[1], paste0(matrix.type, "/Info/Class"))
  binarize <- (matrixClass == "Sparse.Binary.Matrix")

  # Get the SummarizedExperimets with one column per cell and then convert to SingleCellExperiments
  se <- suppressMessages(
    ArchR::getMatrixFromProject(archrProj, matrix.type, binarize = binarize, threads = 1))

  # check that colData colnames match the getCellColData colnames and replace if not matching
  anExpression <- "projColData[, colnames(projColData) %ni% colnames(colData)]"
  if (anExpression %in% colnames(colData(se))) {
    cell.col.data <- ArchR::getCellColData(archrProj)
    missing.colname <- setdiff(colnames(cell.col.data), colnames(colData(se)))
    if (length(missing.colname) == 1 && length(which(colnames(colData(se)) == anExpression)) == 1) {
      colnames(colData(se))[which(colnames(colData(se)) == anExpression)] <- missing.colname
    }
  }

  # Create the rowRanges for the experiment
  if (matrix.type == "TileMatrix") {
    se <- assign.tile.rowranges(se,chrSizes = ArchR::getGenomeAnnotation(archrProj)$chromSizes)
  } else if (is.null(rowRanges(se)) & all(c("start", "end", "strand") %in% names(rowData(se)))) {
    row.data <- rowData(se)
    row.data$start[row.data$strand == 2] <- rowData(se)$end[row.data$strand == 2]
    row.data$end[row.data$strand == 2] <- rowData(se)$start[row.data$strand == 2]
    rowRanges(se) <- GenomeInfoDb::sortSeqlevels(GenomicRanges::GRanges(
      seqnames = row.data$seqnames,
      ranges = IRanges::IRanges(
        start = row.data$start,
        end = row.data$end),
      strand = c("+", "-", "*")[as.numeric(row.data$strand)]
    ))
    # remove some columns
    rowData(se) <- row.data[, !names(row.data) %in% c("start", "end", "seqnames", "strand")]
  }

  sce <- methods::as(se, "SingleCellExperiment")

  # change the assay Name
  assayNames(sce)[grep("GeneExpressionMatrix|GeneIntegrationMatrix", assayNames(sce))] <- "logcounts"
  assayNames(sce)[grep("PeakMatrix|TileMatrix", assayNames(sce))] <- "counts"


  # add reduced dimensions
  reduced.dims.list <- S4Vectors::SimpleList()

  ## determine matrix types
  reduced.dim.matrix.type <- unlist(lapply(
    archrProj@reducedDims,
    function(x) {
      id <- x$useMatrix
      if (is.null(id)) {
        id <- defaultMatrix
      }
      id
    }))


  for (rd.name in names(reduced.dim.matrix.type[reduced.dim.matrix.type == matrix.type])) {
    reduced.dims.list[[rd.name]] <- ArchR::getReducedDims(archrProj, rd.name)
  }

  # add reduced dimensions - Embeddings
  if (matrix.type %in% names(embedding.map)) {
    embedding.list <- S4Vectors::SimpleList()
    for (embedding.name in embedding.map[[matrix.type]]) {
      embedding.list[[embedding.name]] <- ArchR::getEmbedding(archrProj, embedding.name)
    }
    reduced.dims.list <- append(reduced.dims.list, embedding.list)
  }
  if (length(reduced.dims.list) > 0) {
    for (rd.name in names(reduced.dims.list)) {
      if (all(rownames(reduced.dims.list[[rd.name]]) %in% colnames(sce))) {
        reduced.dims.list[[rd.name]] <- reduced.dims.list[[rd.name]][colnames(sce), ]
      }
    }

    reducedDims(sce) <- reduced.dims.list
  }

  # remove <Rle> from colData
  colData(sce)$Sample <- as.character(colData(sce)$Sample)

  return(sce)
}


#' @import SummarizedExperiment
#' @keywords internal
assign.tile.rowranges <- function(se, chrSizes) {
  # Create the GRanges for the tile matrix
  se.row.data <- rowData(se)
  tile.size <- se.row.data$start[2] - se.row.data$start[1]
  tile.size.check <- se.row.data$start[3] - se.row.data$start[2]
  if (tile.size != tile.size.check) {
    stop('not sure of the tile size')
  }
  se.tile.size <- tile.size
  if (length(tile.size) == 1) {
    rowRanges(se) <- GenomicRanges::GRanges(
      seqnames = rowData(se)$seqnames,
      ranges = IRanges::IRanges(
        start = rowData(se)$start + 1,
        end = rowData(se)$start + tile.size
      )
    )
    # adjust the tiles on the end of the chromosomes to match the chromosome ends
    for (chr in GenomeInfoDb::seqlevels(rowRanges(se))) {
      chr.end <- BiocGenerics::end(chrSizes)[as.logical(seqnames(chrSizes) == chr)]
      BiocGenerics::end(rowRanges(se))[as.logical(seqnames(rowRanges(se)) == chr) &
                                         BiocGenerics::end(rowRanges(se)) > chr.end] <- chr.end
    }
    rowRanges(se) <- GenomeInfoDb::sortSeqlevels(rowRanges(se))
    se <- se[order(rowRanges(se)), ]
    rownames(se) <- paste0(seqnames(rowRanges(se)), ":",
                           BiocGenerics::start(rowRanges(se)), "-", BiocGenerics::end(rowRanges(se)))
  }

  return(se)
}

