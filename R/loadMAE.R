
#' check if a file is a valid HDF5 file
#'
#' Check if a file path argument points to an non-corrupt HDF5 file.
#'
#' Compares the first 8 bytes of a file to those of the standard HDF5 file header.
#'
#' @param path path to file to test
#'
#' @return
#' Returns invisible \code{path} if check is successful, otherwise signals an error.
#'
#' @section Further development:
#' The HDF5 file header contains 8 bytes, which hold specific meanings.
#' Currently the function only tests that the header of the file
#' specified by \code{path} is identical to a healthy HDF5 file and signals
#' a general error if that is not the case.
#' Reporting specific types of corruption can be implemented.
#'
#' @references
#' <http://web.ics.purdue.edu/~aai/HDF5/html/H5.format.html#BootBlock>
#'
#' @author Aleksander Chlebowski
assertHDF5 <- function(path) {

    fileHead <- readBin(path, "raw", n = 8L)
    hdf5Head <- as.raw(c(137L, 72L, 68L, 70L, 13L, 10L, 26L, 10L))

    ident <- identical(fileHead, hdf5Head)

    if (ident) {
        return(invisible(path))
    } else {
        stop("Assertion on \"path\" failed: Must be a hdf5 file.", call. = FALSE)
    }
}


#' @import MultiAssayExperiment
#'
loadMAE <- function(file, experiments, verbose) {
    checkmate::assertFileExists(file, access = "r")
    assertHDF5(file)
    checkmate::assertCharacter(experiments)
    checkmate::assertFlag(verbose)

    if (verbose) {
        dataset <- dynGet("dataset",
                          ifnotfound = paste("stored in file",
                                             dynGet("file",
                                                    ifnotfound = stop("resource missing"))))
        message("loading dataset ", dataset)
    }

    # list experiments in file
    fileContents <- rhdf5::h5ls(file)
    expStored <- fileContents[fileContents[["group"]] == "/", "name"]
    lapply(experiments, checkmate::assertChoice, choices = expStored)

    # load files
    if (verbose) message("loading experiments")
    expList <- lapply(experiments, loadExp, file = file, verbose = verbose)
    names(expList) <- experiments

    # build mae
    if (verbose) message("building MultiAssayExperiment")
    mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = expList)

    return(mae)
}

loadExp <- function(file, expName, verbose) {
    checkmate::assertFileExists(file, access = "r")
    assertHDF5(file)
    checkmate::assertString(expName)
    checkmate::assertFlag(verbose)

    # load class
    expClass <- rhdf5::h5read(file, sprintf("%s/class", expName))

    if (verbose) message("loading ", expName)

    # list file contents
    fileContents <- rhdf5::h5ls(file)
    colnames.present <- is.element("colnames",
                                   fileContents[fileContents[["group"]] == sprintf("/%s/properties", expName), "name"])
    colData.present <- is.element("colData",
                                  fileContents[fileContents[["group"]] == sprintf("/%s/properties", expName), "name"])
    rownames.present <- is.element("rownames",
                                   fileContents[fileContents[["group"]] == sprintf("/%s/properties", expName), "name"])
    rowData.present <- is.element("rowData",
                                  fileContents[fileContents[["group"]] == sprintf("/%s/properties", expName), "name"])
    rowRanges.present <- is.element("rowRanges",
                                    fileContents[fileContents[["group"]] == sprintf("/%s/properties", expName), "name"])
    metadata.present <- is.element("metadata",
                                   fileContents[fileContents[["group"]] == sprintf("/%s/properties", expName), "name"])
    reducedDims.present <- is.element("reducedDims",
                                      fileContents[fileContents[["group"]] == sprintf("/%s", expName), "name"])
    altExps.present <- is.element("altExps",
                                  fileContents[fileContents[["group"]] == sprintf("/%s", expName), "name"])


    # load everything in group specified by experiment

    # load properties
    if (verbose) message("\t loading properties")

    # cell information
    if (colnames.present) {
        if (verbose) message("\t ... colnames")
        colnames <- rhdf5::h5read(file, sprintf("%s/properties/colnames", expName))
    }
    if (colData.present) {
        if (verbose) message("\t ... colData")
        colData <- rhdf5::h5read(file, sprintf("%s/properties/colData", expName))
    }

    # feature information
    if (rownames.present) {
        if (verbose) message("\t ... rownames")
        rownames <- rhdf5::h5read(file, sprintf("%s/properties/rownames", expName))
    }
    if (rowData.present) {
        if (verbose) message("\t ... rowData")
        rowData <- rhdf5::h5read(file, sprintf("%s/properties/rowData", expName))
    }

    # load row ranges if saved
    if (rowRanges.present) {
        if (verbose) message("\t ... rowRanges")
        rowRanges <- rhdf5::h5read(file = file, name = sprintf("%s/properties/rowRanges", expName))
        rowRanges <- restoreGR(rowRanges)
    }

    # load metadata if saved
    if (metadata.present) {
        if (verbose) message("\t ... metadata")
        metadata <- rhdf5::h5read(file = file, name = sprintf("%s/properties/metadata", expName))
        metadata <- eval(parse(text = metadata))
    }

    # load assays
    if (verbose) message("\t loading assays")
    assaysStored <- fileContents[fileContents[["group"]] == sprintf("/%s/assays", expName), "name"]
    assays <- lapply(assaysStored, function(ass) {
        if (verbose) message("\t ... ", ass)
        HDF5Array::H5SparseMatrix(filepath = file, group = sprintf("/%s/assays/%s", expName, ass))
    })
    names(assays) <- assaysStored

    # load reduced dimensions if saved
    if (reducedDims.present) {
        if (verbose) message("\t loading reducedDims")
        reducedDimsStored <- fileContents[fileContents[["group"]] == sprintf("/%s/reducedDims", expName), "name"]
        reducedDims <- lapply(reducedDimsStored, function(red) {
            if (verbose) message("\t ... ", red)
            rhdf5::h5read(file = file, name = sprintf("%s/reducedDims/%s", expName, red))
        })
        names(reducedDims) <- reducedDimsStored
    }

    # load alternative experiment
    if (altExps.present) {
        if (verbose) message("\t loading altExps")
        altExpsStored <- fileContents[fileContents[["group"]] == sprintf("/%s/altExps", expName), "name"]
        altExps <- lapply(altExpsStored, function(alt) {
            if (verbose) message("\t ... ", alt)
            loadExp(file = file,
                    expName = sprintf("%s/altExps/%s", expName, alt),
                    verbose = verbose)
        })
        names(altExps) <- altExpsStored
    }


    # rebuild experiment
    if (verbose) message("\t rebuilding experiment")
    ans <- SummarizedExperiment::SummarizedExperiment(assays = assays)
    if (colData.present) {
        if (verbose) message("\t ...colData")
        SummarizedExperiment::colData(ans) <- S4Vectors::DataFrame(colData)
    }
    if (colnames.present) {
        if (verbose) message("\t ...colnames")
        colnames(ans) <- colnames
    }
    if (rowData.present) {
        if (verbose) message("\t ...rowData")
        SummarizedExperiment::rowData(ans) <- S4Vectors::DataFrame(rowData)
    }
    if (rowRanges.present) {
        if (verbose) message("\t ...rowRanges")
        SummarizedExperiment::rowRanges(ans) <- rowRanges
    }
    if (rownames.present) {
        if (verbose) message("\t ...rownames")
        rownames(ans) <- rownames
    }
    if (metadata.present) {
        if (verbose) message("\t ...metadata")
        S4Vectors::metadata(ans) <- metadata
    }

    # if inherits from SCE, convert and add more slots
    if (methods::extends(expClass, "SingleCellExperiment")) {
        # direct conversion from (Ranged)SummarizedExperiment fails
        # due to different requirements for the int_elementMetadata slot
        ans <- methods::as(ans, "SingleCellExperiment")
        if (!methods::is(ans, expClass)) {
            ans <- methods::as(ans, expClass)
        }

        if (reducedDims.present) {
            if (verbose) message("\t ...reducedDims")
            SingleCellExperiment::reducedDims(ans, withDimnames = FALSE) <- reducedDims
        }
        if (altExps.present) {
            if (verbose) message("\t ...altExps")
            SingleCellExperiment::altExps(ans, withDimnames = FALSE) <- altExps
        }
    }

    return(ans)
}

#' @keywords internal
#'
restoreGR <- function(df) {
    checkmate::assertDataFrame(df)
    basic <- c("seqnames", "start", "end", "width", "strand")
    df[basic] <- lapply(df[basic], function(x) methods::as(x, typeof(x)))
    ans <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)

    return(ans)
}
