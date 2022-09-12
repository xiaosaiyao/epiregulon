#' Calculate the per cell activity of master regulators based on a regulon
#'
#' @param sce A SingleCellExperiment object with rows representing genes and columns represent cells.
#' Rownames (either gene symbols or geneID) must be consistent with the naming convention in the regulon.
#' @param regulon  A data frame consisting of tf (regulator) and target in the column names, with additional columns
#' indicating degree of association between tf and target such as "mor" or "corr" obtained from addWeights.
#' @param mode String indicating the name of column to be used as the weights
#' @param method String indicating the method for calculating activity. Available methods are weightedMean or aucell
#' @param ncore Integer specifying the number of cores to be used in AUCell
#' @param assay String specifying the name of the assay to be retrieved from the SingleCellExperiment object. Set to
#' "logcounts" as the default
#' @param genesets A feature set collection in the form of CompressedSplitDataFrameList that contains genes in the
#' first column and weights in the second column. See [genomitory](
#' http://cedar.gene.com/gran/dev/PkgDocumentation/genomitory/uploads.html#from-compressedsplitdataframelists)
#' for more information on compressedSplitDataFramelists. See details
#' @return A matrix of inferred transcription factor (row) activities in single cells (columns)
#' @export
#' @import methods utils
#' @details
#' This function calculates activity score from a regulon that is a data frame consisting of a tf column,
#' a target column and a weight column. Alternatively, instead of a regulon, this function also accepts weighted
#' signature sets that are stored in Genomitory as CompressedSplitDataFrameList where each gene set or signature
#' is a DataFrame. The user has the option of computing signature score by weighted mean of target gene expression or
#' the relative ranking of the target genes computed by AUCell.
#'
#' @examples
#' # create a mock singleCellExperiment object for gene expression matrix
#' example_sce <- scuttle::mockSCE()
#' example_sce <- scuttle::logNormCounts(example_sce)
#' example_sce$cluster <- sample(LETTERS[1:5], ncol(example_sce), replace = TRUE)
#'
#' # create a mock regulon
#' regulon <- data.frame(tf = c(rep("Gene_0001", 5), rep("Gene_0002", 10)),
#'                       target = c(paste0("Gene_000", 2:6), paste0("Gene_00", 11:20)))
#' regulon.w <- addWeights(regulon, example_sce, cluster_factor = "cluster",
#'                         exprs_values = "counts", min_targets = 5)
#'
#' # calculate activity
#' activity <- calculateActivity(example_sce, regulon.w, assay = "logcounts")
#'
#' @examples
#' \dontrun{
#' # Compute signature scores
#' library(genomitory)
#' breast <- getFeatureSetCollection("GMTY194:analysis/breast.gmt.bz2@REVISION-1")
#' names(breast) <- breast@elementMetadata@listData[["name"]]
#' activity <- calculateActivity(GeneExpressionMatrix, genesets = breast)
#'}
#' @author Xiaosai Yao, Shang-yang Chen

calculateActivity <- function (sce,
                              regulon = NULL,
                              mode = "weight",
                              method = "weightedmean",
                              ncore = 1,
                              assay = "logcounts",
                              genesets = NULL) {
  method <- tolower(method)
  scale.mat <- assay(sce, assay)

  #convert delayedMatrix to dgCMatrix
  if (class(scale.mat)[1] == "DelayedMatrix") {
    writeLines("converting DelayedMatrix to dgCMatrix")
    scale.mat <- as(scale.mat, Class = "dgCMatrix")
  }


  #convert genesets to regulon
  if (!is.null(genesets)){
    if ( class(genesets)[1] == "CompressedSplitDFrameList") {
      regulon <- do.call(rbind,lapply(names(genesets), function(x) {data.frame(tf=x, target=genesets[[x]][,"gene_id"], weight=genesets[[x]][,"weights"])}))
    } else if (class(genesets)[1] == "CompressedCharacterList"){
      regulon <- do.call(rbind,lapply(names(genesets), function(x) {data.frame(tf=x, target=genesets[[x]], weight=1)}))
    }

  }


  #remove genes in regulon not found in sce
  regulon <- regulon[which(regulon$target %in% rownames(scale.mat)),]

  #calculate activity
  if (method == "weightedmean") {
    message(paste("calculating TF activity from regulon using "),
            method)

    tf_indexes <- split(seq_len(nrow(regulon)), regulon$tf)
    unique_tfs <- names(tf_indexes)

    pb <- txtProgressBar(min = 0,
                        max = length(unique_tfs),
                        style = 3)

    score <- vector("list", length(unique_tfs))
    counter <- 0

    for (i in seq_along(unique_tfs)) {

      tf <- unique_tfs[i]
      regulon.current <- regulon[ tf_indexes[[tf]], ]
      geneset <- data.frame(regulon.current$target, regulon.current[, mode])
      score[[i]] <- pathwayscoreCoeffNorm(scale.mat,
                                         geneset,
                                         tf)
      Sys.sleep(1 / 100)
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
    }
    score.combine <- do.call(cbind, score)
    score.combine <- Matrix::t(score.combine)
  }
  else if (method == "aucell") {
    message(paste("calculating TF activity from regulon using "),
            method)
    geneSets <- split(regulon$target, regulon$tf)
    writeLines("ranking cells...")
    cells_rankings <- AUCell::AUCell_buildRankings(scale.mat,
                                                  splitByBlocks=TRUE,
                                                  nCores = ncore,
                                                  plotStats = F)
    writeLines("calculating AUC...")
    cells_AUC <- AUCell::AUCell_calcAUC(geneSets,
                                       rankings = cells_rankings,
                                       nCores = ncore)
    #score.combine = data.frame(AUCell::getAUC(cells_AUC))
    score.combine <- AUCell::getAUC(cells_AUC)
  }
  return(score.combine)
}

#' Aggregate the weighted gene expression of genes belonging to a given pathway or a regulon
#'
#' @param scale.mat A matrix of normalized gene expression with genes as the rows and cells as the columns
#' @param geneset A data frame with the first column indicating the name of the genes in the pathway and the second column indicating the weights of genes contributing to the pathway or regulon
#' @param geneset_name String indicating the name of the pathway or the regulon
#'
#' @return A vector of inferred activity scores for every single cell
#' @export
#' @author Xiaosai Yao

pathwayscoreCoeffNorm <- function(scale.mat, geneset, geneset_name) {

  score <- Matrix::crossprod(scale.mat[geneset[,1], , drop=FALSE], geneset[,2])/nrow(geneset)
  colnames(score) <- geneset_name
  return(score)
}
