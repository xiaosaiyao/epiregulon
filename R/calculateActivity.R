#' Calculate the per cell activity of master regulators based on a regulon
#'
#' @param sce A SingleCellExperiment object with rows representing genes and columns represent cells.
#' Rownames (either gene symbols or geneID) must be consistent with the naming convention in the regulon.
#' @param regulon  A data frame consisting of tf (regulator) and target in the column names, with additional columns
#' indicating degree of association between tf and target such as "mor" or "corr" obtained from addWeights.
#' @param normalize Logical indicating whether row means should be substracted from expression matrix
#' @param mode String indicating the name of column to be used as the weights
#' @param method String indicating the method for calculating activity. Available methods are weightedMean or aucell
#' @param ncore Integer specifying the number of cores to be used in AUCell
#' @param assay String specifying the name of the assay to be retrieved from the SingleCellExperiment object. Set to
#' "logcounts" as the default
#' @param genesets A feature set collection in the form of CompressedSplitDataFrameList that contains genes in the
#' first column and weights in the second column. See [genomitory](
#' http://cedar.gene.com/gran/dev/PkgDocumentation/genomitory/uploads.html#from-compressedsplitdataframelists)
#' for more information on compressedSplitDataFramelists. See details
#' @param clusters A vector indicating cluster assignment
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
#' # calculate cluster-specific activity if cluster-specific weights are supplied
#' regulon.w$weight_positive = runif(nrow(regulon.w), -1,1)
#' regulon.w$weight_negative = runif(nrow(regulon.w), -1,1)
#' activity.cluster <- calculateActivity(example_sce, regulon=regulon.w,
#' cluster=example_sce$Mutation_Status, assay = "logcounts")
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
                               normalize = FALSE,
                               mode = "weight",
                               method = "weightedmean",
                               ncore = 1,
                               assay = "logcounts",
                               genesets = NULL,
                               clusters = NULL) {
  method <- tolower(method)
  scale.mat <- assay(sce, assay)

  #convert delayedMatrix to dgCMatrix
  if (checkmate::test_class(scale.mat, classes = "DelayedMatrix")) {
    writeLines("converting DelayedMatrix to dgCMatrix")
    scale.mat <- as(scale.mat, Class = "dgCMatrix")
  }


  #convert genesets to regulon
  if (!is.null(genesets)){
    if ( checkmate::test_class(genesets, classes =  "CompressedSplitDFrameList")) {
      regulon <- do.call(rbind,lapply(names(genesets), function(x) {data.frame(tf=x, target=genesets[[x]][,"gene_id"], weight=genesets[[x]][,"weights"])}))
    } else if ( checkmate::test_class(genesets, classes =  "CompressedCharacterList")) {
      regulon <- do.call(rbind,lapply(names(genesets), function(x) {data.frame(tf=x, target=genesets[[x]], weight=1)}))
    }

  }


  #remove genes in regulon not found in sce
  regulon <- regulon[which(regulon$target %in% rownames(scale.mat)),]

  #normalize genes
  if (normalize){
    rowmeans <- Matrix::rowMeans(scale.mat)
    scale.mat <- apply(scale.mat, 2, function(x){x-rowmeans})
  }

  #calculate activity
  if (method == "weightedmean") {
    message(paste("calculating TF activity from regulon using "), method)


    # if cluster information is provided and if there are cluster-specific weights provided,
    # compute total activity by summation of cluster-specific activity
    if (!is.null(clusters)) {

      tf_target_mat <-
        lapply(sort(unique(clusters)), function(cluster_name) {
          local({
            # convert regulon to a matrix of tf * targets for matrix multiplication
            tf_target_mat <- reshape2::dcast(regulon,
                                             target ~ tf,
                                             fun.aggregate = mean,
                                             value.var = paste0(mode, "_", cluster_name))
            rownames(tf_target_mat) <- tf_target_mat$target
            tf_target_mat <- tf_target_mat[, -1]
            # convert NA to 0
            tf_target_mat[is.na(tf_target_mat)] <- 0
            tf_target_mat <- as(as.matrix(tf_target_mat), "dgCMatrix")
          })
        })
      names(tf_target_mat) <- sort(unique(clusters))


      score.combine <- list()
      for (cluster_name in sort(unique(clusters))){
        score.combine[[cluster_name]] <-
          Matrix::t(scale.mat)[, rownames(tf_target_mat[[cluster_name]]), drop = FALSE] %*%
                            tf_target_mat[[cluster_name]]
        # nullify cells not belonging to this cluster
        score.combine[[cluster_name]][which(clusters != cluster_name),] <- 0

      }
      score.combine <- Reduce("+", score.combine)


    } else {
      # if no cluster information is provided, calculate activity for all cells
      # convert regulon to a matrix of tf * targets for matrix multiplication
      tf_target_mat <- reshape2::dcast(regulon, target ~ tf, fun.aggregate = mean, value.var = mode)
      rownames(tf_target_mat) <- tf_target_mat$target
      tf_target_mat <- tf_target_mat[,-1]
      # convert NA to 0
      tf_target_mat[is.na(tf_target_mat)] <- 0
      tf_target_mat <- as(as.matrix(tf_target_mat), "dgCMatrix")
      # cross product of scale.matrix and tf_target matrix
      score.combine <- Matrix::t(scale.mat)[,rownames(tf_target_mat), drop = FALSE] %*%
        tf_target_mat
      # need to normalize
    }
    score.combine <- Matrix::t(score.combine)
    #normalize by number of targets
    freq <- table(regulon$tf)
    score.combine <- apply(score.combine, 2, function(x) {x/freq[rownames(score.combine)]})

  }
  else if (method == "aucell") {
    message(paste("calculating TF activity from regulon using "), method)
    geneSets <- split(regulon$target, regulon$tf)
    writeLines("ranking cells...")
    cells_rankings <- AUCell::AUCell_buildRankings(scale.mat,
                                                   splitByBlocks = TRUE,
                                                   nCores = ncore,
                                                   plotStats = FALSE)
    writeLines("calculating AUC...")
    cells_AUC <- AUCell::AUCell_calcAUC(geneSets,
                                        rankings = cells_rankings,
                                        nCores = ncore)
    #score.combine = data.frame(AUCell::getAUC(cells_AUC))
    score.combine <- AUCell::getAUC(cells_AUC)
  }
  return(score.combine)
}
