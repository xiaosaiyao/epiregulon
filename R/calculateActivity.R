#' Calculate the per cell activity of master regulators based on a regulon
#'
#' @param expMatrix A SingleCellExperiment object containing gene expression information with rows representing genes and columns represent cells.
#' Rownames (either gene symbols or geneID) must be consistent with the naming convention in the regulon.
#' @param exp_assay String specifying the name of the assay to be retrieved from the SingleCellExperiment object. Set to
#' "logcounts" as the default
#' @param regulon  A DataFrame object consisting of tf (regulator) and target in the column names, with additional columns
#' indicating degree of association between tf and target such as "mor" or "corr" obtained from `addWeights`.
#' @param normalize Logical indicating whether row means should be subtracted from expression matrix. default is FALSE
#' @param mode String indicating the name of column to be used as the weights
#' @param method String indicating the method for calculating activity. Available methods are `weightedMean` or `aucell`
#' @param ncore Integer specifying the number of cores to be used in AUCell
#' @param genesets A list of genesets. Each list element can be a dataframe with the first column indicating the genes and second column indicating the weights.
#' Alternatively, each list element is a character vector corresponding to the genes in the geneset. A feature set collection in the form of CompressedSplitDataFrameList that
#' contains genes in the first column and weights in the second column. See details
#' @param clusters A vector indicating cluster assignment
#' @param FUN function to aggregate the weights
#' @param ... additional arguments for FUN
#' @return A matrix of inferred transcription factor (row) activities in single cells (columns)
#' @export
#' @import methods utils
#' @details
#' This function calculates activity score from a regulon that is a DataFrame consisting of a tf column,
#' a target column and a weight column. Alternatively, instead of a regulon, this function also accepts weighted
#' signature sets where each gene set or signature is a data frame or unweighted signature sets where each gene set is a character vector.
#' The user has the option of computing signature score by weighted mean of target gene expression or
#' the relative ranking of the target genes computed by AUCell.
#'
#' @examples
#' # create a mock singleCellExperiment object for gene expMatrixession matrix
#' set.seed(1000)
#' gene_sce <- scuttle::mockSCE()
#' gene_sce <- scuttle::logNormCounts(gene_sce)
#' rownames(gene_sce) <- paste0("Gene_",1:2000)
#'
#' # create a mock singleCellExperiment object for peak matrix
#' peak_gr <- GRanges(seqnames = "chr1",
#'                    ranges = IRanges(start = seq(from = 1, to = 10000, by = 100), width = 100))
#' peak_counts <- matrix(sample(x = 0:4, size = ncol(gene_sce)*length(peak_gr), replace = TRUE),
#'                       nrow = length(peak_gr), ncol=ncol(gene_sce))
#' peak_sce <- SingleCellExperiment(list(counts = peak_counts), colData = colData(gene_sce))
#' rownames(peak_sce) <- paste0("Peak_",1:100)
#'
#' # create a mock regulon
#' regulon <- data.frame(tf = c(rep("Gene_1",10), rep("Gene_2",10)),
#'                       idxATAC = sample(1:100, 20),
#'                       target = c(paste0("Gene_", sample(3:2000,10)),
#'                                  paste0("Gene_",sample(3:2000,10))))
#'
#' #  # prune regulon
#' pruned.regulon <- pruneRegulon(expMatrix = gene_sce,
#'                                exp_assay = "logcounts",
#'                                peakMatrix = peak_sce,
#'                                peak_assay = "counts",
#'                                regulon = regulon,
#'                                clusters = gene_sce$Treatment,
#'                                regulon_cutoff = 0.5,
#'                                p_adj = TRUE)
#'
#' regulon.w <- addWeights(regulon = regulon,
#'                         expMatrix = gene_sce,
#'                         clusters = gene_sce$Treatment,
#'                         exp_assay = "logcounts",
#'                         min_targets = 5)
#'
#' # calculate activity
#' activity <- calculateActivity(expMatrix = gene_sce,
#'                               regulon = regulon.w,
#'                               exp_assay = "logcounts")
#'
#' # calculate cluster-specific activity if cluster-specific weights are supplied
#' regulon.w$weight <- matrix(runif(nrow(regulon.w)*2, -1,1), nrow(regulon.w),2)
#' colnames(regulon.w$weight) <- c("treat1","treat2")
#'
#' activity.cluster <- calculateActivity(gene_sce,
#' regulon = regulon.w, clusters = gene_sce$Treatment,
#' exp_assay = "logcounts", FUN = colMeans)
#'
#' # compute signature scores from weighted genesets
#' weighted_genesets <- list(set1 = data.frame(genes = c("Gene_1", "Gene_2", "Gene_3"),
#' weights = c(1,2,3)), set2 = data.frame(genes = c("Gene_4", "Gene_5", "Gene_6"), weights = c(4,5,6)))
#'
#' activity <- calculateActivity(gene_sce, genesets = weighted_genesets)
#'
#' # compute signature scores from unweighted genesets
#' unweighted_genesets <- list(set1 = c("Gene_1", "Gene_2", "Gene_3"),
#'                             set2 = c("Gene_4", "Gene_5", "Gene_6"))
#' activity <- calculateActivity(gene_sce, genesets = unweighted_genesets)
#'

#' @author Xiaosai Yao, Shang-yang Chen

calculateActivity <- function (expMatrix = NULL,
                               exp_assay = "logcounts",
                               regulon = NULL,
                               normalize = FALSE,
                               mode = "weight",
                               method = c("weightedmean","aucell"),
                               ncore = 1,
                               genesets = NULL,
                               clusters = NULL,
                               FUN = mean,
                               ...) {
  method <- tolower(method)
  method <- match.arg(method)

  # convert expMatrix to dgCMatrix
  if (checkmate::test_class(expMatrix,classes = "SummarizedExperiment")){
    expMatrix <- assay(expMatrix, exp_assay)
    expMatrix <- as(expMatrix, "dgCMatrix")
  }


  # convert genesets to regulon
  if (!is.null(genesets)){
    if (is.list(genesets)){
      regulon <- genesets2regulon(genesets)
    } else {
      stop("genesets should be a list of data frames or character vectors")
    }
  }

  # check that rownames match regulon
  fraction_genes <- length(which(regulon$target %in% rownames(expMatrix)))/ length(regulon$target)
  if (fraction_genes <  0.01) {
    stop("Less than 1% of target genes in the regulon are found in expression matrix. Check rownames of gene expression matrix ")
  }


  # remove genes in regulons not found in expMatrix
  regulon <- regulon[which(regulon$target %in% rownames(expMatrix)),]

  # calculate activity
  if (method == "weightedmean") {
    message("calculating TF activity from regulon using ", method)


    # if cluster information is provided and if there are cluster-specific weights provided,
    # compute total activity by summation of cluster-specific activity
    if (!is.null(clusters)) {

      # aggregate weights across the same tf-target pairs
      regulon[,mode] <- I(as.matrix(regulon[,mode]))
      aggregated.regulon <- aggregateMatrix(regulon, mode, FUN)

      # create a cluster-specific tf x target weight matrix
      tf_target_mat <-  lapply(sort(unique(clusters)), function(cluster_name) {
        createTfTgMat(aggregated.regulon, mode = paste0("weight.", cluster_name))})

      names(tf_target_mat) <- sort(unique(clusters))

      # if normalize gene expression (taking the mean across all cells)
      if(normalize) meanExpr <- Matrix::rowMeans(expMatrix[rownames(tf_target_mat[[1]]),])


      # Calculate the number of targets per cluster
      freq <- initiateMatCluster(clusters, nrow = length(unique(regulon$tf)))
      rownames(freq) <- unique(regulon$tf)
      for (cluster in colnames(regulon$weight)){
        freq_counts <- table(regulon$tf[which(regulon$weight[, cluster] != 0) ])
        freq[names(freq_counts), cluster] <- freq_counts }


      # Calculating scores
      score.combine <- list()

      for (cluster_name in sort(unique(clusters))){
        score.combine[[cluster_name]] <-
          Matrix::t(expMatrix)[, rownames(tf_target_mat[[cluster_name]]), drop = FALSE] %*%
          tf_target_mat[[cluster_name]]

        #normalize genes
        if(normalize){
          mean_activity <- meanExpr %*% tf_target_mat[[cluster_name]]
          score.combine[[cluster_name]] <- sweep(score.combine[[cluster_name]],
                                                 2, mean_activity, "-")

        }

        # normalize by the number of target genes
        score.combine[[cluster_name]] <- sweep(score.combine[[cluster_name]], 2,
                                               freq[colnames(score.combine[[cluster_name]]), cluster_name], "/")

        # nullify cells not belonging to this cluster
        score.combine[[cluster_name]][which(clusters != cluster_name),] <- 0



      }
      score.combine <- Reduce("+", score.combine)

      score.combine <- Matrix::t(score.combine)


    } else if (is.null(clusters)) {
      # if no cluster information is provided, calculate activity for all cells
      # convert regulon to a matrix of tf * targets for matrix multiplication
      regulon[, mode] <- as.numeric(regulon[, mode])
      aggregated.regulon <- aggregateMatrix(regulon, mode, FUN)

      tf_target_mat <- createTfTgMat(aggregated.regulon, mode)

      # cross product of expMatrix and tf_target matrix
      score.combine <- Matrix::t(expMatrix)[,rownames(tf_target_mat), drop = FALSE] %*%
        tf_target_mat

      # need to normalize
      if(normalize){
        meanExpr <- Matrix::rowMeans(expMatrix[rownames(tf_target_mat),])
        mean_activity <- meanExpr %*% tf_target_mat
        score.combine <- sweep(score.combine, 2, mean_activity, "-")
      }



      #normalize by number of targets
      freq <- table(regulon$tf)
      score.combine <- sweep(score.combine, 2, freq[colnames(score.combine)], "/")

      score.combine <- Matrix::t(score.combine)

    }
  }else if (method == "aucell") {
    message("calculating TF activity from regulon using ", method)
    geneSets <- split(regulon$target, regulon$tf)
    message("ranking cells...")
    cells_rankings <- AUCell::AUCell_buildRankings(expMatrix,
                                                   splitByBlocks = TRUE,
                                                   nCores = ncore,
                                                   plotStats = FALSE)
    message("calculating AUC...")
    cells_AUC <- AUCell::AUCell_calcAUC(geneSets,
                                        rankings = cells_rankings,
                                        nCores = ncore)

    score.combine <- AUCell::getAUC(cells_AUC)
  }
  return(score.combine)
}

genesets2regulon <- function (genesets){
  for (i in seq_len(length(genesets))){
    if ( is(genesets[[i]], "DFrame") | is(genesets[[i]], "data.frame")) {
      genesets[[i]] <- S4Vectors::DataFrame(tf = names(genesets)[i],
                                           target = genesets[[i]][,1],
                                           weight = genesets[[i]][,2])
    } else if (is.vector(genesets[[i]])) {
      genesets[[i]] <- S4Vectors::DataFrame(tf = names(genesets)[i],
                                           target = genesets[[i]],
                                           weight = 1)
    }
  }
  args <- list(make.row.names = FALSE)
  regulon <- do.call(rbind, args = c(genesets,  args))
}


