#' A function to calculate weights for regulons
#'
#' @param regulon a regulon in the form of a tall format matrix consisting of tf(regulator), target and a column indicating degree of association between TF and target such as "mor" or "corr".
#'           example regulon:
#'           tf      target  corr
#' @param sce scater single cell object compute weights on
#' @param cluster_factor cluster labels of the scater single cell object
#' @param block_factor batch labels of the scater single cell object
#' @param exprs_values name of the assay to be retrieved from scater object
#' @param corr whether to calculate weights based on correlation
#' @param MI whether to calculate weights based on mutual information
#' @param multicore number of cores to use
#'
#' @return dataframe with three columns in the form of regulon input. The third column is replaced with the inferred weights
#' @importFrom stats cor
#' @importFrom SummarizedExperiment assays colData
#' @export
#'
#' @examples
#' # load scRNA-seq data to infer activity on
#' sce <- readRDS("/gstore/project/lineage/sam/heme_GRN/scRNA-Granja-2019.rds")
#' # regulon matrix from getRegulon function
#' regulon.w=addWeights(regulon=regulon, sce=sce, cluster_factor="BioClassification", block_factor=NULL, corr=TRUE, MI=FALSE, multicore=TRUE)

addWeights=function(regulon, sce, cluster_factor, block_factor= NULL, exprs_values=NULL,  corr=TRUE, MI=TRUE, multicore=TRUE){

  if (is.null(exprs_values)){
    exprs_values="logcounts"
  }

  if (isTRUE(multicore)){

    BPPARAM <- BiocParallel::MulticoreParam()
  } else {
    BPPARAM <- BiocParallel::MulticoreParam(workers=1)
  }

  #compute average expression across clusters and batches
  message("calculating average expression across clusters...")
  if (is.null(block_factor)){
    averages.se <- scater::sumCountsAcrossCells(sce, exprs_values=exprs_values, ids=colData(sce)[cluster_factor], average=T, BPPARAM=BPPARAM)

  } else {
    groupings = data.frame(cluster=colData(sce)[cluster_factor], block=colData(sce)[block_factor])
    averages.se = scater::sumCountsAcrossCells(sce, exprs_values=exprs_values, groupings, average=T, BPPARAM=BPPARAM)
  }

  expr=assays(averages.se)$average

  if (isTRUE(corr)){
    message("computing correlation of the regulon...")

    unique_tfs = unique(regulon$tf)

    pb = txtProgressBar(min = 0, max = length(unique_tfs), style = 3)

    regulon_list <- list()

    for (i in 1:length(unique_tfs)){

      tf_regulon = subset(regulon, tf==unique_tfs[i])
      tf_index = match(tf_regulon$tf[1], rownames(expr))

      weights = apply(tf_regulon, 1, function(tf_target_pair){

        target_index = match(tf_target_pair['target'], rownames(expr))

        if (is.na(tf_index) | is.na(target_index)){

          return(NA)

        } else {

          corr = cor(expr[tf_index,], expr[target_index,], use = "na.or.complete")
          return(corr)

        }
      })

      tf_regulon$weight = weights
      regulon_list[[i]] <- tf_regulon

      Sys.sleep(1 / 100)
      setTxtProgressBar(pb, i)

    }
  }

  #adding mutual information

  if (isTRUE(MI)) {
    message("computing mutual information of the regulon...")
    regulon$MI = 0

    pb = txtProgressBar(min = 0, max = nrow(regulon), style = 3)
    for (i in 1:nrow(regulon)){
      tf_index = match(regulon$tf[i], rownames(expr))
      target_index = match(regulon$target[i], rownames(expr))
      if (length(unique(expr[tf_index,]))<10 | length(unique(expr[target_index,]))<10 ){
        regulon$MI[i] = NA}
      else{
        y2d = entropy::discretize2d(expr[tf_index,], expr[target_index,], numBins1=10, numBins2=10)
        mi = entropy::mi.empirical(y2d)
        regulon$MI[i] = mi}
      Sys.sleep(1 / 100)
      setTxtProgressBar(pb, i)
    }
  }

  regulon = do.call("rbind", regulon_list)
  return(regulon)
}
