addWeights=function(regulon, sce, cluster_factor, block_factor= NULL, exprs_values=NULL,  corr=TRUE, MI=TRUE, multicore=TRUE){
  require(scater)

  if (is.null(exprs_values)){
    exprs_values="logcounts"
  }

  if (isTRUE(multicore)){
    require(BiocParallel)
    BPPARAM <- MulticoreParam()
  } else {
    BPPARAM <- MulticoreParam(workers=1)
  }

  #compute average expression across clusters and batches
  message("calculating average expression across clusters...")
  if (is.null(block_factor)){
    averages.se <- scater::sumCountsAcrossCells(sce, exprs_values=exprs_values, ids=colData(sce)[cluster_factor], average=T, BPPARAM=BPPARAM)

  } else {
    groupings = DataFrame(cluster=colData(sce)[cluster_factor], block=colData(sce)[block_factor])
    averages.se = scater::sumCountsAcrossCells(sce, exprs_values=exprs_values, groupings, average=T, BPPARAM=BPPARAM)
  }

  expr=assays(averages.se)$average

  if (isTRUE(corr)){
    message("computing correlation of the regulon...")
    # adding correlation
    regulon$corr = 0

    pb = txtProgressBar(min = 0, max = nrow(regulon), style = 3)
    for (i in 1:nrow(regulon)){
      tf_index = match(regulon$tf[i], rownames(expr))
      target_index = match(regulon$target[i], rownames(expr))
      if (is.na(tf_index) | is.na(target_index)){
        regulon$corr[i] = NA}
      else{
        corr = cor(expr[tf_index,], expr[target_index,], use = "na.or.complete")
        regulon$corr[i] = corr}
      Sys.sleep(1 / 100)
      setTxtProgressBar(pb, i)
    }
  }
  #adding mutual information

  if (isTRUE(MI)) {
    message("computing mutual information of the regulon...")
    require(entropy)
    regulon$MI = 0

    pb = txtProgressBar(min = 0, max = nrow(regulon), style = 3)
    for (i in 1:nrow(regulon)){
      tf_index = match(regulon$tf[i], rownames(expr))
      target_index = match(regulon$target[i], rownames(expr))
      if (length(unique(expr[tf_index,]))<10 | length(unique(expr[target_index,]))<10 ){
        regulon$MI[i] = NA}
      else{
        y2d = discretize2d(expr[tf_index,], expr[target_index,], numBins1=10, numBins2=10)
        mi = mi.empirical(y2d)
        regulon$MI[i] = mi}
      Sys.sleep(1 / 100)
      setTxtProgressBar(pb, i)
    }
  }
  return(regulon)
}
