#' A function that calculates the per cell activity of master regulons based on a regulon
#'
#' @param scale.mat a matrix of genes by cells. Rows represent genes and columns represent cells. Rownames (either gene symbols or geneID) must be consistent with the naming convention in the regulon.
#' @param regulon a regulon in the form of a tall format matrix consisting of tf(regulator), target and a column indicating degree of association between TF and target such as "mor" or "corr".
#'           example regulon:
#'           tf      target  corr
#'          Esr1    Pgr     0.56
#' @param mode a string indicating the mode of regulon such as "mor" or "corr" when weightedMean is the chosen method
#' @param method method for calculating activity. Available methods are weightedMean or aucell. The parameter used for the weights in weightedMean is further specified by mode.
#' @param ncore number of cores to use
#'
#' @return a matrix of inferred transcription factor (row) activities in single cells (columns)
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom magrittr %>%
#'
#' @examples 1+1
calculateActivity=function(scale.mat, regulon, mode, method=NULL, ncore=NULL){
  method=tolower(method)

  message(method)
  if (is.null(ncore)){
    ncore=1
  }
  if (method == "weightedmean"){
    message(paste("calculating TF activity from regulon using "), method)

    TFs.found = unique(regulon$tf)
    pb = txtProgressBar(min = 0, max = length(TFs.found), style = 3)
    score=list()
    for (i in 1:length(TFs.found)) {
      regulon.current=regulon[regulon$tf==TFs.found[i],]
      geneset=data.frame(regulon.current$target,regulon.current[,mode])
      score[[TFs.found[i]]]=pathwayscoreCoeffNorm(scale.mat, rownames(scale.mat), geneset, TFs.found[i])
      Sys.sleep(1 / 100)
      setTxtProgressBar(pb, i)
    }

    score.combine=dplyr::bind_cols(score) %>% t() %>% as.data.frame()

  } else if (method == "aucell"){
    #requireNamespace(AUCell)
    message(paste("calculating TF activity from regulon using "), method)

    geneSets = split(regulon$target, regulon$tf)
    message("ranking cells...")
    cells_rankings = AUCell::AUCell_buildRankings(scale.mat, nCores=ncore, plotStats=F)
    message("calculating AUC...")
    cells_AUC = AUCell::AUCell_calcAUC(geneSets, rankings=cells_rankings, nCores = ncore)
    score.combine=data.frame(AUCell::getAUC(cells_AUC))
  }

  return(score.combine)


}

pathwayscoreCoeffNorm=function(mat_scale, colnames, pathway, tf.name){
  pathway_index=match(pathway[,1], colnames)
  if (length(pathway_index)==1){
    score=as.data.frame(mat_scale[pathway_index,] * pathway[,2],na.rm=T)
  } else {
    score=as.data.frame(colMeans(mat_scale[pathway_index,] * pathway[,2],na.rm=T))
  }
  colnames(score) <- tf.name
  return(score)
}
