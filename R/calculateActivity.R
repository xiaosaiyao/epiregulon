#' A function that calculates the per cell activity of master regulators based on a regulon
#'
#' @param sce A SingleCellExperiment object with rows representing genes and columns represent cells. Rownames (either gene symbols or geneID) must be consistent with the naming convention in the regulon.
#' @param regulon  A data frame consisting of tf (regulator) and target in the column names, with additional columns indicating degree of association between tf and target such as "mor" or "corr" obtained from addWeights.
#' @param mode String indicating the name of column to be used as the weights
#' @param method String indicating the method for calculating activity. Available methods are weightedMean or aucell
#' @param ncore Integer specifying the number of cores to be used in AUCell
#' @param assay String specifying the name of the assay to be retrieved from the SingleCellExperiment object. Set to "logcounts" as the default
#' @return a matrix of inferred transcription factor (row) activities in single cells (columns)
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @examples
#' # regulon matrix from addRegulon or addWeights functinos
#' score.combine <- calculateActivity(sce, regulon.w, "weight", method="weightedMean", assay="logcounts")

calculateActivity = function (sce,
                              regulon,
                              mode = "weight",
                              method = "weightedmean",
                              ncore=1,
                              assay = "logcounts") {
  method = tolower(method)
  scale.mat = as.matrix(assay(sce, assay))

  if (method == "weightedmean") {
    message(paste("calculating TF activity from regulon using "),
            method)
    unique_tfs = unique(regulon$tf)
    pb = txtProgressBar(min = 0,
                        max = length(unique_tfs),
                        style = 3)
    score = list()
    counter = 0
    for (tf in unique_tfs) {
      regulon.current = regulon[regulon$tf ==  tf, ]
      geneset = data.frame(regulon.current$target, regulon.current[, mode])
      score[[tf]] = pathwayscoreCoeffNorm(scale.mat,
                                          rownames(scale.mat),
                                          geneset,
                                          tf)
      Sys.sleep(1 / 100)
      counter = counter + 1
      setTxtProgressBar(pb, counter)
    }
    score.combine = do.call(cbind, score)
    score.combine = as.data.frame(t(score.combine))
  }
  else if (method == "aucell") {
    message(paste("calculating TF activity from regulon using "),
            method)
    geneSets = split(regulon$target, regulon$tf)
    message("ranking cells...")
    cells_rankings = AUCell::AUCell_buildRankings(scale.mat,
                                                  nCores = ncore, plotStats = F)
    message("calculating AUC...")
    cells_AUC = AUCell::AUCell_calcAUC(geneSets, rankings = cells_rankings,
                                       nCores = ncore)
    score.combine = data.frame(AUCell::getAUC(cells_AUC))
  }
  return(score.combine)
}

#' A function to aggregate the weighted gene expression of genes belonging to a given pathway or a regulon
#'
#' @param mat_scale A matrix of normalized gene expression with genes as the rows and cells as the columns
#' @param genenames A vector of characters corresponding to the gene names of mat_scale. It usually corresponds to the row names of the mat_matrix
#' @param pathway A data frame with the first column indicating the name of the genes in the pathway and the second column indicating the weights of genes contributing to the pathway or regulon
#' @param geneset_name String indicating the name of the pathway or the regulon
#'
#' @return A vector of inferred activity scores for every single cell
#' @export
pathwayscoreCoeffNorm = function(mat_scale, genenames, pathway, geneset_name) {
  pathway_index = match(pathway[, 1], genenames)

  if (length(pathway_index) == 1) {
    score = as.data.frame(mat_scale[pathway_index, ] * pathway[, 2], na.rm = T)
  } else {
    score = as.data.frame(colMeans(mat_scale[pathway_index, ] * pathway[, 2], na.rm = T))
  }
  colnames(score) <- geneset_name
  return(score)
}
