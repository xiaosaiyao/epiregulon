#' A function to test for differential TF activity between pairs of single cell clusters/groups
#'
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param groups A vector of cluster or group labels for single cells
#' @param test.type The type of statistical tests to be passed to scran::findMarkers, can be "t", "wilcox". or "binom"
#' @param pval.type A string specifying how p-values are to be combined across pairwise comparisons for a given group/cluster.
#' @param direction A string specifying direction of differential TF activitie, can be "up" or "down"
#' @param ... Further arguments to pass to scran::findMarkers
#'
#' @return A named list of dataframes containing differential TF activity test results for each cluster/group
#' @export
#'
#' @examples 1+1 = 2
findDifferentialActivity <- function(activity_matrix, groups, test.type= "t", pval.type="some", direction="up", ...){

  tf_markers <- scran::findMarkers(activity_matrix, groups, test.type= "t", pval.type="some", direction="up", ...)
  return(tf_markers)

}

#' A function to compile and summarize the output from findDifferentialActivity function
#'
#' @param da_list List of dataframes from running findDifferentialActivity
#' @param fdr_cutoff cutoff for FDR value, default is 0.05
#'
#' @return A compiled dataframe of TFs with differential activities across clusters/groups
#' @export
#'
#' @examples
getSigGenes <- function(da_list, fdr_cutoff = 0.05){

  classes <- names(da_list)

  top.list <- lapply(seq_along(da_list), function(i){
    da_genes <- as.data.frame(da_list[[i]])
    da_genes <- da_genes[,c("p.value","FDR","summary.logFC")]
    da_genes <- subset(da_genes, FDR < fdr_cutoff)
    da_genes$class <- classes[[i]]
    da_genes$tf <- rownames(da_genes); rownames(da_genes) <- NULL
    da_genes <- da_genes[order(-da_genes$summary.logFC),]
    return(da_genes)
  })

  return(do.call(rbind, top.list))

}
