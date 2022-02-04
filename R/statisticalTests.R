findDifferentialActivity <- function(activity_matrix, groups, test.type= "t", pval.type="some", direction="up", ...){

  tf_markers <- scran::findMarkers(activity_matrix, groups, test.type= "t", pval.type="some", direction="up", ...)
  return(tf_markers)

}

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
