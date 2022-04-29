#' A function to test for differential TF activity between pairs of single cell clusters/groups
#'
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param groups A character or integer vector of cluster or group labels for single cells
#' @param test.type String indicating the type of statistical tests to be passed to scran::findMarkers, can be "t", "wilcox". or "binom"
#' @param pval.type A string specifying how p-values are to be combined across pairwise comparisons for a given group/cluster.
#' @param direction A string specifying direction of differential TF activity, can be "up" or "down"
#' @param ... Further arguments to pass to scran::findMarkers
#'
#' @return A named list of dataframes containing differential TF activity test results for each cluster/group
#' @export
#'
#' @examples
#' set.seed(1)
#' score.combine=cbind(matrix(runif(2000,0,2), 20,100), matrix(runif(2000,0,10), 20,100))
#' rownames(score.combine)=paste0("TF",1:20)
#' colnames(score.combine)=paste0("cell",1:200)
#' cluster=c(rep(1,100),rep(2,100))
#' markers <- findDifferentialActivity(score.combine, cluster, pval.type="some", direction="up", test.type= "t")
#' sig.genes <- getSigGenes(markers, fdr_cutoff = 1, logFC_cutoff=0.1)

findDifferentialActivity <- function(activity_matrix, groups, test.type= "t", pval.type="some", direction="up", ...){

  activity_matrix=na.omit(as.matrix(activity_matrix))
  tf_markers <- scran::findMarkers(activity_matrix, groups, test.type= "t", pval.type="some", direction="up", ...)
  return(tf_markers)

}

#' A function to compile and summarize the output from findDifferentialActivity function
#'
#' @param da_list List of dataframes from running findDifferentialActivity
#' @param fdr_cutoff A numeric scalar to specify the cutoff for FDR value. Default is 0.05
#' @param logFC_cutoff A numeric scalar to specify the cutoff for log fold change.
#' @param topgenes A integer scalar to indicate the number of top ordered genes to include in output
#'
#' @return A compiled dataframe of TFs with differential activities across clusters/groups
#' @export
#'
#' @examples
#' set.seed(1)
#' score.combine=cbind(matrix(runif(2000,0,2), 20,100), matrix(runif(2000,0,10), 20,100))
#' rownames(score.combine)=paste0("TF",1:20)
#' colnames(score.combine)=paste0("cell",1:200)
#' cluster=c(rep(1,100),rep(2,100))
#' markers <- findDifferentialActivity(score.combine, cluster, pval.type="some", direction="up", test.type= "t")
#' sig.genes <- getSigGenes(markers, fdr_cutoff = 1, logFC_cutoff=0.1)
#' head(sig.genes)

getSigGenes=function(da_list, fdr_cutoff = 0.05, logFC_cutoff = NULL, topgenes = NULL){

  classes <- names(da_list)

  top.list <- lapply(seq_along(da_list), function(i){
    da_genes <- as.data.frame(da_list[[i]])
    da_genes <- da_genes[,c("p.value","FDR","summary.logFC")]

    if (is.null(logFC_cutoff)){
      logFC_cutoff = round(quantile(da_genes$summary.logFC, 0.95), digits = 1)
    }else {
      logFC_cutoff= logFC_cutoff
    }
    message ("Using a logFC cutoff of ", logFC_cutoff, " for class ", classes[i])
    da_genes <- da_genes[which(da_genes[,"FDR"] < fdr_cutoff & da_genes[, "summary.logFC"] > logFC_cutoff), ]

    if (nrow(da_genes) != 0){
      da_genes$class <- classes[[i]]
      da_genes$tf <- rownames(da_genes); rownames(da_genes) <- NULL
    }

    if (is.null(topgenes)){
      da_genes <- da_genes[with(da_genes, order(FDR,-summary.logFC)),]
    } else {
      da_genes <- da_genes[head(with(da_genes, order(FDR, -summary.logFC)), topgenes),]
    }
    #print(da_genes)


    return(da_genes)
  })

  #print(top.list)

  return(do.call(rbind, top.list))

}

regulonEnrich_ <- function(TF, regulon, corr, corr_cutoff, genesets){

  regulon.TF=unique(regulon$target[which(regulon$tf == TF & regulon[, corr] >corr_cutoff)])
  if (is.null(regulon.TF)) {
    results=data.frame(p.adjust=NA, Description=NA, GeneRatio=0, Odds.Ratio=NA)
  } else {
  enrichresults = clusterProfiler::enricher(regulon.TF, TERM2GENE = genesets)
  results=enrichresults@result
  results$GeneRatio=(as.numeric(lapply(strsplit(results$GeneRatio, split = "/"), "[",1)))/(as.numeric(lapply(strsplit(results$GeneRatio, split = "/"), "[",2)))
  results$BgRatio=(as.numeric(lapply(strsplit(results$BgRatio, split = "/"), "[",1)))/(as.numeric(lapply(strsplit(results$BgRatio, split = "/"), "[",2)))
  results$Odds.Ratio=results$GeneRatio/results$BgRatio
  results=results[order(results$p.adjust),]
  results$Description = factor(as.character(results$Description), level = unique(as.character(results$Description[nrow(results):1])))
  }
  return(results)
}


#' A function to perform geneset enrichment of user-defined regulons
#'
#' @param TF  A character vector of TF names
#' @param regulon A matrix of weighted regulon consisting of tf, targets, corr and weight
#' @param corr String indicating the column name that should be used to filter target genes for geneset enrichment. Default is "weight".
#' @param corr_cutoff A numeric scalar to indicate the cutoff to filter on the column specified by corr. Default is 0.5.
#' @param genesets A dataframe with the first column being the name of the geneset and the second column being the name of the genes
#'
#' @return A dataframe showing the significantly enriched pathways
#' @export
#'
#' @examples
#' \dontrun{
#' #retrieve genesets
#' H <- EnrichmentBrowser::getGenesets(org = "hsa", db = "msigdb", cat = "H", gene.id.type = "SYMBOL" )
#' C6 <- EnrichmentBrowser::getGenesets(org = "hsa", db = "msigdb", cat = "C6", gene.id.type = "SYMBOL" )
#'
#' #combine genesets and convert genesets to be compatible with enricher
#' gs <- c(H,C6)
#' gs.list <- do.call(rbind,lapply(names(gs), function(x) {data.frame(gs=x, genes=gs[[x]])}))
#'
#' > head(gs.list)
#' gs   genes
#' 1 M5890_HALLMARK_TNFA_SIGNALING_VIA_NFKB   ABCA1
#' 2 M5890_HALLMARK_TNFA_SIGNALING_VIA_NFKB   ACKR3
#' 3 M5890_HALLMARK_TNFA_SIGNALING_VIA_NFKB    AREG
#'
#' enrichment_results <- regulonEnrich(c("NKX2-1","GATA6","AR"), regulon=regulon.w, genesets=gs.list)
#'}
#'
regulonEnrich <- function(TF, regulon, corr="weight", corr_cutoff=0.5, genesets){
  regulonls <- lapply(TF, function(x) {
    regulonEnrich_(x, regulon, corr, corr_cutoff, genesets)
  })
  names(regulonls)=TF
  return(regulonls)
}
