#' Test for differential TF activity between pairs of single cell clusters/groups
#'
#' @param activity_matrix A matrix of TF activities inferred from calculateActivity
#' @param groups `r lifecycle::badge("deprecated")` A character or integer vector of cluster or group labels for single cells
#' @param clusters A character or integer vector of cluster or group labels for single cells
#' @param test.type String indicating the type of statistical tests to be passed to scran::findMarkers, can be "t", "wilcox". or "binom"
#' @param pval.type A string specifying how p-values are to be combined across pairwise comparisons for a given group/cluster.
#' @param direction A string specifying direction of differential TF activity, can be "up" or "down"
#' @param ... Further arguments to pass to scran::findMarkers
#'
#' @return A named list of dataframes containing differential TF activity test results for each cluster/group
#'
#' @importFrom lifecycle deprecated
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' score.combine <- cbind(matrix(runif(2000,0,2), 20,100), matrix(runif(2000,0,10), 20,100))
#' rownames(score.combine) <- paste0("TF",1:20)
#' colnames(score.combine) <- paste0("cell",1:200)
#' cluster <- c(rep(1,100),rep(2,100))
#' markers <- findDifferentialActivity(
#' activity_matrix = score.combine,
#' clusters = cluster,
#' pval.type = "some",
#' direction = "up",
#' test.type = "t")
#' sig.genes <- getSigGenes(markers, fdr_cutoff = 1, logFC_cutoff = 0.1)
#' @author Xiaosai Yao, Shang-yang Chen

findDifferentialActivity <- function(activity_matrix,
                                     clusters,
                                     test.type = "t",
                                     pval.type = "some",
                                     direction="up",
                                     groups=deprecated(),
                                     ...){

  if(lifecycle::is_present(groups)){
      lifecycle::deprecate_warn( "1.0.0", "findDifferentialActivity(groups)",
                                 "findDifferentialActivity(clusters)")
    clusters <- groups
  }

  activity_matrix <- stats::na.omit(as.matrix(activity_matrix))
  tf_markers <- scran::findMarkers(activity_matrix, clusters, test.type=test.type,
                                   pval.type=pval.type, direction=direction, ...)
  return(tf_markers)

}

#' Compile and summarize the output from findDifferentialActivity function
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
#' score.combine <- cbind(matrix(runif(2000,0,2), 20,100), matrix(runif(2000,0,10), 20,100))
#' rownames(score.combine) <- paste0("TF",1:20)
#' colnames(score.combine) <- paste0("cell",1:200)
#' cluster <- c(rep(1,100),rep(2,100))
#' markers <- findDifferentialActivity(score.combine, cluster, pval.type = "some", direction = "up",
#' test.type = "t")
#' sig.genes <- getSigGenes(markers, fdr_cutoff = 1, logFC_cutoff = 0.1)
#' utils::head(sig.genes)
#' @author Xiaosai Yao, Shang-yang Chen

getSigGenes <- function(da_list,
                        fdr_cutoff = 0.05,
                        logFC_cutoff = NULL,
                        topgenes = NULL){

  classes <- names(da_list)

  top.list <- lapply(seq_along(da_list), function(i){
    da_genes <- as.data.frame(da_list[[i]])
    da_genes <- da_genes[,c("p.value","FDR","summary.logFC")]

    if (is.null(logFC_cutoff)){
      logFC_cutoff <- round(stats::quantile(da_genes$summary.logFC, 0.95, na.rm=TRUE), digits = 1)
    }else {
      logFC_cutoff <- logFC_cutoff
    }
    message ("Using a logFC cutoff of ", logFC_cutoff, " for class ", classes[i])
    da_genes <- da_genes[which(da_genes[,"FDR"] < fdr_cutoff & da_genes[, 3] > logFC_cutoff), ]

    if (nrow(da_genes) != 0){
      da_genes$class <- classes[[i]]
      da_genes$tf <- rownames(da_genes); rownames(da_genes) <- NULL
    }

    if (is.null(topgenes)){

      da_genes <- da_genes[order(da_genes$FDR, -(da_genes[, 3])),]
    } else {
      da_genes <- da_genes[head(order(da_genes$FDR, -(da_genes[, 3])), topgenes),]
    }
    #print(da_genes)


    return(da_genes)
  })

  #print(top.list)

  return(do.call(rbind, top.list))

}

regulonEnrich_ <- function(TF,
                           regulon,
                           corr,
                           corr_cutoff,
                           genesets){
  message(TF)
  regulon.TF <- unique(regulon$target[which(regulon$tf == TF & regulon[, corr] > corr_cutoff)])
  if (length(regulon.TF) < 3) {
    results <- data.frame(p.adjust = NA, Description = NA, GeneRatio = 0, Odds.Ratio = NA)
  } else {
  enrichresults <- clusterProfiler::enricher(regulon.TF, TERM2GENE = genesets)

  results <- enrichresults@result
  results$GeneRatio <- (as.numeric(lapply(strsplit(results$GeneRatio, split = "/"), "[",1)))/
    (as.numeric(lapply(strsplit(results$GeneRatio, split = "/"), "[",2)))
  results$BgRatio <- (as.numeric(lapply(strsplit(results$BgRatio, split = "/"), "[",1)))/
    (as.numeric(lapply(strsplit(results$BgRatio, split = "/"), "[",2)))
  results$Odds.Ratio <- results$GeneRatio/results$BgRatio
  results <- results[order(results$p.adjust),]
  results$Description <- factor(as.character(results$Description), levels = unique(as.character(results$Description[nrow(results):1])))
  }
  return(results)
}


#' Perform geneset enrichment of user-defined regulons
#'
#' @param TF  A character vector of TF names
#' @param regulon A matrix of weighted regulon consisting of tf, targets, corr and weight
#' @param weight String indicating the column name that should be used to filter target genes for geneset enrichment. Default is "weight".
#' @param weight_cutoff A numeric scalar to indicate the cutoff to filter on the column specified by corr. Default is 0.5.
#' @param genesets A dataframe with the first column being the name of the geneset and the second column being the name of the genes
#'
#' @return A dataframe showing the significantly enriched pathways
#' @export
#'
#' @examples
#' #retrieve genesets
#' H <- EnrichmentBrowser::getGenesets(org = "hsa", db = "msigdb",
#' cat = "H", gene.id.type = "SYMBOL" )
#' C6 <- EnrichmentBrowser::getGenesets(org = "hsa", db = "msigdb",
#' cat = "C6", gene.id.type = "SYMBOL" )
#'
#' #combine genesets and convert genesets to be compatible with enricher
#' gs <- c(H,C6)
#' gs.list <- do.call(rbind,lapply(names(gs), function(x) {
#' data.frame(gs=x, genes=gs[[x]])}))
#'
#' head(gs.list)
#'
#' #get regulon
#' library(dorothea)
#' data(dorothea_hs, package = "dorothea")
#' regulon <- dorothea_hs
#' enrichment_results <- regulonEnrich(c("ESR1","AR"), regulon = regulon, weight = "mor",
#' genesets = gs.list)
#'
#' @author Xiaosai Yao

regulonEnrich <- function(TF,
                          regulon,
                          weight = "weight",
                          weight_cutoff = 0.5,
                          genesets) {
    regulonls <- lapply(TF, function(x) {
      regulonEnrich_(x, regulon, weight, weight_cutoff, genesets)
    })
    names(regulonls) <- TF
    return(regulonls)
  }
