
#' Find Differential Genes
#' 
#' Functions adapted from the deprecated scran::findMarkers to identify the most differential genes through
#' pairwise t-tests with p-values combined via various methods. 
#' 
#' @param x A matrix or SummarizedExperiment object
#' @param groups A vector of length equal to the number of cells, indicating the group of each cell
#' @param direction A string specifying the direction of log-fold changes to be considered in the alternative hypothesis.
#' See `scran::pairwiseTTests`
#' @param pval.type A string specifying how p-values are to be combined across pairwise comparisons for a given group. 
#' See `scran::combineMarkers`
#' @param min.prop Numeric scalar specifying the minimum proportion of significant comparisons per gene.
#' Defaults to 0.5 when pval.type="some", otherwise defaults to zero.  See `scran::combineMarkers`
#' @param assay.type A string specifying which assay values to use. Default is `logcounts`
#' @param combined logical indicating pairwise p-values should be combined
#' 
#' @return 
#' If combined is FALSE, this function returns a list of 4 elements:
#' \itemize{
#' \item `log.p.value`: A list of data.frame, of length equal to the number of groups, indicating the 
#' log p-values of each group member vs versus all other group members
#' \item `log.FDR`: A list of data.frame, of length equal to the number of groups, indicating the 
#' log FDR of each group member vs versus all other group members
#' \item `logFC`: A list of data.frame, of length equal to the number of groups, indicating the 
#' log fold changes of each group member vs versus all other group members
#' \item `all.groups`: a vector corresponding to the unique group names
#' }
#' If combined is TRUE, output is a list of data.frame objects of length equal to the number of groups,
#' indicating the combined log p-values, log FDR values and log fold changes.
#' @export
#' @importFrom scrapper modelGeneVariances
#' @importFrom metapod combineParallelPValues
#' @examples
#' set.seed(1000)
#' sce <- scuttle::mockSCE()
#' sce <- scrapper::normalizeRnaCounts.se(sce)
#' rownames(sce) <- paste0('Gene_',1:2000)
#' simple_output <- findMarkersSimple(sce, sce$Treatment, combined=FALSE, 
#'     assay.type="logcounts", direction="any")
#' simple_output_combined <- findMarkersSimple(sce, sce$Treatment, combined=TRUE, 
#'     assay.type="logcounts", direction="any")
findMarkersSimple <- function(x, 
                              groups,
                              direction=c("any", "up", "down"),
                              pval.type=c("any", "some", "all"), 
                              min.prop=NULL,
                              assay.type="logcounts",
                              combined=TRUE) {
  
  pval.type <- match.arg(pval.type)
  direction <- match.arg(direction)
  
  if (is(x, "SummarizedExperiment")){
    mat <- assay(x, assay.type)
  } else {
    mat <- x
  }
  
  
  
  compiled <- pairwiseTTestsSimple(mat, groups, direction)
  
  if (combined) {
    combined <- combineMarkersSimple(compiled, groups, pval.type=pval.type, min.prop=min.prop)
  } else {
    compiled
  }
  
}

#' @importFrom stats pt
pairwiseTTestsSimple <- function(mat, groups, direction) {
  stats <- modelGeneVariances(mat, block=groups, block.average.policy="none", fit.trend=FALSE)$per.block
  group.sizes <- table(groups)
  
  compiled.p <- list()
  
  compiled.q <- list()
  compiled.logFC <- list()
  all.groups <- names(stats)
  for (g1 in all.groups) {
    left.mean <- stats[[g1]]$means
    left.var <- stats[[g1]]$variances
    left.n <- group.sizes[[g1]] 
    left.df <- max(0L, left.n - 1L)
    
    all.p <- list()
    for (g2 in all.groups) {
      if (g1 == g2) {
        next
      }
      
      right.mean <- stats[[g2]]$means
      right.var <- stats[[g2]]$variances
      right.n <- group.sizes[[g2]]
      right.df <- max(0L, right.n - 1L)
      
      # Perform Welch's t-test here.
      left.err <- left.var / left.n
      right.err <- right.var / right.n
      cur.err <- left.err + right.err
      cur.df <- cur.err^2 / (left.err^2 / left.df + right.err^2 / right.df)
      
      cur.lfc <- left.mean - right.mean
      cur.t <- cur.lfc / sqrt(cur.err)
      
      if (direction == "up") {
        p <- pt(cur.t, df=cur.df, lower.tail=FALSE, log.p=TRUE)
      } else if (direction == "down") {
        p <- pt(cur.t, df=cur.df, lower.tail=TRUE, log.p=TRUE)
      } else {
        p <- log(2) + pt(abs(cur.t), df=cur.df, lower.tail=FALSE, log.p=TRUE)
      }
      
      compiled.logFC[[g1]][[paste0("logFC.", g2)]] <- cur.lfc
      compiled.p[[g1]][[paste0("log.p.value.", g2)]] <- p
      compiled.q[[g1]][[paste0("log.FDR.", g2)]] <- .logBH(p)
    }
    
    
    compiled.logFC[[g1]] <- as.data.frame(compiled.logFC[[g1]])
    rownames(compiled.logFC[[g1]]) <- rownames(mat)
    compiled.p[[g1]] <- as.data.frame(compiled.p[[g1]])
    rownames(compiled.p[[g1]]) <- rownames(mat)
    compiled.q[[g1]] <- as.data.frame(compiled.q[[g1]])
    rownames(compiled.q[[g1]]) <- rownames(mat)
  }
  compiled <- list(log.p.value=compiled.p, log.FDR=compiled.q, logFC=compiled.logFC, all.groups=names(stats))
}


combineMarkersSimple <- function(compiled, groups, pval.type, min.prop) {
  
  method <- switch(pval.type, any="simes", some="holm-min", all="berger")
  if (is.null(min.prop))  {
    min.prop <- if (pval.type=="any") 0 else 0.5
  }
  
  
  all.groups <- compiled$all.groups
  combined <- list()
  
  for (g1 in all.groups){
    df <- data.frame(matrix(NA, ncol=3, nrow=nrow(compiled[[1]][[1]])))
    colnames(df) <- c("log.p.value", "log.FDR", "logFC")
    combined[[g1]] <- df
    combined[[g1]]$`log.p.value` <- combineParallelPValues(compiled[["log.p.value"]][[g1]], 
                                                       method=method, 
                                                       min.prop=min.prop, 
                                                       log.p=TRUE)$p.value  
    combined[[g1]]$`log.FDR` <- .logBH(combined[[g1]]$`log.p.value`)
    combined[[g1]]$`logFC` <- .choose_logFC(all.p=compiled[["log.p.value"]][[g1]], 
                                          all.logFC=compiled[["logFC"]][[g1]], 
                                          pval.type=pval.type, 
                                          min.prop=min.prop)
    rownames(combined[[g1]]) <- rownames(compiled[["log.p.value"]][[g1]])
     
  }
  combined
}

.logBH <- function(log.p.val) {
  o <- order(log.p.val)
  repval <- log.p.val[o] + log(length(o)/seq_along(o))
  repval <- rev(cummin(rev(repval)))
  repval[o] <- repval
  repval
}

#' @importFrom matrixStats rowRanks
.choose_logFC <- function(all.p, all.logFC, pval.type, min.prop){
  if (pval.type == "any"){
    col_index <- apply(all.p, 1, which.min)
  } else if (pval.type == "all") {
    col_index <- apply(all.p, 1, which.max)
  } else {
    p_rank <- rowRanks(as.matrix(all.p), ties.method="first")
    chosen_rank <- floor(min.prop*(ncol(all.p)-1) + 1)
    col_index <- vapply(seq_len(nrow(p_rank)), function(i){
      which(p_rank[i,] == chosen_rank)}, 1)
  }
  
  chosen_logFC <- vapply(seq_along(col_index), function(i){
    all.logFC[i, col_index[i]]}, 1)
}


