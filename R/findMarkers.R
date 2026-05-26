#' @importFrom scrapper modelGeneVariances
#' @importFrom metapod combineParallelPValues
findMarkersSimple <- function(x, 
                              groups,
                              direction=c("any", "up", "down"),
                              pval.type=c("any", "some", "all"), 
                              min.prop=NULL,
                              assay.type="logcounts",
                              combined=TRUE) {
  
  pval.type <- match.arg(pval.type)
  direction <- match.arg(direction)
  
  mat <- assay(x, assay.type)
  
  
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
        p <- pt(cur.t, df = cur.df, lower.tail = FALSE, log.p = TRUE)
      } else if (direction == "down") {
        p <- pt(cur.t, df = cur.df, lower.tail = TRUE, log.p = TRUE)
      } else {
        p <- log(2) + pt(abs(cur.t), df = cur.df, lower.tail = FALSE, log.p = TRUE)
      }
      
      compiled.logFC[[g1]][[paste0("logFC.", g2)]] <- cur.lfc
      compiled.p[[g1]][[paste0("p.value.", g2)]] <- p
      compiled.q[[g1]][[paste0("FDR.", g2)]] <- .logBH(p)
    }
    
    
    compiled.logFC[[g1]] <- as.data.frame(compiled.logFC[[g1]])
    rownames(compiled.logFC[[g1]]) <- rownames(mat)
    compiled.p[[g1]] <- as.data.frame(compiled.p[[g1]])
    rownames(compiled.p[[g1]]) <- rownames(mat)
    compiled.q[[g1]] <- as.data.frame(compiled.q[[g1]])
    rownames(compiled.q[[g1]]) <- rownames(mat)
  }
  compiled <- list(p.value=compiled.p, FDR=compiled.q, logFC=compiled.logFC, all.groups=names(stats))
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
    colnames(df) <- c("p.value", "FDR", "logFC")
    combined[[g1]] <- df
    combined[[g1]]$`p.value` <- combineParallelPValues(compiled[["p.value"]][[g1]], 
                                                       method = method, 
                                                       min.prop = min.prop, 
                                                       log.p = TRUE)$p.value  
    combined[[g1]]$`FDR` <- .logBH(combined[[g1]]$`p.value`)
    combined[[g1]]$logFC <- .choose_logFC(all.p = compiled[["p.value"]][[g1]], 
                                          all.logFC = compiled[["logFC"]][[g1]], 
                                          pval.type = pval.type, 
                                          min.prop = min.prop)
     
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
    p_rank <- rowRanks(as.matrix(all.p), ties.method = "first")
    chosen_rank <- floor(min.prop*(ncol(all.p)-1) + 1)
    col_index <- vapply(seq_len(nrow(p_rank)), function(i){
      which(p_rank[i,] == chosen_rank)}, 1)
  }
  
  chosen_logFC <- vapply(seq_along(col_index), function(i){
    all.logFC[i, col_index[i]]}, 1)
}


