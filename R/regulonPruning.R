#' Prune regulons for true transcription factor - regulatory elements - target genes relationships
#'
#' @param regulon A dataframe informing the gene regulatory relationship with the ```tf``` column
#' representing transcription factors, ```idxATAC``` corresponding to the index in the peakMatrix and
#'  ```target``` column representing target genes
#' @param expMatrix A SingleCellExperiment object or matrix containing gene expression with
#' genes in the rows and cells in the columns
#' @param peakMatrix A SingleCellExperiment object or matrix containing peak accessibility with
#' peaks in the rows and cells in the columns
#' @param chromvarMatrix A SingleCellExperiment object or matrix containing averaged accessibility at the TF
#' binding sites with tfs in the rows and cells in the columns. This can be used as an alternative to TF expression
#' @param exp_assay String indicating the name of the assay in expMatrix for gene expression
#' @param peak_assay String indicating the name of the assay in peakMatrix for chromatin accessibility
#' @param chromvar_assay String indicating the name of the assay in chromvarMatrix for chromatin accessibility
#' @param test String indicating whether `binom` or `chi.sq` test should be performed
#' @param clusters A vector corresponding to the cluster labels of the cells if
#' cluster-specific joint probabilities are also required. If left ```NULL```, joint probabilities
#' are calculated for all cells
#' @param exp_cutoff A scalar indicating the minimum gene expression above which
#' gene is considered active. Default value is 1. Applied to both transcription
#' factors and target genes.
#' @param peak_cutoff A scalar indicating the minimum peak accessibility above which peak is
#' considered open. Default value is 0
#' @param chromvar_cutoff A scalar indicating the minimum chromvar values for a tf to be
#' considered active. Default value is 0
#' @param regulon_cutoff A scalar indicating the maximal value for p-value for a tf-idxATAC-target trio
#' to be retained in the pruned regulon.
#' @param p_adj A logical indicating whether p adjustment should be performed
#' @param prune_value String indicating whether to filter regulon based on `pval` or `padj`.
#' @param aggregate A logical indicating whether to collapse the regulatory elements of the
#' same genes and use them as a whole. Note that checking with the `peak_cutoff`
#' is made before the collapse.
#' @param BPPARAM A BiocParallelParam object specifying whether calculation should be parallelized.
#' Default is set to BiocParallel::MulticoreParam()
#
#' @return A dataframe of a pruned regulon with p-values indicating the probability of independence
#' either for all cells or for individual clusters
#'
#' @details
#' The function prunes the network by performing tests of independence on the observed number of cells
#' jointly expressing transcription factor (TF), regulatory element (RE) and target gene (TG) vs
#' the expected number of cells if TF/RE and TG are independently expressed.
#'
#' In other words, if no regulatory relationship exists, the expected probability of cells expressing all
#' three elements is P(TF, RE) * P(TG), that is, the product of (1) proportion of cells both expressing transcription factor
#' and having accessible corresponding regulatory element, and (2) proportion of cells expressing
#' target gene. The expected number of cells expressing all three elements is therefore n*P(TF, RE)*P(TG),
#' where n is the total number of cells. However, if a TF-RE-TG relationship exists,
#' we expect the observed number of cells jointly having all three elements (TF, RE, TG) to deviate from
#' the expected number of cells predicted from an independent relationship.
#'
#' If the user provides cluster assignment, the tests of independence are performed on a per-cluster basis
#' in addition to providing all cells statistics. This enables pruning by cluster, and thus yields cluster-specific
#' gene regulatory relationships.
#'
#' We implement two tests, the binomial test and the chi-square test.
#'
#' In the binomial test, the expected probability is P(TF, RE) * P(TG), and the number of trials is the number of cells,
#' and the observed successes is the number of cells jointly expressing all three elements.
#'
#' In the chi-square test, the expected probability for having all 3 elements active is also P(TF, RE) * P(TG) and
#' the probability otherwise is 1- P(TF, RE) * P(TG). The observed cell count for the active category is the number of cells
#' jointly expressing all three elements, and the cell count for the inactive category is n - n_triple.
#'
#'
#'
#' @export
#' @import utils SingleCellExperiment
#'
#' @examples
#' # create a mock singleCellExperiment object for gene expMatrixession matrix
#' set.seed(1000)
#' gene_sce <- scuttle::mockSCE()
#' gene_sce <- scuttle::logNormCounts(gene_sce)
#' rownames(gene_sce) <- paste0("Gene_",1:2000)
#'
#' # create a mock singleCellExperiment object for peak matrix
#' peak_gr <- GRanges(seqnames = "chr1",
#'                   ranges = IRanges(start = seq(from = 1, to = 10000, by = 100), width = 100))
#' peak_counts <- matrix(sample(x = 0:4, size = ncol(gene_sce)*length(peak_gr), replace = TRUE),
#'                      nrow = length(peak_gr), ncol=ncol(gene_sce))
#' peak_sce <- SingleCellExperiment(list(counts = peak_counts), colData = colData(gene_sce))
#' rownames(peak_sce) <- paste0("Peak_",1:100)
#'
#' # create a mock regulon
#' regulon <- data.frame(tf = c(rep("Gene_1",10), rep("Gene_2",10)),
#'                      idxATAC = sample(1:100, 20),
#'                      target = c(paste0("Gene_", sample(3:2000,10)),
#'                                 paste0("Gene_",sample(3:2000,10))))
#'
#  # prune regulon
#' pruned.regulon <- pruneRegulon(expMatrix = gene_sce,
#' exp_assay = "logcounts", peakMatrix = peak_sce, peak_assay = "counts",
#' regulon = regulon, clusters = gene_sce$Treatment, aggregate = FALSE, regulon_cutoff = 0.5)

#'
#' @author Xiaosai Yao, Tomasz Wlodarczyk


pruneRegulon <- function(regulon,
                         expMatrix = NULL,
                         peakMatrix = NULL,
                         chromvarMatrix = NULL,
                         exp_assay = "logcounts",
                         peak_assay = "PeakMatrix",
                         chromvar_assay = NULL,
                         test = c("chi.sq","binom"),
                         clusters = NULL,
                         exp_cutoff = 1,
                         peak_cutoff = 0,
                         chromvar_cutoff = 0,
                         regulon_cutoff = 0.05,
                         p_adj = TRUE,
                         prune_value = "pval",
                         aggregate = FALSE,
                         BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)){

  # choose test method
  test <- match.arg(test)
  message("pruning network with ", test, " tests using a regulon cutoff of ", prune_value, "<", regulon_cutoff)


  # extracting assays from SE
  if (checkmate::test_class(expMatrix,classes = "SummarizedExperiment")){
    expMatrix <- assay(expMatrix, exp_assay)
  }

  if (checkmate::test_class(peakMatrix,classes = "SummarizedExperiment")){
    peakMatrix <- assay(peakMatrix, peak_assay)
  }

  if (!is.null (chromvar_assay)){
    if (checkmate::test_class(chromvarMatrix,classes = "SummarizedExperiment")){
      chromvarMatrix <- assay(chromvarMatrix, chromvar_assay)
    }
  }



  # clean up regulons by removing tf and targets not found in regulons
  regulon <- regulon[regulon$tf %in% rownames(expMatrix),]
  regulon <- regulon[regulon$target %in% rownames(expMatrix),]
  regulon <- regulon[order(regulon$tf),]


  # binarize peak and expression matrices according to cutoff
  message("binarizing matrices")
  peakMatrix.bi <- binarize_matrix(peakMatrix, peak_cutoff)
  expMatrix.bi <- tfMatrix.bi <- binarize_matrix(expMatrix, exp_cutoff)


  unique_clusters <- c("all", as.character(sort(unique(clusters))))

  res <- list()
  regulon.split <- split(regulon, regulon$tf)

  message("pruning regulons")
  if (test == "binom") {
    # Perform binomial test
    res <- BiocParallel::bplapply(
      X = seq_len(length(regulon.split)),
      FUN = binom_bp,
      regulon.split,
      expMatrix.bi,
      peakMatrix.bi,
      tfMatrix.bi,
      clusters,
      unique_clusters,
      BPPARAM = BPPARAM
      )



  } else if (test == "chi.sq") {
    # Perform chi-square test
    res <- BiocParallel::bplapply(
      X = seq_len(length(regulon.split)),
      FUN = chisq_bp,
      regulon.split,
      expMatrix.bi,
      peakMatrix.bi,
      tfMatrix.bi,
      clusters,
      unique_clusters,
      BPPARAM = BPPARAM
      )

  } else {

    stop("test must be either binom or chi.sq")
  }


  res <- do.call("rbind", res)


  # append test stats to regulon
  regulon.combined  <- cbind(regulon, res[,grep("pval_",colnames(res)), drop = FALSE], res[,grep("stats_",colnames(res)), drop = FALSE])


  # add p-value adjustment

  if (p_adj){
    pval_columns <- grepl("pval_", colnames(regulon.combined))
    qvalue <- apply(regulon.combined[,pval_columns, drop = FALSE], 2,
                    function(x) {stats::p.adjust(x, method = "holm", n = nrow(regulon.combined))})
    qval_columns <- gsub("pval_", "padj_", colnames(regulon.combined)[pval_columns])
    colnames(qvalue) <- qval_columns
    regulon.combined <- cbind(regulon.combined, qvalue)
  }

  # prune by p-value
  regulon.prune_value <- regulon.combined[,grepl(prune_value, colnames(regulon.combined)), drop = FALSE]
  prune_value_min <- apply(regulon.prune_value, 1, function (x){
    if (sum(is.na(x)) == length(x))
     Inf
    else
      min(x, na.rm = TRUE)
  })
  regulon.combined <- regulon.combined[which(prune_value_min < regulon_cutoff),]

  function(x) {if (length(x)>0) min(x) else Inf}

  # if aggregate is true, collapse regulatory elements to have regulons containing tf and target
  if (aggregate == TRUE){
    message("aggregating regulons...")
    aggregate_by <- colnames(regulon.combined)[grep("stats|pval|padj",
                                                    colnames(regulon.combined))]

    regulon.combined <- stats::aggregate(regulon.combined[aggregate_by],
                                         by = regulon.combined[c("tf", "target")],
                                         FUN = mean, na.rm = TRUE)


  }


  return(regulon.combined)

}




binarize_matrix <- function(matrix_obj, cutoff){
  if (is(matrix_obj, "dgCMatrix")) {
    matrix_obj@x <- as.double(matrix_obj@x > cutoff)
    matrix_obj
  } else {
    matrix_obj.bi.index <- Matrix::which(matrix_obj > cutoff, arr.ind = TRUE)
    Matrix::sparseMatrix(x = rep(1, nrow(matrix_obj.bi.index)),
                         i = matrix_obj.bi.index[,1],
                         j =  matrix_obj.bi.index[,2],
                         dims = dim(matrix_obj),
                         dimnames = dimnames(matrix_obj))
  }
}


binom_bp <- function(n,
                     regulon.split,
                     expMatrix.bi,
                     peakMatrix.bi,
                     tfMatrix.bi,
                     clusters,
                     unique_clusters,
                     BPPARAM = BPPARAM
                     ){

  full_ncells <- ncol(peakMatrix.bi)
  has_tf <- tfMatrix.bi[regulon.split[[n]]$tf[1],] == 1
  expMatrix.bi <- expMatrix.bi[regulon.split[[n]]$target,, drop=FALSE]
  expMatrix.tf.bi <- expMatrix.bi[, has_tf, drop=FALSE]
  peakMatrix.bi <- peakMatrix.bi[regulon.split[[n]]$idxATAC, has_tf, drop=FALSE]

  triple.bi <- peakMatrix.bi * expMatrix.tf.bi
  tf_re.bi <- peakMatrix.bi

  res <- list()

  for (selected_cluster in unique_clusters){
    if(selected_cluster != "all"){
      is_current_cluster <- clusters == selected_cluster
      expCurrent <- as.vector(expMatrix.bi %*% is_current_cluster)
      n_cells <- sum(is_current_cluster)

      #subset is_current_cluster to only cells with tf greater than cutoff
      is_current_cluster <- is_current_cluster[has_tf]
      tf_reCurrent <- as.vector(tf_re.bi %*% is_current_cluster)
      n_triple <- as.vector(triple.bi %*% is_current_cluster)

    } else {
      expCurrent <- Matrix::rowSums(expMatrix.bi)
      n_triple <- Matrix::rowSums(triple.bi)
      tf_reCurrent <- Matrix::rowSums(tf_re.bi)
      n_cells <- full_ncells
    }

    null_probability <- tf_reCurrent * expCurrent / n_cells ^ 2
    p.value <- binomTest(n_triple, n_cells, p=null_probability)
    z_score <- stats::qnorm(p.value/2) * sign(null_probability - n_triple/n_cells)

    res[[selected_cluster]] <- cbind(p.value, z_score)
    colnames(res[[selected_cluster]]) <- c(paste0("pval_",selected_cluster), paste0("stats_", selected_cluster))
  }

  res <- do.call("cbind", res)


}

binomTest <- function (k, size, p) {
  if (size >= 10000) {
    e1 <- p * size
    e2 <- size - e1
    chi <- (k - e1)^2 / e1 + (size - k - e2)^2 / e2
    return(stats::pchisq(chi, df = 1, lower.tail = FALSE))
  }

  p.value <- rep_len(1, length(k))
  for (ip in unique(p)) {
      current <- p == ip
      d <- stats::dbinom(0:size, prob = ip, size = size)
      o <- order(d)
      cumsump <- cumsum(d[o])[order(o)]
      p.value[current] <- cumsump[k[current] + 1]
  }
  p.value
}



chisq_bp <- function(n,
                     regulon.split,
                     expMatrix.bi,
                     peakMatrix.bi,
                     tfMatrix.bi,
                     clusters,
                     unique_clusters,
                     BPPARAM = BPPARAM
                     ){

  full_ncells <- ncol(peakMatrix.bi)
  has_tf <- tfMatrix.bi[regulon.split[[n]]$tf[1],] == 1
  expMatrix.bi <- expMatrix.bi[regulon.split[[n]]$target,, drop=FALSE]
  expMatrix.tf.bi <- expMatrix.bi[, has_tf, drop=FALSE]
  peakMatrix.bi <- peakMatrix.bi[regulon.split[[n]]$idxATAC, has_tf, drop=FALSE]

  triple.bi <- peakMatrix.bi * expMatrix.tf.bi
  tf_re.bi <- peakMatrix.bi

  res <- list()

  for (selected_cluster in unique_clusters){
    if(selected_cluster != "all"){
      is_current_cluster <- clusters == selected_cluster
      expCurrent <- as.vector(expMatrix.bi %*% is_current_cluster)
      n_cells <- sum(is_current_cluster)

      #subset is_current_cluster to only cells with tf greater than cutoff
      is_current_cluster <- is_current_cluster[has_tf]
      tf_reCurrent <- as.vector(tf_re.bi %*% is_current_cluster)
      n_triple <- as.vector(triple.bi %*% is_current_cluster)

    } else {
      expCurrent <- Matrix::rowSums(expMatrix.bi)
      n_triple <- Matrix::rowSums(triple.bi)
      tf_reCurrent <- Matrix::rowSums(tf_re.bi)
      n_cells <- full_ncells
    }

    null_probability <- tf_reCurrent * expCurrent / n_cells ^ 2
    res[[selected_cluster]] <- chisqTest(k = n_triple, size = n_cells, p = null_probability)
    colnames(res[[selected_cluster]]) <- c(paste0("pval_",selected_cluster), paste0("stats_", selected_cluster))
  }

  res <- do.call("cbind", res)


}

chisqTest <- function (k, size, p) {

  e1 <- p * size
  e2 <- size - e1
  chi <- (k - e1)^2 / e1 + (size - k - e2)^2 / e2
  df <- cbind(p = stats::pchisq(chi, df = 1, lower.tail = FALSE), stats = chi)
  return(df)


}



