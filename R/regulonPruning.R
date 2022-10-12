#' Calculate joint probability of linked TF, RE and targets
#'
#' @param expMatrix A SingleCellExperiment object or matrix containing gene expression with
#' with genes in the rows and cells in the columns
#' @param exp_assay String indicating the name of the assay in expMatrix for gene expression
#' @param exp_cutoff A scalar indicating the minimum gene expression above which gene is considered
#' active. Default value is 1
#' @param peakMatrix A SingleCellExperiment object or matrix containing peak accessibility with
#' peaks in the rows and cells in the columns
#' @param peak_assay String indicating the name of the assay in peakMatrix for chromatin accessibility
#' @param peak_cutoff A scalar indicating the minimum peak accessibility above which peak is
#' considered open. Default value is 0
#' @param chromvarMatrix A SingleCellExperiment object or matrix containing averaged accessibility at the TF
#' binding sites with tfs in the rows and cells in the columns. This can be used as an alternative to TF expression
#' @param chromvar_assay String indicating the name of the assay in chromvarMatrix for chromatin accessibility
#' @param chromvar_cutoff A scalar indicating the minimum chromvar values for a tf to be
#' considered active. Default value is 0
#' @param regulon A dataframe informing the gene regulatory relationship with the ```tf``` column
#' representing transcription factors, ```idxATAC```
#' corresponding to the index in the peakMatrix and ```target``` column
#' representing target genes
#' @param regulon_prop_cutoff A scalar indicating the minimum value for the joint probability of
#' a tf-idxATAC-target trio to be retained in the pruned regulon.
#' @param regulon_p_cutoff A scalar indicating the maximal value for p-value for a tf-idxATAC-target trio
#' to be retained in the pruned regulon.
#' @param clusters A vector corresponding to the cluster labels of the cells if
#' cluster-specific joint probabilities are also required. If left ```NULL```, joint probabilities
#' are calculated for all cells
#' @param aggregate A logical indicating whether to collapse the regulatory elements of the
#' same genes. If ```TRUE```, the output will only contain tf and target. If ```FALSE```, the output
#' will contain tf, idxATAC and target.
#' @param aggregate_by A string indicating the name of the columns to aggregate targets by
#' @param triple_prop A logical indicating whether number of cells with identified tf-re-tg triple
#' should be included in output
#' @param p_val_corr logical indicating whether p value Bonferroni correction for multiple comparison
#' has to be done
#' @param BPPARAM A BiocParallelParam object specifying whether calculation should be parallelized.
#' Default is set to BiocParallel::MulticoreParam()
#'
#' @return A dataframe of a pruned regulon containing joint probabilities for tf-idxATAC-target trios
#' either for all cells or for individual clusters
#' @details This function calculates the joint probability for each of the TF-peak-target trios to be
#' active - that is, out of all the cells, how many cells have the TF and target expression exceed
#' ```exp_cutoff``` and chromatin accessibility exceed ```peak_cutoff``` simultaneously.
#' The joint probability can be used to prune the networks since a true regulatory relationship
#' likely requires cells to express the transcription factor, have an accessible peak region and
#' expressing the target gene simultaneously. While there could be time delays between tf binding,
#' chromatin accessibility and target gene expression, requiring baseline expression of all three
#' components greatly enhances the likelihood that this regulatory relationship holds true.
#'
#' This function can also compute cluster-specific joint probabilities. The output can be filtered to
#' generate cluster-specific networks which can be fed into differential network analysis (to be continued).
#'
#' The aggregate function outputs either a bipartite network of the form TF-target (```aggregate = TRUE```)
#' or a tripartite network of the form TF-RE-target (```aggregate = FALSE```).
#
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
#' #' create a mock singleCellExperiment object for peak matrix
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
#' # calculate joint probability for all cells
#' pruned.regulon <- calculateJointProbability(expMatrix = gene_sce,
#' exp_assay = "logcounts", peakMatrix = peak_sce,peak_assay = "counts",
#' regulon = regulon, regulon_prop_cutof = 0.5)
#'
#' #calculate joint probability for each cluster
#' pruned.regulon <- calculateJointProbability(expMatrix = gene_sce,
#' exp_assay = "logcounts",peakMatrix = peak_sce,peak_assay = "counts",
#' regulon = regulon,clusters = gene_sce$Treatment,
#' regulon_prop_cutof = 0.5,aggregate = FALSE)
#'
#' @author Xiaosai Yao



calculateJointProbability <- function(expMatrix,
                                      exp_assay = "logcounts",
                                      exp_cutoff = 1,
                                      peakMatrix,
                                      peak_assay = "PeakMatrix",
                                      peak_cutoff = 0,
                                      chromvarMatrix = NULL,
                                      chromvar_assay = NULL,
                                      chromvar_cutoff = 0,
                                      regulon,
                                      regulon_prop_cutoff = 0,
                                      regulon_p_cutoff = 1,
                                      clusters = NULL,
                                      aggregate = TRUE,
                                      aggregate_by = "triple_prop_all",
                                      triple_prop = TRUE,
                                      p_val_corr = FALSE,
                                      BPPARAM=BiocParallel::MulticoreParam()
                                      ) {

  # #convert delayedMatrix to dgCMatrix
  # if (checkmate::test_class(expMatrix,classes = "DelayedMatrix")) {
  #   writeLines("converting DelayedMatrix to dgCMatrix")
  #   expMatrix <- as(expMatrix, Class = "dgCMatrix")
  # }
  # if (checkmate::test_class(peakMatrix,classes = "DelayedMatrix")){
  #   writeLines("converting DelayedMatrix to dgCMatrix")
  #   peakMatrix <- as(peakMatrix, Class = "dgCMatrix")
  # }


  if (checkmate::test_class(expMatrix,classes = "SummarizedExperiment")){
    expMatrix <- assay(expMatrix, exp_assay)
  }

  if (checkmate::test_class(peakMatrix,classes = "SummarizedExperiment")){
    peakMatrix <- assay(peakMatrix, peak_assay)
  }

  if (!is.null (chromvar_assay)){
    if (checkmate::test_class(chromvarMatrix,classes = "SummarizedExperiment")){
      chromvarMatrix <- assay(chromvar_assay, peak_assay)
    }
  }



  if (is.null(clusters)) {
    uniq_clusters <- "all"
  } else {
    uniq_clusters <- c("all", sort(unique(clusters)))
  }

  #clean up regulon
  regulon <- regulon[regulon$tf %in% rownames(expMatrix),]
  regulon <- regulon[regulon$target %in% rownames(expMatrix),]

  #order by TF
  regulon <- split(regulon, regulon$tf)

  total_cell <- ncol(expMatrix)

  writeLines("compute joint probability for all trios")


  prob_matrix_tf <- BiocParallel::bplapply(X = regulon,
                                           FUN = calculateJointProbability_bp,
                                           expMatrix,
                                           exp_cutoff,
                                           peakMatrix,
                                           peak_cutoff,
                                           chromvarMatrix,
                                           chromvar_cutoff,
                                           clusters,
                                           uniq_clusters,
                                           triple_prop,
                                           total_cell = total_cell,
                                           BPPARAM = BPPARAM)
  prob_matrix <- do.call("rbind", prob_matrix_tf)
  regulon <- do.call(rbind, regulon)

  #add probability matrix to original regulon
  regulon.combined <- cbind(regulon, prob_matrix)
  if (p_val_corr){
    p_val_columns <- grepl("p_val", colnames(regulon.combined))
    q_value <- apply(regulon.combined[,p_val_columns], 2,
                     function(x) {stats::p.adjust(x, method = "holm", n = nrow(regulon))})
    q_val_columns <- gsub("p_val", "padj_val", colnames(regulon.combined)[p_val_columns])
    colnames(q_value) <- q_val_columns
    regulon.combined <- cbind(regulon.combined, q_value)
  }

  # pruning
  # by triple_prop_all
  regulon.combined <- regulon.combined[which(regulon.combined["triple_prop_all"] > regulon_prop_cutoff), ]

  # by p-value
  p_value <- regulon.combined[,grepl("p_val", colnames(regulon.combined)), drop = FALSE]
  p_value_min <- apply(p_value, 1, min)
  regulon.combined <- regulon.combined[which(p_value_min < regulon_p_cutoff),]

  # aggregation
  # if aggregate is true, collapse regulatory elements to have regulons containing tf and target
  if (aggregate == TRUE){

    aggr_formula <- eval(parse(text = paste0(aggregate_by, "~ tf + target")))
    regulon.combined <- stats::aggregate(aggr_formula, data = regulon.combined,
                                         FUN = mean, na.rm = TRUE)
  }

  return(regulon.combined)


}


calculateJointProbability_bp <- function(regulon,
                                          expMatrix,
                                          exp_cutoff,
                                          peakMatrix,
                                          peak_cutoff,
                                          chromvarMatrix,
                                          chromvar_cutoff,
                                          clusters,
                                          uniq_clusters,
                                          total_cell,
                                          triple_prop){
  message(regulon$tf[1])


  # filter cells that did not pass tf cutoff either by expression or chromvar
  if (is.null(chromvarMatrix)){
    tf.bi.index <- Matrix::which(expMatrix[regulon$tf[1],,drop = FALSE] > exp_cutoff, arr.ind = TRUE)
  } else {
    tf.bi.index <- Matrix::which(chromvarMatrix[regulon$tf[1],,drop = FALSE] > chromvar_cutoff, arr.ind = TRUE)
  }
  tf.bi <- Matrix::sparseMatrix(x = rep(1,nrow(tf.bi.index)),
                                i = tf.bi.index[,1],
                                j = tf.bi.index[,2],
                                dims = c(1, total_cell))



  # create target and peak matrices
  target.exp <- expMatrix[regulon$target,,drop = FALSE]
  re.peak <- peakMatrix[regulon$idxATAC,,drop = FALSE]
  # 1s represent cells that pass threshold and 0s represent cells that fail threshold


  target.bi.index <- Matrix::which(target.exp > exp_cutoff, arr.ind = TRUE)
  target.bi <- Matrix::sparseMatrix(x = rep(1,nrow(target.bi.index)),
                                      i = target.bi.index[,1],
                                      j = target.bi.index[,2],
                                      dims = c(nrow(target.exp), ncol(target.exp)) )


  peak.bi.index <- Matrix::which(re.peak > peak_cutoff, arr.ind = TRUE)
  peak.bi <- Matrix::sparseMatrix(x = rep(1,nrow(peak.bi.index)),
                                    i = peak.bi.index[,1],
                                    j =  peak.bi.index[,2],
                                    dims = c(nrow(re.peak), ncol(re.peak)))

  # identify cells with tf being expressed and chromatin of corresponding re accessible
  tf_re.bi <- Matrix::t(Matrix::t(peak.bi)*as.vector(tf.bi))


  # identify cells with tf-re-tg triples
  triple.bi <- tf_re.bi * target.bi

  res_list <- lapply(as.list(uniq_clusters), function(x) test_triple(x,
                                                                     clusters,
                                                                     tf_re.bi,
                                                                     target.bi,
                                                                     triple.bi,
                                                                     total_cell,
                                                                     triple_prop))
  res_matrix <- BiocGenerics::Reduce(cbind, res_list)
  if(triple_prop)
    colnames(res_matrix) <- paste0(rep(c("p_val_", "z_score_", "triple_prop_"), length(uniq_clusters)),
                                   rep(uniq_clusters, each =3))

  else
    colnames(res_matrix) <- paste0(c("p_val_", "z_score_"), rep(uniq_clusters, each = 2))

  return(res_matrix)
}


test_triple <- function(selected_cluster, clusters, tf_re.bi, target.bi, triple.bi, n_cells, triple_prop){
  if (selected_cluster == "all"){
    n_tf_re <- Matrix::rowSums(tf_re.bi)
    n_target <- Matrix::rowSums(target.bi)
    n_triple <- Matrix::rowSums(triple.bi)
  }
  else{
    n_tf_re <- Matrix::rowSums(tf_re.bi[,clusters==selected_cluster, drop = FALSE])
    n_target <- Matrix::rowSums(target.bi[,clusters==selected_cluster, drop = FALSE])
    n_triple <- Matrix::rowSums(triple.bi[,clusters==selected_cluster, drop = FALSE])
    n_cells <- sum(clusters==selected_cluster)
  }
  res <- t(mapply(binom_test, n_triple, n_cells, n_tf_re, n_target))
  if (triple_prop)
    return(cbind(res, n_triple/n_cells))
  res
}

binom_test <- function(n_triple, n_cells, n_tf_re, n_target){
  null_probability <- n_tf_re*n_target/n_cells^2
  binom_res <- stats::binom.test(n_triple, n_cells, null_probability)
  z_score <- stats::qnorm(binom_res$p.value/2)*sign(null_probability - binom_res$estimate)
  c(binom_res$p.value, z_score)
}


