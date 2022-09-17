#' Calculate joint probability of linked TF, RE and targets
#'
#' @param expMatrix A SingleCellExperiment object or matrix containing gene expression with
#' with genes in the rows and cells in the columns
#' @param exp_assay String indicating the name of the assay in expMatrix for gene expression
#' @param exp_cutoff A scalar indicating the minimum gene expression for a gene to be considered
#' active. Default value is 1
#' @param peakMatrix A SingleCellExperiment object or matrix containing peak accessibility with
#' peaks in the rows and cells in the columns
#' @param peak_assay String indicating the name of the assay in peakMatrix for chromatin accessibility
#' @param peak_cutoff A scalar indicating the minimum peak accessibility for a peak to be
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
#' @param regulon_cutoff A scalar indicating the minimum value for the joint probability of
#' a tf-idxATAC-target trio to be retained in the pruned regulon.
#' @param clusters A vector corresponding to the cluster labels of the cells if
#' cluster-specific joint probabilities are also required. If left ```NULL```, joint probabilities
#' are calculated for all cells
#' @param aggregate A logical indicating whether to collapse the regulatory elements of the
#' same genes. If ```TRUE```, the output will only contain tf and target. If ```FALSE```, the output
#' will contain tf, idxATAC and target.
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
#' regulon = regulon, regulon_cutoff = 0.5)
#'
#' #calculate joint probability for each cluster
#' pruned.regulon <- calculateJointProbability(expMatrix = gene_sce,
#' exp_assay = "logcounts",peakMatrix = peak_sce,peak_assay = "counts",
#' regulon = regulon,clusters = gene_sce$Treatment,
#' regulon_cutoff = 0.5,aggregate = FALSE)
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
                                      regulon_cutoff = 0,
                                      clusters = NULL,
                                      aggregate = TRUE,
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
  regulon <- regulon[order(regulon$tf),]
  tf_uniq <- unique(regulon$tf)


  total_cell <- ncol(expMatrix)



  prob_matrix <- c()

  writeLines("compute joint probability for all trios")

  pb <- txtProgressBar(min = 0,
                       max = length(tf_uniq),
                       style = 3)

  counter <- 0



  prob_matrix_tf <- BiocParallel::bplapply(X = tf_uniq,
                                           FUN = calculateJointProbability_bp,
                                           regulon,
                                           expMatrix,
                                           exp_cutoff,
                                           peakMatrix,
                                           peak_cutoff,
                                           chromvarMatrix,
                                           chromvar_cutoff,
                                           clusters,
                                           uniq_clusters,
                                           BPPARAM = BPPARAM)
  prob_matrix <- do.call("rbind", prob_matrix_tf)

  # calculate the number of cells in each cluster
  clusters_freq <- table(clusters)
  clusters_freq <- c(sum(clusters_freq),clusters_freq)

  # normalize by total cell counts
  prob_matrix <- sweep(prob_matrix, MARGIN=2, STATS = clusters_freq, FUN = "/")



  #add probability matrix to original regulon
  regulon.combined <- cbind(regulon, prob_matrix)

  #if aggregate is true, collapse regulatory elements to have regulons containing tf and target
  if (aggregate == TRUE){
    regulon.combined <- stats::aggregate(prob_matrix ~ tf + target, data = regulon.combined,
                                  FUN = mean, na.rm = TRUE)
  }

  regulon.combined <- regulon.combined[regulon.combined[,"all"] > regulon_cutoff, ]
  return(regulon.combined)


}


calculateJointProbability_bp <- function(tf,
                                          regulon,
                                          expMatrix,
                                          exp_cutoff,
                                          peakMatrix,
                                          peak_cutoff,
                                          chromvarMatrix,
                                          chromvar_cutoff,
                                          clusters,
                                          uniq_clusters){
  message(tf)
  regulon_tf <- regulon[regulon$tf == tf,]

  # filter cells that did not pass tf cutoff either by expression or chromvar
  if (is.null(chromvarMatrix)){
    cells_sel_tf <- Matrix::which(expMatrix[tf,] > exp_cutoff)
  } else {
    cells_sel_tf <- Matrix::which(chromvarMatrix[tf,] > chromvar_cutoff)
  }


  #initiate a probability matrix to keep track of the number of cells that fulfill threshold for all cutoffs
  prob_matrix_tf <- matrix(0, nrow = nrow(regulon_tf), ncol = length(uniq_clusters))
  colnames(prob_matrix_tf) <- uniq_clusters

  if (length(cells_sel_tf) != 0){



    #track new clusters because cells got filtered
    new_clusters <- clusters[cells_sel_tf]


    target.short <- expMatrix[regulon_tf$target, cells_sel_tf, drop=FALSE]
    peak.short <- peakMatrix[regulon_tf$idxATAC, cells_sel_tf, drop=FALSE]

    # create target and peak matrix that pass threshold.
    # 1s represent cells that pass threshold and 0s represent cells that fail threshold

    target.bi.index <- Matrix::which(target.short > exp_cutoff, arr.ind = TRUE)
    target.bi <- Matrix::sparseMatrix(x = rep(1,nrow(target.bi.index)),
                                      i = target.bi.index[,1],
                                      j = target.bi.index[,2],
                                      dims = c(nrow(target.short), ncol(target.short)) )


    peak.bi.index <- Matrix::which(peak.short > peak_cutoff, arr.ind = TRUE)
    peak.bi <- Matrix::sparseMatrix(x = rep(1,nrow(peak.bi.index)),
                                    i = peak.bi.index[,1],
                                    j =  peak.bi.index[,2],
                                    dims = c(nrow(peak.short), ncol(peak.short)) )

    # perform element by element matrix multiplication to get cells that pass both thresholds
    target.peak.bi <- target.bi * peak.bi


    prob_matrix_tf[,"all" ] <- Matrix::rowSums(target.peak.bi)

    # also computes joint probabilities by cluster
    for (cluster in unique(clusters)){
      cluster_index = which(new_clusters  == cluster)
      if (length(cluster_index) >0) {
        prob_matrix_tf[, cluster] <- Matrix::rowSums(target.peak.bi[,cluster_index, drop=FALSE])
      }

    }
  }
  return(prob_matrix_tf)
}








