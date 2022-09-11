#' A function to calculate joint probability of linked TF, RE and targets
#'
#' @param expMatrix A SingleCellExperiment or matrix containing gene expression with
#' with genes in the rows and cells in the columns
#' @param exp_assay String indicating the name of the assay in expMatrix for gene expression
#' @param exp_cutoff A scalar indicating the minimum gene expression for a gene to be considered
#' active. Default value is 1
#' @param peakMatrix A SingleCellExperiment or matrix containing peak accessibility with
#' peaks in the rows and cells in the columns
#' @param peak_assay String indicating the name of the assay in peakMatrix for chromatin accessibility
#' @param peak_cutoff A scalar indicating the minimum peak accessibility for a peak to be
#' considered open. Default value is 0
#' @param regulon A dataframe informing the gene regulatory relationship with the "tf" column
#' representing transcription factors, "idxATAC"
#' corresponding to the index in the peakMatrix and "target" column
#' representing target genes
#' @param regulon_cutoff A scalar indicating the minimum value for the joint probability of
#' a tf-idxATAC-target trio to be retained in the pruned regulon.
#' @param clusters A vector corresponding to the cluster labels of the cells if
#' cluster-specific joint probabilities are also required. If left null, joint probabilities
#' are calculated for all cells
#' @param aggregate A logical indicating whether to collapse the regulatory elements of the
#' same genes. If TRUE, the output will only contain tf and target. If FALSE, the output
#' will contain tf, idxATAC and target.
#'
#' @return A dataframe of a pruned regulon containing joint probabilities for tf-idxATAC-target trios
#' either for all cells or in a cluster-specific manner
#' @export
#' @import utils SingleCellExperiment stats
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
                                      regulon,
                                      regulon_cutoff = 0,
                                      clusters = NULL,
                                      aggregate = TRUE) {

  if (checkmate::test_class(expMatrix,classes = "SummarizedExperiment")){
    expMatrix <- assay(expMatrix, exp_assay)
  }

  if (checkmate::test_class(peakMatrix,classes = "SummarizedExperiment")){
    peakMatrix <- assay(peakMatrix, peak_assay)
  }


  total_cell <- ncol(expMatrix)

  exp_index <- apply(expMatrix, 1, function(x) which(x > exp_cutoff))
  names(exp_index) <- rownames(expMatrix)
  peak_index <- apply(peakMatrix, 1, function(x) which(x > peak_cutoff))
  names(peak_index) <- rownames(peakMatrix)

  if (is.null(clusters)) {
    uniq_clusters <- "all"
  } else {
    uniq_clusters <- c("all", unique(clusters))
  }

  #clean up regulon
  regulon <- regulon[regulon$tf %in% rownames(expMatrix),]
  regulon <- regulon[regulon$target %in% rownames(expMatrix),]

  #initiate conditional probability matrix
  prob_matrix <- matrix(NA, nrow = nrow(regulon), ncol = length(uniq_clusters))
  colnames(prob_matrix) = uniq_clusters

  #initiate count vector
  all_count <- list()


  writeLines("compute joint probability for all trios")
  pb <- txtProgressBar(min = 0,
                       max = nrow(regulon),
                       style = 3)

  for (i in 1:nrow(regulon)){

    # count the number of cells satisfying all 3 conditions
    all_count[[i]] <- Reduce(intersect, list(tf = exp_index[[regulon$tf[i]]],
                                                target = exp_index[[regulon$target[i]]],
                                                peakMatrix = peak_index[[regulon$idxATAC[i]]]))

    prob_matrix[i, "all"] <- length(all_count[[i]])/total_cell


    Sys.sleep(1 / 100)

    setTxtProgressBar(pb, i)
  }



  # Prune network by removing trios with probability = 0
  wanted_trios <- which(prob_matrix[,"all"] > regulon_cutoff)
  prob_matrix <- prob_matrix[wanted_trios,]
  regulon <- regulon[wanted_trios,]
  all_count <- all_count[wanted_trios]

  writeLines("compute cluster-specific joint probability")

  for (cluster in unique(clusters)){
    cluster_index <- which(clusters == cluster)
    cluster_length <- length(cluster_index)

    cluster_count <- sapply(1:length(all_count), function (x) length(intersect(all_count[[x]], cluster_index)))

    prob_matrix[, cluster] <- cluster_count/cluster_length


  }

  regulon.combined <- cbind(regulon, prob_matrix)

  if (aggregate == TRUE){
    regulon.combined <- aggregate(prob_matrix ~ tf + target, data = regulon.combined,
                             FUN = mean, na.rm = TRUE)
  }
 return(regulon.combined)

}
