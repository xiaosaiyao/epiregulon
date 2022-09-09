#' Title
#'
#' @param expr
#' @param atac
#' @param regulon
#' @param clusters
#' @param expr_cutoff
#' @param atac_cutoff
#'
#' @return
#' @export
#'
#' @examples
#' create a mock singleCellExperiment object for gene expression matrix
#' set.seed(1000)
#' gene_sce <- scuttle::mockSCE()
#' gene_sce <- scuttle::logNormCounts(gene_sce)
#' gene_gr <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr3","chr4"), nrow(gene_sce)/4),
#'                    ranges = IRanges(start = seq(from = 1, length.out=nrow(gene_sce), by = 1000),
#'                    width = 100))
#' rownames(gene_sce) <- paste0("Gene_", 1:ncol(gene_sce))
#' gene_gr$name <- rownames(gene_sce)
#' rowRanges(gene_sce) <- gene_gr
#'
#' # create a mock singleCellExperiment object for peak matrix
#' peak_gr <- GRanges(seqnames = "chr1",
#'                    ranges = IRanges(start = seq(from = 1, to = 10000, by = 1000), width = 100))
#' peak_counts <- matrix(sample(x = 0:4, size = ncol(gene_sce)*length(peak_gr), replace = TRUE),
#'                       nrow = length(peak_gr), ncol=ncol(gene_sce))
#' peak_sce <- SingleCellExperiment(list(counts = peak_counts), colData = colData(gene_sce))
#' rowRanges(peak_sce) <- peak_gr
#' rownames(peak_sce) <- paste0("peak",1:10)
#' #'
#'
#' # create a mock regulon
#' regulon <- data.frame(tf = c(rep("Gene_1",5), rep("Gene_2",10)),
#'                       target = c(paste0("Gene_000",2:6), paste0("Gene_00",11:20)))
#'

calculateJointProbability <- function(expr,
                                 expr_assay = "logcounts",
                                 atac,
                                 atac_assay = "PeakMatrix",
                                 regulon,
                                 regulon_cutoff = 0,
                                 clusters = NULL,
                                 cutoff_cluster = "all",
                                 expr_cutoff = 1,
                                 atac_cutoff = 0,
                                 aggregate = TRUE) {

  if (class(expr) %in% c("SingleCellExperiment", "SummarizedExperiment", "RangedSummarizedExperiment")){
    expr <- assay(expr, expr_assay)
  }

  if (class(atac) %in% c("SingleCellExperiment", "SummarizedExperiment", "RangedSummarizedExperiment")){
    atac <- assay(atac, atac_assay)
  }

  expr_index <- apply(expr, 1, function(x) which(x > expr_cutoff))
  names(expr_index) <- rownames(expr)
  atac_index <- apply(atac, 1, function(x) which(x > atac_cutoff))
  names(atac_index) <- rownames(atac)

  if (is.null(clusters)) {
    uniq_clusters <- "all"
  } else {
    uniq_clusters <- c("all", unique(clusters))
  }

  #add all clusters
  #uniq_clusters

  #initiate conditional probability matrix
  prob_matrix <- matrix(NA, nrow = nrow(regulon), ncol = length(uniq_clusters))
  colnames(prob_matrix) = uniq_clusters

  for (cluster in uniq_clusters){

    if (cluster == "all"){
      cluster_index <- 1:ncol(expr)
      } else {
      cluster_index <- which(clusters == cluster)
      }

    writeLines(paste("cluster", cluster))

    pb <- txtProgressBar(min = 0,
                         max = nrow(regulon),
                         style = 3)

    cluster_length <- length(cluster_index)

    #compute conditional probability for all trios
    for (i in 1:nrow(regulon)){

      # count the number of cells satisfying all 3 conditions
      all_count <- length (Reduce(intersect, list(tf = expr_index[[regulon$tf[i]]],
                             target = expr_index[[regulon$target[i]]],
                             atac = atac_index[[regulon$idxATAC[i]]],
                             cluster = cluster_index)))

      prob_matrix[i, cluster] <- all_count/cluster_length


      Sys.sleep(1 / 100)

      setTxtProgressBar(pb, i)
    }

  }
  regulon.combined <- cbind(regulon, prob_matrix)
  regulon.pruned <- regulon.combined[which(regulon.combined[,cutoff_cluster] > regulon_cutoff),]

  if (aggregate == TRUE){
    regulon.pruned <- aggregate(all ~ tf + target, data = regulon.pruned,
                             FUN = mean, na.rm = TRUE)
  }
return(regulon.pruned)

}
