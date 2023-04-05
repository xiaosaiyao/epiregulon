# set up matrices and regulon
peakMatrix <- matrix(rbinom(1000*1000,1,0.01), 1000, 1000)
rownames(peakMatrix) <- paste("peak", 1:1000, sep="_")
colnames(peakMatrix) <- paste("cell", 1:1000, sep="_")
peakMatrix <- as(peakMatrix, "dgCMatrix")

expMatrix <- matrix(rnorm(1000*1000,0.3,1), 1000, 1000)
rownames(expMatrix) <- paste("gene", 1:1000, sep="_")
colnames(expMatrix) <- paste("cell", 1:1000, sep="_")
expMatrix <- as(expMatrix, "dgCMatrix")

clusters <- sample(c("A","B","C"), 1000, replace = TRUE)

regulon <- data.frame(tf = sample(paste("gene", 1:50, sep = "_"), 1000, replace = TRUE),
                      idxATAC = sample(1:1000, 1000, replace = TRUE),
                      target = sample(paste("gene", 1:1000, sep = "_"), 1000, replace = TRUE))


regulon <- regulon[order(regulon$tf), ]

if (is.null(clusters)) {
  cluster_id <- factor(integer(ncol(peakMatrix)))
} else {
  cluster_id <- factor(clusters)
}

######## C++
# if cluster information is provided
stats_fast <- countCells(regulon, expMatrix, peakMatrix, cluster_id, peak_cutoff=0, exp_cutoff=1, clusters)



######## old R code
peakMatrix.bi <- binarize_matrix(peakMatrix, cutoff = 0)
expMatrix.bi <- tfMatrix.bi <- binarize_matrix(expMatrix, cutoff = 1)
regulon.split <- split(regulon, regulon$tf)

chisq_bp_test <- function (n,
                           regulon.split,
                           expMatrix.bi,
                           peakMatrix.bi,
                           tfMatrix.bi,
                           clusters,
                           unique_clusters
){
  full_ncells <- ncol(peakMatrix.bi)
  has_tf <- tfMatrix.bi[regulon.split[[n]]$tf[1],] == 1
  expMatrix.bi <- expMatrix.bi[regulon.split[[n]]$target,, drop=FALSE]
  expMatrix.tf.bi <- expMatrix.bi[, has_tf, drop=FALSE]
  peakMatrix.bi <- peakMatrix.bi[regulon.split[[n]]$idxATAC, has_tf, drop=FALSE]

  triple.bi <- peakMatrix.bi * expMatrix.tf.bi
  tf_re.bi <- peakMatrix.bi

  emptyMatrix <- matrix(NA, nrow=nrow(regulon.split[[n]]), ncol=length(unique_clusters))
  colnames(emptyMatrix) <- unique_clusters
  stats <- list(triple = emptyMatrix, peak = emptyMatrix, target = emptyMatrix)

  for (selected_cluster in unique_clusters){

    if(selected_cluster != "all"){
      is_current_cluster <- as.logical(clusters == selected_cluster)
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

    stats[["triple"]][, selected_cluster] <- n_triple
    stats[["peak"]][, selected_cluster] <- tf_reCurrent
    stats[["target"]][, selected_cluster] <- expCurrent

  }


  stats

}

stats_slow <- lapply(
  X = seq_len(length(regulon.split)),
  FUN = chisq_bp_test,
  regulon.split,
  expMatrix.bi,
  peakMatrix.bi,
  tfMatrix.bi,
  cluster_id,
  c("all", levels(cluster_id))
)

triple <- list()
peak <-list()
target <-list()
for (i in seq_len(length(regulon.split))){
 triple[[i]] <-  stats_slow[[i]][["triple"]]
 peak[[i]] <-  stats_slow[[i]][["peak"]]
 target[[i]] <-  stats_slow[[i]][["target"]]
}

triple <- do.call(rbind, triple)
peak <- do.call(rbind, peak)
target <- do.call(rbind, target)

combined_stats_slow <- list(triple=triple, peak=peak, target=target)

dimnames(combined_stats_slow$triple) <- NULL
dimnames(combined_stats_slow$peak) <- NULL
dimnames(combined_stats_slow$target) <- NULL


test_that("triple", {
  expect_identical(as.matrix(combined_stats_slow$triple), as.matrix(stats_fast$triple))
})

test_that("peak", {
  expect_identical(as.matrix(combined_stats_slow$triple), as.matrix(stats_fast$triple))
})

test_that("target", {
  expect_identical(as.matrix(combined_stats_slow$triple), as.matrix(stats_fast$triple))
})


############ C++
# if cluster information is not provided
clusters <- NULL
if (is.null(clusters)) {
  cluster_id <- factor(integer(ncol(peakMatrix)))
} else {
  cluster_id <- factor(clusters)
}
stats_fast_no_cluster <- countCells(regulon, expMatrix, peakMatrix, cluster_id, peak_cutoff=0, exp_cutoff=1, clusters)
head(stats_fast_no_cluster$triple)
