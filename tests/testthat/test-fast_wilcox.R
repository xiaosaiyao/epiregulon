calculate_U <- function(x,y){
  U <- 0
  for(i in seq_along(x)){
    U <- U + sum(x[i]>y) + 0.5*sum(x[i]==y)
  }
  U
}

# set up matrices and regulon
set.seed(1001)
peakMatrix <- matrix(rbinom(1000*1000,1,0.01), 1000, 1000)
rownames(peakMatrix) <- paste("peak", 1:1000, sep="_")
colnames(peakMatrix) <- paste("cell", 1:1000, sep="_")
peakMatrix <- as(peakMatrix, "dgCMatrix")

unique_exp_values <- runif(5000)*3
# make ties
exp_values <- sample(unique_exp_values,1e6, replace = TRUE)
#exp_values <- runif(1e6)
exp_values[sample(1:1e6, 1e6*0.25)] <- 0

expMatrix <- matrix(exp_values, 1000, 1000)
rownames(expMatrix) <- paste("gene", 1:1000, sep="_")
colnames(expMatrix) <- paste("cell", 1:1000, sep="_")
expMatrix <- as(expMatrix, "dgCMatrix")

#clusters <- NULL
clusters <- sample(LETTERS[1:3], 1000, replace = TRUE)

regulon <- data.frame(tf = sample(paste("gene", 1:50, sep = "_"), 1000, replace = TRUE),
                      idxATAC = sample(1:1000, 1000, replace = TRUE),
                      target = sample(paste("gene", 1:1000, sep = "_"), 1000, replace = TRUE))

# prepare regulon to being inputted to fast_wilcox
peakMatrix <- binarize_matrix(peakMatrix, cutoff = NULL)
copy <- regulon
all.targets <- sort(unique(regulon$target))
all.tfs <- sort(unique(regulon$tf))
copy$tf <- match(copy$tf, all.tfs)
copy$target <- match(copy$target, all.targets)
if (!is.null(clusters)){
  expMatrix_tfs_clusters <- expMatrix[all.tfs,,drop = FALSE]
  for(cluster in unique(clusters)){
    cluster_ind <- which(clusters == cluster)
    expMatrix_tfs_clusters[,cluster_ind, drop = FALSE] <- binarize_matrix(expMatrix_tfs_clusters[,cluster_ind, drop = FALSE],
                                                                          cutoff = NULL)
  }
}
expMatrix_tfs <- binarize_matrix(expMatrix[all.tfs,,drop = FALSE], cutoff = NULL)
exprs_trans_target <- Matrix::t(expMatrix[all.targets,,drop=FALSE])
exprs_trans_tf <- Matrix::t(expMatrix_tfs)
all.peaks <- sort(unique(copy$idxATAC))
copy$idxATAC <- match(copy$idxATAC, all.peaks)
peak_trans <- Matrix::t(peakMatrix[all.peaks,,drop=FALSE])
if (!is(peak_trans, "dgCMatrix")) {
  peak_trans <- as(peak_trans, "dgCMatrix")
}

reg.order <- order(copy$target, copy$tf, copy$idxATAC)
copy <- copy[reg.order,,drop=FALSE]

ties <- matrix(nrow=4, ncol = 1000)
total0 <- matrix(nrow=4, ncol = 1000)
total1 <- matrix(nrow=4, ncol = 1000)
U <- matrix(nrow=4, ncol = 1000)

for (cluster in c("all","A","B","C")){
  res_row_numb <- which(c("all","A","B","C") == cluster)
  if(cluster == "all"){
    selected_cells <- 1:nrow(exprs_trans_target)
  }
  else
    selected_cells <- which(clusters == cluster)
  exprs_trans_cluster <- exprs_trans_target[selected_cells,]
  peak_trans_cluster <- peak_trans[selected_cells,]
  expMatrix_tfs_cluster <- binarize_matrix(expMatrix[all.tfs,selected_cells,drop = FALSE], cutoff = NULL)
  exprs_trans_tf_cluster <- Matrix::t(expMatrix_tfs_cluster)
  for(i in seq_len(nrow(copy))){
    target_expr <- exprs_trans_cluster[,copy$target[i]]
    tf_cutoff <- mean(exprs_trans_tf_cluster[,copy$tf[i]])
    peak_cutoff <- mean(peak_trans_cluster[,copy$idxATAC[i]])
    # identify cells with transcription factor being expressed and regulatory element being accessible
    tf_re_ind <- which(exprs_trans_tf_cluster[,copy$tf[i]]>tf_cutoff & peak_trans_cluster[,copy$idxATAC[i]]>peak_cutoff)
    pos_group_expr <- target_expr[tf_re_ind]
    # identify cells without transcription factor being expressed or regulatory element being accessible
    neg_ind <- setdiff(1:nrow(exprs_trans_cluster), tf_re_ind)
    neg_group_expr <- target_expr[neg_ind]
    # index of the first elements in the arrays of tied cells
    tie_pos_first <- which(diff(sort(target_expr))==0)
    tie_pos_last <- rev(length(target_expr)- which(diff(sort(target_expr,
                                                             decreasing = TRUE))==0) +1)
    mid_elements_ind <- intersect(tie_pos_first, tie_pos_last)
    tie_pos_first <- tie_pos_first[!tie_pos_first %in% mid_elements_ind]
    tie_pos_last <- tie_pos_last[!tie_pos_last %in% mid_elements_ind]
    # determine length of the tied elements arrays by subtracting first position from last position + 1
    ties[res_row_numb,i] <- as.numeric(sum((tie_pos_last-tie_pos_first+1)^3) - sum(tie_pos_last-tie_pos_first+1))
    total0[res_row_numb,i] <- as.numeric(length(neg_ind))
    total1[res_row_numb,i] <- as.numeric(length(tf_re_ind))
    U[res_row_numb,i] <- calculate_U(pos_group_expr, neg_group_expr)
  }
}

fclusters <- factor(clusters)
iclusters <- as.integer(fclusters)


output <- fast_wilcox(
  exprs_x = exprs_trans_target@x,
  exprs_i = exprs_trans_target@i,
  exprs_p = exprs_trans_target@p,
  exprs_tf_x = as.logical(exprs_trans_tf@x),
  exprs_tf_i = exprs_trans_tf@i,
  exprs_tf_p = exprs_trans_tf@p,
  peak_x = peak_trans@x,
  peak_i = peak_trans@i,
  peak_p = peak_trans@p,
  target_id = copy$target - 1L,
  tf_id = copy$tf - 1L,
  peak_id = copy$idxATAC - 1L,
  clusters = numeric(0),
  cell_numb = nrow(exprs_trans_target)
)

if(!is.null(clusters)){
  exprs_trans_tf_clusters <- Matrix::t(expMatrix_tfs_clusters)
  fclusters <- factor(clusters)
  fclusters_order <- order(levels(fclusters))
  iclusters <- as.integer(fclusters)
  output_clusters <- fast_wilcox(
    exprs_x = exprs_trans_target@x,
    exprs_i = exprs_trans_target@i,
    exprs_p = exprs_trans_target@p,
    exprs_tf_x = as.logical(exprs_trans_tf_clusters@x),
    exprs_tf_i = exprs_trans_tf_clusters@i,
    exprs_tf_p = exprs_trans_tf_clusters@p,
    peak_x = peak_trans@x,
    peak_i = peak_trans@i,
    peak_p = peak_trans@p,
    target_id = copy$target - 1L,
    tf_id = copy$tf - 1L,
    peak_id = copy$idxATAC - 1L,
    clusters = iclusters - 1L,
    cell_numb = nrow(exprs_trans_target)
  )
  output <- mapply(function(x,y) rbind(x,y), output, output_clusters, SIMPLIFY = FALSE)
}

test_that("fast_wilcox function calculates ties correctly", {
  expect_identical(output$ties, ties)
})

test_that("fast_wilcox function calculates group sizes correctly", {
  expect_identical(output$total0, total0)
  expect_identical(output$total1, total1)
})

test_that("fast_wilcox function calculates U statistic correctly", {
  expect_identical(output$auc, U)
})
