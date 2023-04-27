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
exp_values[sample(1:1e6, 1e6*0.25)] <- 0
#exp_values <- runif(1e6)
expMatrix <- matrix(exp_values, 1000, 1000)
rownames(expMatrix) <- paste("gene", 1:1000, sep="_")
colnames(expMatrix) <- paste("cell", 1:1000, sep="_")
expMatrix <- as(expMatrix, "dgCMatrix")

clusters <- sample(c("A","B","C"), 1000, replace = TRUE)

regulon <- data.frame(tf = sample(paste("gene", 1:50, sep = "_"), 1000, replace = TRUE),
                      idxATAC = sample(1:1000, 1000, replace = TRUE),
                      target = sample(paste("gene", 1:1000, sep = "_"), 1000, replace = TRUE))

# prepare regulon to being inputted to fast_wilcox
copy <- regulon
all.genes <- union(unique(regulon$target), unique(regulon$tf))
copy$tf <- match(copy$tf, all.genes)
copy$target <- match(copy$target, all.genes)
exprs_trans <- Matrix::t(expMatrix[all.genes,,drop=FALSE])

peak_trans <- Matrix::t(peakMatrix)
if (!is(peak_trans, "dgCMatrix")) {
  peak_trans <- as(peak_trans, "dgCMatrix")
}
if (!is(exprs_trans, "dgCMatrix")) {
  exprs_trans <- as(exprs_trans, "dgCMatrix")
}

reg.order <- order(copy$target, copy$tf, copy$idxATAC)
copy <- copy[reg.order,,drop=FALSE]

ties <- matrix(nrow=4, ncol = 1000)
total0 <- matrix(nrow=4, ncol = 1000)
total1 <- matrix(nrow=4, ncol = 1000)
U <- matrix(nrow=4, ncol = 1000)

for (cluster in c("A","B","C", "all")){
  res_row_numb <- which(c("A","B","C", "all") == cluster)
  if(cluster == "all"){
    selected_cells <- 1:nrow(exprs_trans)
  }
  else
    selected_cells <- which(clusters == cluster)

  exprs_trans_cluster <- exprs_trans[selected_cells,]
  peak_trans_cluster <- peak_trans[selected_cells,]
  for(i in seq_len(nrow(copy))){
    target_expr <- exprs_trans_cluster[,copy$target[i]]
    # identify cells with transcription factor being expressed and regulatory element being accessible
    tf_re_ind <- which(exprs_trans_cluster[,copy$tf[i]]>0 & peak_trans_cluster[,copy$idxATAC[i]]>0)
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

output <- epiregulon:::fast_wilcox(
  exprs_x = exprs_trans@x,
  exprs_i = exprs_trans@i,
  exprs_p = exprs_trans@p,
  peak_x = peak_trans@x,
  peak_i = peak_trans@i,
  peak_p = peak_trans@p,
  target_id = copy$target - 1L,
  tf_id = copy$tf - 1L,
  peak_id = copy$idxATAC - 1L,
  clusters = iclusters - 1L
)

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
