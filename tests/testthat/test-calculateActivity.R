regulon <- S4Vectors::DataFrame(tf = rep(LETTERS[1:5], times = 1:5),
                      target = LETTERS [6:20],
                      weight = seq(0,1,length.out = 15))

# account for scenario when weight==NA
regulon$weight[2] <- NA

geneExpressionMatrix <- matrix(0, nrow = 15, ncol = 100, dimnames = list(LETTERS[6:20], NULL))
geneExpressionMatrix[1:(88*17)] <- rep(c(0, 0, 1, 2, 3, 0, 0, 1, 2, 0, 2, 0, 0, 0, 1, 3, 0), 88)
geneExpressionMatrix <- as(geneExpressionMatrix, "sparseMatrix")


### test createTfTgMat
tf_tg_matrix <- matrix(0, nrow = 15, ncol = 5, dimnames = list(LETTERS[6:20], LETTERS[1:5]))
tf_tg_matrix <- as(tf_tg_matrix, "sparseMatrix")

for (i in seq_len(nrow(regulon))){
  tf_tg_matrix[regulon[i,"target"], regulon[i,"tf"]] <- regulon[i, "weight"]
}

# replace NA with 0 in tf_tg_matrix
tf_tg_matrix[is.na(tf_tg_matrix)] <- 0

tf_tg_matrix2 <- createTfTgMat(regulon, "weight", clusters=NULL)

test_that("createTfTgMat works correctly", {
  expect_equal(as.matrix(tf_tg_matrix2), as.matrix(tf_tg_matrix))
})


### test calculateScore
activity_matrix <- as.matrix(Matrix::t(Matrix::t(geneExpressionMatrix) %*%
                                         tf_tg_matrix))
activity_matrix2 <- calculateScore(geneExpressionMatrix, tf_tg_matrix2)

test_that("calculateScore works correctly", {
  expect_equal(as.matrix(Matrix::t(activity_matrix2)), activity_matrix)
})

### test calculateFrequency
# divide by the number of target genes
freq <- c(1,3,4,5)
names(freq) <- c("B","C","D", "E") #no weight != 0 in A, 1 NA weight changed to 0 in B
freq2 <- calculateFrequency(regulon=regulon, mode="weight")
freq3 <- as(freq2, "double")
names(freq3) <- names(freq2)

test_that("calculateFrequency works correctly", {
  expect_equal(freq3, freq)
})

### test normalizeByFrequency
activity_matrix.norm <- activity_matrix/c(1,1,3,4,5)
activity_matrix.norm2 <- normalizeByFrequency(activity_matrix2, freq2)
test_that("normalizeByFrequency works correctly", {
  expect_equal( as.matrix(Matrix::t(activity_matrix.norm2)), activity_matrix.norm)
})

### test calculateActivity
sce <- SingleCellExperiment(assay = list(counts = geneExpressionMatrix))
activity_matrix.3 <- calculateActivity(sce, regulon = regulon, exp_assay = "counts")

test_that("calculateActivity works correctly", {
  expect_equal(as.matrix(activity_matrix.3),
               activity_matrix.norm[Matrix:::rowSums(activity_matrix.norm)!=0,]
               )
})

### subtract mean gene expression (centering at zero)
geneExpressionMatrix.mean <- sweep(geneExpressionMatrix, 1,
                              Matrix::rowMeans(geneExpressionMatrix),"-")
activity_matrix.center <- as.matrix(Matrix::t(Matrix::t(geneExpressionMatrix.mean)
                                              %*% tf_tg_matrix))

activity_matrix.center.norm <- activity_matrix.center/c(1,1,3,4,5)

test_that("calculateActivity normalizes correctly (zero-centering)", {
  expect_equal(as.matrix(calculateActivity(sce,
                                           regulon = regulon,
                                           exp_assay = "counts",
                                           normalize = TRUE)),
               activity_matrix.center.norm[Matrix:::rowSums(activity_matrix.center.norm)!=0,])
})

###################################################################################
# test with cluster information
clusters <- rep(c("C1", "C2"), each = 50)
regulon$weight_C1 <- seq(0.9,0,length.out = 15)
regulon$weight_C2 <- seq(0,0.7,length.out = 15)

regulon$weight <- cbind(regulon$weight, regulon$weight_C1, regulon$weight_C2)
colnames(regulon$weight) <- c("all", "C1", "C2")

geneExpr_C1 <- geneExpressionMatrix[,1:50]
geneExpr_C2 <- geneExpressionMatrix[,51:100]


### test createTfTgMat

tf_tg_matrix_C1 <- matrix(0, nrow = 15, ncol = 5,
                          dimnames = list(LETTERS[6:20], LETTERS[1:5]))

for (i in seq_len(nrow(regulon))){
  tf_tg_matrix_C1[regulon[i,"target"], regulon[i,"tf"]] <- regulon[i, "weight_C1"]
}

tf_tg_matrix_C2 <- matrix(0, nrow = 15, ncol = 5,
                          dimnames = list(LETTERS[6:20], LETTERS[1:5]))

for (i in seq_len(nrow(regulon))){
  tf_tg_matrix_C2[regulon[i,"target"], regulon[i,"tf"]] <- regulon[i, "weight_C2"]
}

tf_tg_matrix3 <- list(C1=tf_tg_matrix_C1, C2=tf_tg_matrix_C2)

tf_tg_matrix4 <- createTfTgMat(regulon, "weight", clusters=clusters)

tf_tg_matrix4 <- lapply(tf_tg_matrix4, as.matrix)
test_that("createTfTgMat works correctly with clusters", {
  expect_equal(tf_tg_matrix4, tf_tg_matrix3)
})


### test calculateScore

activity_matrix3 <- as(cbind(Matrix::t(Matrix::t(geneExpr_C1) %*% tf_tg_matrix_C1),
                         Matrix::t(Matrix::t(geneExpr_C2) %*% tf_tg_matrix_C2)), "CsparseMatrix")

activity_matrix4 <- matrix(0, nrow=ncol(geneExpressionMatrix), ncol=length(unique(regulon$tf)))
rownames(activity_matrix4) <- colnames(geneExpressionMatrix)
colnames(activity_matrix4) <- colnames(tf_tg_matrix4[[1]])
activity_matrix4 <- as(activity_matrix4, "CsparseMatrix")

activity_matrix4 <- calculateScore(geneExpressionMatrix, tf_tg_matrix4, clusters=clusters, activity_matrix4)

test_that("calculateScore works correctly with clusters", {
  expect_equal(Matrix::t(activity_matrix4), activity_matrix3)
})


### test calculateFrequency
freq3 <- freq4 <- initiateMatCluster(clusters, nrow = length(unique(regulon$tf)), 1)
rownames(freq3) <- rownames(freq4) <- unique(regulon$tf)
freq3[,"C1"] <- c(1,2,3,4,4)
freq3[,"C2"] <- 1:5
freq3[,"all"] <- c(1,1,3,4,5)

freq6 <- calculateFrequency(freq4, regulon, mode="weight")

test_that("calculateFrequency works correctly with clusters", {
  expect_equal(freq6, freq3)
})

### test normalizeByFrequency
activity_matrix.norm3 <- activity_matrix3
activity_matrix.norm3[,1:50] <- activity_matrix3[,1:50]/freq3[,"C1"]
activity_matrix.norm3[,51:100] <- activity_matrix3[,51:100]/freq3[,"C2"]
activity_matrix.norm4 <- normalizeByFrequency(activity_matrix4, freq6, clusters)

test_that("normalizeByFrequency works correctly with clusters", {
  expect_equal( as.matrix(Matrix::t(activity_matrix.norm4)), as.matrix(activity_matrix.norm3))
})



### test calculateActivity
activity_matrix.norm5 <- calculateActivity(sce,
                                           regulon = regulon,
                                           exp_assay = "counts",
                                           clusters = clusters,
                                           FUN = "mean")
test_that("calculateActivity works correctly with clusters", {
  expect_equal(activity_matrix.norm5,
               activity_matrix.norm3[Matrix:::rowSums(activity_matrix.norm3)!=0,])
})






