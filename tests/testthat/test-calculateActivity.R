regulon <- data.frame(tf = rep(LETTERS[1:5], each = 3),
                      target = LETTERS [6:20],
                      weight = seq(0.01,1,length.out = 15))

geneExpressionMatrix <- matrix(0, nrow = 15, ncol = 100, dimnames = list(LETTERS[6:20], NULL))
geneExpressionMatrix[1:(88*17)] <- rep(c(0, 0, 1, 2, 3, 0, 0, 1, 2, 0, 2, 0, 0, 0, 1, 3, 0), 88)
geneExpressionMatrix <- as(geneExpressionMatrix, "sparseMatrix")

tf_tg_matrix <- matrix(0, nrow = 15, ncol = 5, dimnames = list(LETTERS[6:20], LETTERS[1:5]))
tf_tg_matrix <- as(tf_tg_matrix, "sparseMatrix")

for (i in seq_len(nrow(regulon))){
  tf_tg_matrix[regulon[i,"target"], regulon[i,"tf"]] <- regulon[i, "weight"]
}

activity_matrix <- as.matrix(Matrix::t(Matrix::t(geneExpressionMatrix) %*%
                                         tf_tg_matrix))

# divide by the number of target genes
activity_matrix <- activity_matrix/3

sce <- SingleCellExperiment(assay = list(counts = geneExpressionMatrix))

test_that("calculateActivity works correctly", {
  expect_equal(as.matrix(calculateActivity(sce,
                                 regulon = regulon,
                                 exp_assay = "counts")),
               activity_matrix)
})


clusters <- rep(c("C1", "C2"), each = 50)
regulon$weight_C1 <- seq(0.9,0.1,length.out = 15)
regulon$weight_C2 <- seq(0.2,0.7,length.out = 15)
geneExpr_C1 <- geneExpressionMatrix[,1:50]
geneExpr_C2 <- geneExpressionMatrix[,51:100]

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

activity_matrix <- cbind(Matrix::t(Matrix::t(geneExpr_C1) %*% tf_tg_matrix_C1),
                         Matrix::t(Matrix::t(geneExpr_C2) %*% tf_tg_matrix_C2))

# divide by the number of target genes
activity_matrix <- as.matrix(activity_matrix/3)

regulon2 <- DataFrame(regulon[,c("tf","target")])
regulon2$weight <- regulon[,c("weight_C1","weight_C2")]
colnames(regulon2$weight) <- c("C1","C2")

test_that("calculateActivity works correctly with clusters", {
  expect_equal(as.matrix(calculateActivity(sce,
                                 regulon = regulon2,
                                 exp_assay = "counts",
                                 clusters = clusters,
                                 FUN = colSums)),
               activity_matrix)
})

# subtract mean gene expression (centering at zero)
geneExpressionMatrix <- sweep(geneExpressionMatrix, 1,
                              Matrix::rowMeans(geneExpressionMatrix),"-")
activity_matrix <- as.matrix(Matrix::t(Matrix::t(geneExpressionMatrix)
                                       %*% tf_tg_matrix))

# divide by the number of target genes
activity_matrix <- activity_matrix/3

test_that("calculateActivity normalizes correctly (zero-centering)", {
  expect_equal(as.matrix(calculateActivity(sce,
                                 regulon = regulon,
                                 exp_assay = "counts",
                                 normalize = TRUE)),
               activity_matrix)
})





