regulon <- data.frame(tf= rep(LETTERS[1:5], each = 3), target = LETTERS [6:20],
                      weight = seq(0.01,1,length.out = 15))

geneExpressionMatrix <- matrix(0, nrow = 15, ncol = 100, dimnames = list(LETTERS[6:20], NULL))
geneExpressionMatrix[1:(88*17)] <- rep(c(0, 0, 1, 2, 3, 0, 0, 1, 2, 0, 2, 0, 0, 0, 1, 3, 0), 88)

tf_tg_matrix <- matrix(0, nrow = 15, ncol = 5, dimnames = list(LETTERS[6:20], LETTERS[1:5]))

for (i in seq_len(nrow(regulon))){
  tf_tg_matrix[regulon[i,"target"], regulon[i,"tf"]] <- regulon[i, "weight"]
}

activity_matrix <- t(t(geneExpressionMatrix) %*% tf_tg_matrix)

# divide by the number of target genes
activity_matrix <- activity_matrix/3

sce<-SingleCellExperiment(assay=list(counts = geneExpressionMatrix))

test_that("calculateActivity works correctly", {
  expect_equal(calculateActivity(sce, regulon, assay = "counts"), activity_matrix)
})

# subtract mean gene expression (centering at zero)
geneExpressionMatrix <- sweep(geneExpressionMatrix, 1, rowMeans(geneExpressionMatrix),"-")
activity_matrix <- t(t(geneExpressionMatrix) %*% tf_tg_matrix)
# divide by the number of target genes
activity_matrix <- activity_matrix/3

test_that("calculateActivity normalizes correctly (zero-centering)", {
  expect_equal(calculateActivity(sce, regulon, assay = "counts", normalize = TRUE), activity_matrix)
})





