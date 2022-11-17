# prepare dataset for testing Wilcoxon method
regulon <- data.frame(tf = rep(LETTERS[1:10], each =2),
                      idxATAC = 1:20,
                      target = c(rep(letters[1:6],3), letters[4:5]))
expMatrix <- matrix(rep(c(0.2,3,6,0,0,2.3,5,0,0,0,1.1,2, 6,8,0,3),20), ncol = 20, nrow=16,
                    dimnames = list(c(LETTERS[1:10], letters[1:6]), NULL), byrow = TRUE)
expMatrix <- as(expMatrix, "sparseMatrix")
regulon$idxATAC <- 1:20
regulon <- regulon[order(regulon$tf, regulon$target),]
peakMatrix <- matrix(rep(c(0,0,0,1), 100), ncol = 20, nrow=20,byrow=TRUE)
peakMatrix <- as(peakMatrix, "sparseMatrix")

# mark cells with tf expression above cutoff
tfMatrix <- expMatrix[regulon$tf,,drop = FALSE] > 1
targetMatrix <- expMatrix[regulon$target,]

# label cells with both tf being expressed and re being open
tf_re <- tfMatrix * peakMatrix
weights <- numeric(20)
for(i in 1:20){
  weights[i] <- coin::wilcox_test(targetMatrix[i,]~factor(tf_re[i,], levels = c(1,0)))@statistic@teststatistic
}


# transform z-scores to effect sizes
weights <- weights/sqrt(20)
regulon$weight <- weights
expMatrix <- SingleCellExperiment(list(logcounts = expMatrix))
test_that("addWeights works correctly using Wilcoxon test", {
  regulon.w <- addWeights(regulon = regulon,
                          expMatrix = expMatrix,
                          method = "wilcoxon",
                          peakMatrix = peakMatrix,
                          min_targets = 0)
  expect_identical(regulon.w$weight, regulon$weight)
})

# prepare expected data for logFC method
weights <- numeric(20)
for(i in 1:20){
  group1 <- Matrix::rowSums(targetMatrix[i,,drop=FALSE]*tf_re[i,,drop=FALSE])
  group0 <- Matrix::rowSums(targetMatrix[i,,drop=FALSE]*(1-tf_re[i,,drop=FALSE]))
  weights[i] <- group1/sum(tf_re[i,,drop=FALSE]) - group0/sum(1-tf_re[i,,drop=FALSE])
}
regulon$weight <- weights
test_that("addWeights works correctly using 'logFC' method", {
  regulon.w <- addWeights(regulon = regulon,
                          expMatrix = expMatrix,
                          method = "logFC",
                          peakMatrix = peakMatrix,
                          min_targets = 1)
  expect_identical(regulon.w$weight, regulon$weight)
})

# prepare dataset for testing MI method
set.seed(1010)
regulon <- data.frame(tf= rep(LETTERS[1:5], each = 3), target = LETTERS [6:20])
expMatrix <- matrix(abs(rnorm(1e4)), nrow = 20, ncol = 1e4/20, dimnames = list(LETTERS[1:20], NULL))
expMatrix[sample(seq_len(1e4), 1e4/2)] <- 0
expMatrix <- as(expMatrix, "sparseMatrix")

groupings <- rep(paste0("C", 1:20), each = 25)

averages.se <- scuttle::sumCountsAcrossCells(
  expMatrix,
  ids = groupings,
  average = TRUE,
)

expMatrix.av <- assays(averages.se)$average

regulon$weight <- NA
for(i in seq_len(nrow(regulon))){
  y2d <-  entropy::discretize2d(expMatrix.av[regulon$tf[i],],
                               expMatrix.av[regulon$target[i],],
                               numBins1 =  max(10, unique(expMatrix.av[regulon$tf[i],])),
                               numBins2 =  max(10, expMatrix.av[regulon$target[i],]))
  regulon$weight[i] <- suppressWarnings(entropy::mi.empirical(y2d))
}

regulon <- regulon[order(regulon$tf, regulon$target),]

expMatrix.sce <- SingleCellExperiment(list(logcounts = expMatrix))

test_that("addWeights works correctly using mutual information method", {
  regulon.w <- suppressWarnings(addWeights(regulon = regulon,
                                           expMatrix = expMatrix.sce,
                                           method = "MI",
                                           clusters = groupings,
                                           min_targets = 0))
  expect_identical(regulon.w$weight, regulon$weight)
})

# prepare expected data for corr method
geneExpr <- expMatrix.av[regulon$target,]
tfMatrix <- expMatrix.av[regulon$tf,]
regulon$weight <- diag(stats::cor(t(geneExpr), t(tfMatrix)))
test_that("addWeights works correctly using 'corr' method", {
  regulon.w <- addWeights(regulon = regulon, expMatrix = expMatrix.sce, method = "corr",
                          clusters = groupings, min_targets = 0)
  expect_identical(regulon.w$weight, regulon$weight)
})



