regulon <- data.frame(tf = rep(LETTERS[1:10], each =2), idxATAC = 1:20, target = c(rep(letters[1:6],3), letters[4:5]))
expMatrix <- matrix(rep(c(0.2,3,6,0,0,2.3,5,0,0,0,1.1,2, 6,8,0,3),20), ncol = 20, nrow=16,
                    dimnames = list(c(LETTERS[1:10], letters[1:6]), NULL), byrow = TRUE)
regulon$idxATAC <- 1:20
peakMatrix <- matrix(rep(c(0,0,0,1), 100), ncol = 20, nrow=20,byrow=TRUE)
tfMatrix <- expMatrix[regulon$tf,]>1
targetMatrix <- expMatrix[regulon$target,]
tf_re <- tfMatrix * peakMatrix
weights <- numeric(20)
for(i in 1:20){
  weights[i] <- coin::wilcox_test(targetMatrix[i,]~factor(tf_re[i,], levels = c(1,0)))@statistic@teststatistic
}
weights <- weights/sqrt(20)
regulon.w <- regulon[,c("tf", "target")]
regulon.w$weight <- weights
sce <- SingleCellExperiment(list(logcounts = expMatrix))
test_that("addWeights works correctly using Wilcoxon test", {
  regulon <- addWeights(regulon = regulon, sce=sce, method = "wilcoxon",
                        peakMatrix = peakMatrix)
  expect_identical(regulon$weight, regulon.w[order(regulon.w$tf, regulon.w$target),]$weight)
})
