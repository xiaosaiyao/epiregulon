geneM <- matrix(0, nrow = 4, ncol = 100, dimnames = list(LETTERS[1:4], NULL))
clusters <- rep(c("C1", "C2"), each = 50)
geneM[c("A", "D"),1:50] <- rep(c(1.1,2,3,4,5),10)
geneM[c("A", "C"), 51:100] <- rep(c(1.1,0,0,2,0), 10)
geneM["C", 1:50] <- rep(c(2,3,0,0,0),10)
geneM["B",] <- rep(c(0,1.1,0,0,0,2,3,0,0,0),10)
geneM["A",][c(1:30)*3] <- 0

peakM <- matrix(0, nrow = 4, ncol =100, dimnames = list(LETTERS[1:4], NULL))
peakM[c(1,2),] <- rep(c(0,1), 50)
peakM[3,] <- rep(c(1,0,0,0),25)
peakM[2, 1:50] <- rep(c(0,1,1,0,1,1,0,0,0,0),5)
peakM[2, 51:100] <- rep(c(0,4,0,0,0,0,0,0,4,0),5)
peakM[4, 1:50] <- rep(c(1,1,1,1,0), 10)
peakM[4, 51:100] <- rep(c(0,0,1,1,0),10)

regulon <- data.frame(tf = c("A","A","A","C"), target = c("B", "B", "C", "B"), idxATAC = c(1:4))

peakM.b <- peakM > 0
geneM.b <- geneM > 1
tf_re.b <- geneM.b[regulon$tf,]*peakM.b[regulon$idxATAC,]
triplets.b <- tf_re.b*geneM.b[regulon$target,]
target.b <- geneM.b[regulon$target,]
triplet_n <- rowSums(triplets.b)
tf_re_prob <- rowSums(tf_re.b)/100
target_prob <- rowSums(target.b)/100
pvals <- c()
for(i in 1:4){
  pvals[i]<-binom.test(triplet_n[i], 100, target_prob[i]*tf_re_prob[i])$p.val
}


test_that("pruneRegulon calculates p-values in binomial test correctly", {
  regulon.p <- pruneRegulon(regulon = regulon, expMatrix = geneM, peakMatrix = peakM, regulon_cutoff = 2)
  expect_identical(regulon.p$pval_all, pvals)
})


peakM.b <- peakM.b[,1:50]
geneM.b <- geneM.b[,1:50]
tf_re.b <- geneM.b[regulon$tf,]*peakM.b[regulon$idxATAC,]
triplets.b <- tf_re.b*geneM.b[regulon$target,]
target.b <- geneM.b[regulon$target,]
triplet_n <- rowSums(triplets.b)
tf_re_prob <- rowSums(tf_re.b)/50
target_prob <- rowSums(target.b)/50
pvals <- c()
for(i in 1:4){
  pvals[i]<-binom.test(triplet_n[i], 50, target_prob[i]*tf_re_prob[i])$p.val
}

test_that("pruneRegulon calculates cluster p-values in binomial test correctly", {
  regulon.p <- pruneRegulon(regulon = regulon, expMatrix = geneM, peakMatrix = peakM, regulon_cutoff = 2,
                            clusters = rep(c("C1", "C2"), each = 50))
  expect_identical(regulon.p$pval_C1, pvals)
})
