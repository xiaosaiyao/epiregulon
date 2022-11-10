geneM <- matrix(0, nrow = 4, ncol = 100, dimnames = list(LETTERS[1:4], NULL))
clusters <- rep(c("C1", "C2"), each = 50)
geneM[c("A", "D"),1:50] <- rep(c(1.1,2,3,4,5),10)
geneM[c("A", "C"), 51:100] <- rep(c(1.1,0,0,2,0), 10)
geneM["C", 1:50] <- rep(c(2,3,0,0,0),10)
geneM["B",] <- rep(c(0,1.1,0,0,0,2,3,0,0,0),10)
geneM["A",][c(1:30)*3] <- 0
geneM <- as(geneM, "sparseMatrix")

peakM <- matrix(0, nrow = 4, ncol =100, dimnames = list(LETTERS[1:4], NULL))
peakM[c(1,2),] <- rep(c(0,1), 50)
peakM[3,] <- rep(c(1,0,0,0),25)
peakM[2, 1:50] <- rep(c(0,1,1,0,1,1,0,0,0,0),5)
peakM[2, 51:100] <- rep(c(0,4,0,0,0,0,0,0,4,0),5)
peakM[4, 1:50] <- rep(c(1,1,1,1,0), 10)
peakM[4, 51:100] <- rep(c(0,0,1,1,0),10)
peakM <- as(peakM, "sparseMatrix")

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


test_that("pruneRegulon calculates p-values with binomial test correctly", {
  regulon.p <- pruneRegulon(regulon = regulon, expMatrix = geneM, peakMatrix = peakM, regulon_cutoff = 2)
  expect_identical(regulon.p$pval_all, pvals)
})

pvals <- c()
for(i in 1:4){
  pvals[i]<-suppressWarnings(chisq.test(c(triplet_n[i], 100-triplet_n[i]),
                       p =c(target_prob[i]*tf_re_prob[i] ,
                            1- target_prob[i]*tf_re_prob[i]))$p.val)
}

test_that("pruneRegulon calculates p-values with chi-square test correctly", {
  regulon.p <- pruneRegulon(regulon = regulon, expMatrix = geneM, peakMatrix = peakM,
                            regulon_cutoff = 2, test ="chi.sq")
  expect_identical(regulon.p$pval_all, pvals[is.finite(pvals)])
})

# calculating results for the cluster 1
peakM.b_C1 <- peakM.b[,1:50]
geneM.b_C1 <- geneM.b[,1:50]
tf_re.b_C1 <- geneM.b_C1[regulon$tf,]*peakM.b_C1[regulon$idxATAC,]
triplets.b_C1 <- tf_re.b_C1*geneM.b_C1[regulon$target,]
target.b_C1 <- geneM.b_C1[regulon$target,]
triplet_n_C1 <- rowSums(triplets.b_C1)
tf_re_prob_C1 <- rowSums(tf_re.b_C1)/50
target_prob_C1 <- rowSums(target.b_C1)/50
pvals_C1 <- c()
for(i in 1:4){
  pvals_C1[i]<-binom.test(triplet_n_C1[i], 50, target_prob_C1[i]*tf_re_prob_C1[i])$p.val
}

# calculating results for the cluster 2
peakM.b_C2 <- peakM.b[,51:100]
geneM.b_C2 <- geneM.b[,51:100]
tf_re.b_C2 <- geneM.b_C2[regulon$tf,]*peakM.b_C2[regulon$idxATAC,]
triplets.b_C2 <- tf_re.b_C2*geneM.b_C2[regulon$target,]
target.b_C2 <- geneM.b_C2[regulon$target,]
triplet_n_C2 <- rowSums(triplets.b_C2)
tf_re_prob_C2 <- rowSums(tf_re.b_C2)/50
target_prob_C2 <- rowSums(target.b_C2)/50
pvals_C2 <- c()
for(i in 1:4){
  pvals_C2[i]<-binom.test(triplet_n_C2[i], 50, target_prob_C2[i]*tf_re_prob_C2[i])$p.val
}

test_that("pruneRegulon calculates cluster p-values with binomial test correctly", {
  regulon.p <- pruneRegulon(regulon = regulon, expMatrix = geneM, peakMatrix = peakM, regulon_cutoff = 2,
                            clusters = rep(c("C1", "C2"), each = 50))
  expect_identical(regulon.p$pval_C1, pvals_C1)
  expect_identical(regulon.p$pval_C2, pvals_C2)
})

pvals_C1 <- c()
for(i in 1:4){
  pvals_C1[i]<-suppressWarnings(chisq.test(c(triplet_n_C1[i], 50-triplet_n_C1[i]),
                                        p =c(target_prob_C1[i]*tf_re_prob_C1[i] ,
                                             1- target_prob_C1[i]*tf_re_prob_C1[i]))$p.val)
}

pvals_C2 <- c()
for(i in 1:4){
  pvals_C2[i]<-suppressWarnings(chisq.test(c(triplet_n_C2[i], 50-triplet_n_C2[i]),
                                           p =c(target_prob_C2[i]*tf_re_prob_C2[i] ,
                                                1- target_prob_C2[i]*tf_re_prob_C2[i]))$p.val)
}

# create data frame with
pvals <- data.frame(C1 = pvals_C1, C2 = pvals_C2)

# remove rows with NaN values
pvals <- suppressWarnings(pvals[is.finite(apply(pvals,1,function(x) min(x, na.rm = TRUE))),])

test_that("pruneRegulon calculates cluster p-values with chi-square test correctly", {
  regulon.p <- pruneRegulon(regulon = regulon, expMatrix = geneM, peakMatrix = peakM, regulon_cutoff = 2,
                            clusters = rep(c("C1", "C2"), each = 50), test = "chi.sq")
  expect_identical(regulon.p$pval_C1, pvals$C1)
  expect_identical(regulon.p$pval_C2, pvals$C2)
})
