geneM <- matrix(0, nrow = 4, ncol = 100, dimnames = list(LETTERS[1:4], NULL))
clusters <- rep(c("C1", "C2"), each = 50)
geneM[c("A", "D"),1:50] <- rep(c(1.1,2,3,4,5),10)
geneM[c("A", "C"), 51:100] <- rep(c(1.1,0,0,2,0), 10)
geneM["C", 1:50] <- rep(c(2,3,0,0,0),10)
geneM["B",] <- rep(c(0,1.1,0,0,0,2,3,0,0,0),10)
geneM["A",][c(1:30)*3] <- 0
geneM["B", c(6,10,19,20)] <- 1.5 # decrease p-val for A-2-B triplet, cluster 1
geneM <- as(geneM, "sparseMatrix")

peakM <- matrix(0, nrow = 4, ncol = 100, dimnames = list(LETTERS[1:4], NULL))
peakM[c(1,2),] <- rep(c(0,1), 50)
peakM[3,] <- rep(c(1,0,0,0),25)
peakM[2, 1:50] <- rep(c(0,1,1,0,1,1,0,0,0,0),5)
peakM[2, 51:100] <- rep(c(0,4,0,0,0,0,0,0,4,0),5)
peakM[4, 1:50] <- rep(c(1,1,1,1,0), 10)
peakM[4, 51:100] <- rep(c(0,0,1,1,0),10)
peakM[2, c(10,17,19,20, 47)] <- 2 # decrease p-val for A-2-B triplet, cluster 1
peakM <- as(peakM, "sparseMatrix")


regulon <- data.frame(tf = c("A","A","A","C"), target = c("B", "B", "C", "B"), idxATAC = c(1:4))

peakM.b <- peakM > 0
geneM.b <- geneM > 1
tf_re.b <- geneM.b[regulon$tf,]*peakM.b[regulon$idxATAC,]
triplets.b <- tf_re.b*geneM.b[regulon$target,]
target.b <- geneM.b[regulon$target,]
triplet_n <- Matrix::rowSums(triplets.b)
tf_re_n <- Matrix::rowSums(tf_re.b)
target_n <- Matrix::rowSums(target.b)
null_probs <-target_n*tf_re_n/100^2
pvals <- c()
for(i in 1:4){
  pvals[i] <- stats::binom.test(triplet_n[i], 100, null_probs[i])$p.val
}


test_that("pruneRegulon calculates p-values with binomial test correctly", {
  regulon.p <- pruneRegulon(regulon = regulon,
                            expMatrix = SingleCellExperiment(assays= list(logcounts=geneM)),
                            peakMatrix = SingleCellExperiment(assays= list(PeakMatrix=peakM)),
                            regulon_cutoff = 2,
                            test = "binom")
  expect_identical(unname(regulon.p$pval[,"all"]), pvals, tolerance = 1e-3)
})

pvals <- c()
for(i in 1:4){
  pvals[i] <- suppressWarnings(chisq.test(c(triplet_n[i], 100-triplet_n[i]),
                                          p = c(null_probs[i],1- null_probs[i]))$p.val)
  pvals[i][is.na(pvals[i])] <- 1
}

test_that("pruneRegulon calculates p-values with chi-square test correctly", {
  regulon.p <- pruneRegulon(regulon = regulon,
                            expMatrix = SingleCellExperiment(assays= list(logcounts=geneM)),
                            peakMatrix = SingleCellExperiment(assays= list(PeakMatrix=peakM)),
                            regulon_cutoff = 2,
                            test ="chi.sq")

  expect_identical(unname(stats::na.omit(regulon.p$pval[,"all"])),
                   pvals[is.finite(pvals)], tolerance = 1e-3)
})

# calculating results for the cluster 1
peakM.b_C1 <- peakM.b[,1:50]
geneM.b_C1 <- geneM.b[,1:50]
tf_re.b_C1 <- geneM.b_C1[regulon$tf,]*peakM.b_C1[regulon$idxATAC,]
triplets.b_C1 <- tf_re.b_C1*geneM.b_C1[regulon$target,]
target.b_C1 <- geneM.b_C1[regulon$target,]
triplet_n_C1 <- Matrix::rowSums(triplets.b_C1)
tf_re_n_C1 <- Matrix::rowSums(tf_re.b_C1)
target_n_C1 <- Matrix::rowSums(target.b_C1)
null_probs_C1 <- tf_re_n_C1*target_n_C1/50^2
pvals_C1 <- c()
for(i in 1:4){
  pvals_C1[i]<- binom.test(triplet_n_C1[i], 50, null_probs_C1[i])$p.val
}

# calculating results for the cluster 2
peakM.b_C2 <- peakM.b[,51:100]
geneM.b_C2 <- geneM.b[,51:100]
tf_re.b_C2 <- geneM.b_C2[regulon$tf,]*peakM.b_C2[regulon$idxATAC,]
triplets.b_C2 <- tf_re.b_C2*geneM.b_C2[regulon$target,]
target.b_C2 <- geneM.b_C2[regulon$target,]
triplet_n_C2 <- Matrix::rowSums(triplets.b_C2)
tf_re_n_C2 <- Matrix::rowSums(tf_re.b_C2)
target_n_C2 <- Matrix::rowSums(target.b_C2)
null_probs_C2 <- tf_re_n_C2*target_n_C2/50^2
pvals_C2 <- c()
for(i in 1:4){
  pvals_C2[i] <- binom.test(triplet_n_C2[i], 50, null_probs_C2[i])$p.val
}

test_that("pruneRegulon calculates cluster p-values with binomial test correctly", {
  regulon.p <- pruneRegulon(regulon = regulon,
                            expMatrix = SingleCellExperiment(assays= list(logcounts=geneM)),
                            peakMatrix = SingleCellExperiment(assays= list(PeakMatrix=peakM)),
                            regulon_cutoff = 2,
                            clusters = rep(c("C1", "C2"), each = 50),
                            test = "binom")
  expect_identical(unname(regulon.p$pval[,"C1"]), pvals_C1, tolerance = 1e-3)
  expect_identical(unname(regulon.p$pval[,"C2"]), pvals_C2, tolerance = 1e-3)
})

pvals_C1 <- c()
for(i in 1:4){
  pvals_C1[i]<-suppressWarnings(chisq.test(c(triplet_n_C1[i], 50-triplet_n_C1[i]),
                                           p = c(null_probs_C1[i], 1- null_probs_C1[i]))$p.val)
}

pvals_C2 <- c()
for(i in 1:4){
  pvals_C2[i]<- suppressWarnings(chisq.test(c(triplet_n_C2[i], 50-triplet_n_C2[i]),
                                           p = c(null_probs_C2[i], 1- null_probs_C2[i]))$p.val)
}

pvals_C1[is.na(pvals_C1)] <- 1
pvals_C2[is.na(pvals_C2)] <- 1

# create data frame with
pvals <- data.frame(C1 = pvals_C1, C2 = pvals_C2)

# remove rows with NaN values
pvals <- suppressWarnings(pvals[is.finite(apply(pvals,1,function(x) min(x, na.rm = TRUE))),])

test_that("pruneRegulon calculates cluster p-values with chi-square test correctly", {
  regulon.p <- pruneRegulon(regulon = regulon,
                            expMatrix = SingleCellExperiment(assays= list(logcounts=geneM)),
                            peakMatrix = SingleCellExperiment(assays= list(PeakMatrix=peakM)),
                            regulon_cutoff = 2,
                            clusters = rep(c("C1", "C2"), each = 50),
                            test = "chi.sq")
  expect_identical(unname(regulon.p$pval[,"C1"]), pvals$C1, tolerance = 1e-8)
  expect_identical(unname(regulon.p$pval[,"C2"]), pvals$C2, tolerance = 1e-8)
})


selected_rows <- apply(pvals, 1, function(x) min(x, na.rm = TRUE)<0.05)
pvals <- pvals[selected_rows,]

test_that("pruneRegulon correctly applies 'regulon_cutoff'", {
  regulon.p <- pruneRegulon(regulon = regulon,
                            expMatrix = SingleCellExperiment(assays= list(logcounts=geneM)),
                            peakMatrix = SingleCellExperiment(assays= list(PeakMatrix=peakM)),
                            regulon_cutoff = 0.05,
                            clusters = rep(c("C1", "C2"), each = 50),
                            test = "chi.sq")
  expect_identical(unname(regulon.p$pval[,"C1"]), pvals$C1, tolerance = 1e-6)
  expect_identical(unname(regulon.p$pval[,"C2"]), pvals$C2, tolerance = 1e-6)
})


# applying moving cutoff
peakM.b <- peakM > rowMeans(as.matrix(peakM))
geneM.b <- geneM > rowMeans(as.matrix(geneM))
tf_re.b <- geneM.b[regulon$tf,]*peakM.b[regulon$idxATAC,]
triplets.b <- tf_re.b*geneM.b[regulon$target,]
target.b <- geneM.b[regulon$target,]
triplet_n <- Matrix::rowSums(triplets.b)
tf_re_n <- Matrix::rowSums(tf_re.b)
target_n <- Matrix::rowSums(target.b)
null_probs <-target_n*tf_re_n/100^2
pvals <- c()
for(i in 1:4){
  pvals[i] <- stats::binom.test(triplet_n[i], 100, null_probs[i])$p.val
}


test_that("pruneRegulon calculates p-values with binomial test correctly when using moving cutoffs", {
  regulon.p <- pruneRegulon(regulon = regulon,
                            expMatrix = SingleCellExperiment(assays= list(logcounts=geneM)),
                            peakMatrix = SingleCellExperiment(assays= list(PeakMatrix=peakM)),
                            regulon_cutoff = 2,
                            test = "binom",
                            peak_cutoff = NULL,
                            exp_cutoff = NULL)
  expect_identical(unname(regulon.p$pval[,"all"]), pvals, tolerance = 1e-6)
})


# calculating results for the cluster 1
peakM.b_C1 <- peakM[,1:50] > Matrix::rowMeans(peakM[,1:50])
geneM.b_C1 <- geneM[,1:50] > Matrix::rowMeans(geneM[,1:50])
tf_re.b_C1 <- geneM.b_C1[regulon$tf,]*peakM.b_C1[regulon$idxATAC,]
triplets.b_C1 <- tf_re.b_C1*geneM.b_C1[regulon$target,]
target.b_C1 <- geneM.b_C1[regulon$target,]
triplet_n_C1 <- Matrix::rowSums(triplets.b_C1)
tf_re_n_C1 <- Matrix::rowSums(tf_re.b_C1)
target_n_C1 <- Matrix::rowSums(target.b_C1)
null_probs_C1 <- tf_re_n_C1*target_n_C1/50^2
pvals_C1 <- c()
for(i in 1:4){
  pvals_C1[i]<-suppressWarnings(chisq.test(c(triplet_n_C1[i], 50-triplet_n_C1[i]),
                                           p = c(null_probs_C1[i], 1- null_probs_C1[i]))$p.val)
}

# calculating results for the cluster 2

peakM.b_C2 <- peakM[,51:100] > Matrix::rowMeans(peakM[,51:100])
geneM.b_C2 <- geneM[,51:100] > Matrix::rowMeans(geneM[,51:100])
tf_re.b_C2 <- geneM.b_C2[regulon$tf,]*peakM.b_C2[regulon$idxATAC,]
triplets.b_C2 <- tf_re.b_C2*geneM.b_C2[regulon$target,]
target.b_C2 <- geneM.b_C2[regulon$target,]
triplet_n_C2 <- Matrix::rowSums(triplets.b_C2)
tf_re_n_C2 <- Matrix::rowSums(tf_re.b_C2)
target_n_C2 <- Matrix::rowSums(target.b_C2)
null_probs_C2 <- tf_re_n_C2*target_n_C2/50^2
pvals_C2 <- c()
for(i in 1:4){
  pvals_C2[i]<-suppressWarnings(chisq.test(c(triplet_n_C2[i], 50-triplet_n_C2[i]),
                                           p = c(null_probs_C2[i], 1- null_probs_C2[i]))$p.val)
}

# calculating results for all cells

peakM.b <- peakM > Matrix::rowMeans(peakM)
geneM.b <- geneM > Matrix::rowMeans(geneM)
tf_re.b <- geneM.b[regulon$tf,]*peakM.b[regulon$idxATAC,]
triplets.b <- tf_re.b*geneM.b[regulon$target,]
target.b <- geneM.b[regulon$target,]
triplet_n <- Matrix::rowSums(triplets.b)
tf_re_n <- Matrix::rowSums(tf_re.b)
target_n <- Matrix::rowSums(target.b)
null_probs <- tf_re_n*target_n/100^2
pvals_all <- c()
for(i in 1:4){
  pvals_all[i]<-suppressWarnings(chisq.test(c(triplet_n[i], 100-triplet_n[i]),
                                            p = c(null_probs[i], 1- null_probs[i]))$p.val)
}

pvals_C1[is.na(pvals_C1)] <- 1
pvals_C2[is.na(pvals_C2)] <- 1
pvals_all[is.na(pvals_all)] <- 1
pvals <- data.frame(C1 = pvals_C1, C2 = pvals_C2, all = pvals_all)
# remove rows with NaN values
pvals <- suppressWarnings(pvals[is.finite(apply(pvals,1,function(x) min(x, na.rm = TRUE))),])

test_that("pruneRegulon calculates cluster p-values with chi-square test correctly when moving cutoffs are used", {
  regulon.p <- pruneRegulon(regulon = regulon,
                            expMatrix = SingleCellExperiment(assays= list(logcounts=geneM)),
                            peakMatrix = SingleCellExperiment(assays= list(PeakMatrix=peakM)),
                            regulon_cutoff = 2,
                            clusters = rep(c("C1", "C2"), each = 50),
                            test = "chi.sq",
                            exp_cutoff = NULL,
                            peak_cutoff = NULL)
  expect_identical(unname(regulon.p$pval[,"C1"]), pvals$C1, tolerance = 1e-8)
  expect_identical(unname(regulon.p$pval[,"C2"]), pvals$C2, tolerance = 1e-8)
  expect_identical(unname(regulon.p$pval[,"all"]), pvals$all, tolerance = 1e-8)
})
