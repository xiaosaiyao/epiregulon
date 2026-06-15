set.seed(1000)
mat <- matrix(rnorm(10000), ncol=100)
groups <- sample(LETTERS[1:5], ncol(mat), replace=TRUE)
sce <- SingleCellExperiment(list(logcounts=mat), colData=data.frame(Treatment=groups))
gene_names <- paste0('Gene_', 1:100)
rownames(sce) <- gene_names 

mat <- assay(sce, "logcounts")
groups <- sce$Treatment
pval.type <- "any"
min.prop <- NULL
method <- switch(pval.type, any="simes", some="holm-min", all="berger")

# calculate pvalues using t.test and p.adj; combine pvalues using metapod
my.groups <- sort(unique(groups))
ref.p <- list()
ref.q <- list()
cur.p.all <- list()
cur.q.all <- list()
cur.p.all.unlog <- list()
for (g1 in my.groups) {
  mat1 <- mat[,groups == g1]
  
  for (g2 in my.groups) {
    if (g1 == g2) {
      next
    }
    mat2 <- mat[,groups == g2]
    
    cur.p <- numeric(nrow(mat))
    for (r in seq_len(nrow(mat))) {
      cur.p[r] <- t.test(mat1[r,], mat2[r,], alternative="greater")$p.value
    }
    cur.p.all[[g1]][[g2]] <- log(cur.p)
    cur.q.all[[g1]][[g2]] <- log(p.adjust(cur.p, method="BH"))
    cur.p.all.unlog[[g1]][[g2]] <- cur.p
  }
  
  # combine pvalues
  ref.p.unlog <- combineParallelPValues(cur.p.all.unlog[[g1]], method=method, min.prop=min.prop)$p.value
  ref.p[[g1]] <- log(ref.p.unlog)
  ref.q[[g1]] <- log(p.adjust(ref.p.unlog, method="BH"))
  
  # compile pvalues from different groups 
  cur.p.all[[g1]] <- as.data.frame(cur.p.all[[g1]])
  rownames(cur.p.all[[g1]]) <- gene_names 
  colnames(cur.p.all[[g1]]) <- paste0("log.p.value.", colnames(cur.p.all[[g1]]))
  
  # compile FDR values from different groups
  cur.q.all[[g1]] <- as.data.frame(cur.q.all[[g1]])
  rownames(cur.q.all[[g1]]) <- gene_names 
  colnames(cur.q.all[[g1]]) <- paste0("log.FDR.", colnames(cur.q.all[[g1]]))
  
}

# calculate pvalues and FDR using pairwiseTestsSimple
compiled <- pairwiseTTestsSimple(mat, groups, direction="up")

test_that("pvalue in pairwiseTTestsSimple works correctly", {
  expect_equal(compiled$log.p.value, expected=cur.p.all)
})

test_that("FDR in pairwiseTTestsSimple works correctly", {
  expect_equal(compiled$log.FDR, expected=cur.q.all)
})


# calculate combined pvalues and FDR using combineMarkersSimple
combined <- combineMarkersSimple(compiled, groups, pval.type, min.prop)
combined.df <- as.data.frame(lapply(combined, function(x) {x[,"log.p.value"]}))
ref.p.df <- as.data.frame(ref.p)
combined.df <- combined.df[, colnames(ref.p.df)]

test_that("pvalue in combineMarkersSimple works correctly", {
  expect_equal(combined.df, expected=ref.p.df )
})


combined.FDR.df <- as.data.frame(lapply(combined, function(x) {x[,"log.FDR"]}))
ref.q.df <- as.data.frame(ref.q)
combined.FDR.df <- combined.FDR.df[, colnames(ref.q.df)]

test_that("FDR in combineMarkersSimple works correctly", {
  expect_equal(combined.FDR.df, expected=ref.q.df )
})



# compare findMarkersSimple vs scran::findMarkers
simple_output <- findMarkersSimple(sce, sce$Treatment, combined=FALSE, assay.type="logcounts", direction="any")
scran_output <- scran::findMarkers(x=sce, groups=sce$Treatment, test.type="t", 
                                   pval.type="any", full.stats=TRUE, sorted=FALSE, log.p=TRUE, direction="any")


scran_pval_A <- data.frame(log.p.value.B=scran_output$A$stats.B$log.p.value, 
                           log.p.value.C=scran_output$A$stats.C$log.p.value,
                           log.p.value.D=scran_output$A$stats.D$log.p.value,
                           log.p.value.E=scran_output$A$stats.E$log.p.value)
rownames(scran_pval_A) <- gene_names

scran_FDR_A <- data.frame(log.FDR.B=scran_output$A$stats.B$log.FDR, 
                          log.FDR.C=scran_output$A$stats.C$log.FDR,
                          log.FDR.D=scran_output$A$stats.D$log.FDR,
                          log.FDR.E=scran_output$A$stats.E$log.FDR)
rownames(scran_FDR_A) <- gene_names

scran_logFC_A <- data.frame(logFC.B=scran_output$A$stats.B$logFC, 
                            logFC.C=scran_output$A$stats.C$logFC,
                            logFC.D=scran_output$A$stats.D$logFC,
                            logFC.E=scran_output$A$stats.E$logFC)
rownames(scran_logFC_A) <- gene_names



test_that("findMarkersSimple gives the same pvalues as scran::findMarkers", {
  expect_equal(simple_output$log.p.value$A, expected=scran_pval_A )
})

test_that("findMarkersSimple gives the same FDR as scran::findMarkers", {
  expect_equal(simple_output$log.FDR$A, expected=scran_FDR_A )
})

test_that("findMarkersSimple gives the same logFC as scran::findMarkers", {
  expect_equal(simple_output$logFC$A, expected=scran_logFC_A )
})



simple_output2 <- findMarkersSimple(sce, sce$Treatment, combined=TRUE, pval.type="some",
                                    assay.type="logcounts", direction="any")
scran_output2 <- scran::findMarkers(x=sce, groups=sce$Treatment, test.type="t", 
                                   pval.type="some", full.stats=FALSE, sorted=FALSE, log.p=TRUE, direction="any")


test_that("findMarkersSimple in combined mode gives the same pvalues as scran::findMarkers", {
  expect_equal(simple_output2$A$log.p.value, expected=scran_output2$A$log.p.value)
})

test_that("findMarkersSimple in combined mode gives the same FDR as scran::findMarkers", {
  expect_equal(simple_output2$A$log.FDR, expected=scran_output2$A$log.FDR)
})

test_that("findMarkersSimple in combined mode gives the same logFC as scran::findMarkers", {
  expect_equal(simple_output2$A$logFC, expected=scran_output2$A$summary.logFC)
})

##########

set.seed(1000)
gene_sce <- scuttle::mockSCE()
gene_sce <- scrapper::normalizeRnaCounts.se(gene_sce)
rownames(gene_sce) <- paste0('Gene_',1:2000)

# create a mock regulon
regulon <- data.frame(tf=c(rep('Gene_1',10), rep('Gene_2',10)),
                      idxATAC=sample(1:100, 20),
                      target=c(paste0('Gene_', sample(3:2000,10)),
                                 paste0('Gene_',sample(3:2000,10))))

# filter regulon
pruned.regulon <- addLogFC(expMatrix=gene_sce, 
                           clusters=gene_sce$Treatment,
                           regulon=regulon,
                           sig_type="log.p.value", 
                           direction="any",
                           pval.type="any")

# test logFC of all conditions 
diff_exp <- scran::findMarkers(x=gene_sce, 
                               groups=gene_sce$Treatment, 
                               test.type="t", 
                               pval.type="any", 
                               full.stats=TRUE, 
                               sorted=FALSE,
                               log.p=TRUE)
diff_exp_df <- data.frame(matrix(data=NA, nrow=2000, ncol=4))
colnames(diff_exp_df) <- c("treat1.vs.rest.log.p.value", "treat1.vs.rest.logFC", 
                           "treat2.vs.rest.log.p.value", "treat2.vs.rest.logFC")
rownames(diff_exp_df) <- paste0("Gene_",1:2000)

diff_exp_df$treat1.vs.rest.log.p.value <- diff_exp$treat1$log.p.value
diff_exp_df$treat1.vs.rest.logFC <- diff_exp$treat1$summary.stats

diff_exp_df$treat2.vs.rest.log.p.value <- diff_exp$treat2$log.p.value
diff_exp_df$treat2.vs.rest.logFC <- diff_exp$treat2$summary.stats

rownames(regulon) <- regulon$target
combined_diff_exp_df <- cbind(regulon, diff_exp_df[regulon$target,])

test_that("addLogFC works correctly", {
  expect_equal(pruned.regulon, combined_diff_exp_df, tolerance=1e-10)
})


# test logFC of specific conditions

pruned.regulon2 <- addLogFC(expMatrix=gene_sce, 
                            clusters=gene_sce$Treatment,
                            regulon=regulon,
                            sig_type="log.p.value", 
                            pval.type="any", 
                            logFC_condition="treat1",
                            logFC_ref="treat2")

combined_diff_exp_df2 <- combined_diff_exp_df[,c("tf", "idxATAC", "target", "treat1.vs.rest.log.p.value", "treat1.vs.rest.logFC")]
colnames(combined_diff_exp_df2) <- colnames(pruned.regulon2)

test_that("addLogFC works correctly", {
  expect_identical(pruned.regulon2, combined_diff_exp_df2, tolerance=1e-10)
})
