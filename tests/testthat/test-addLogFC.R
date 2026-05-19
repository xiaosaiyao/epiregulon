set.seed(1000)
gene_sce <- scuttle::mockSCE()
gene_sce <- scrapper::normalizeRnaCounts.se(gene_sce)
rownames(gene_sce) <- paste0('Gene_',1:2000)

# create a mock regulon
regulon <- data.frame(tf = c(rep('Gene_1',10), rep('Gene_2',10)),
                     idxATAC = sample(1:100, 20),
                     target = c(paste0('Gene_', sample(3:2000,10)),
                                paste0('Gene_',sample(3:2000,10))))

# filter regulon
pruned.regulon <- addLogFC(expMatrix = gene_sce, 
                           clusters = gene_sce$Treatment,
                           regulon = regulon,
                           sig_type = "p.value", 
                           pval.type = "any")

# test logFC of all conditions 
diff_exp <- scran::findMarkers(x=gene_sce, groups=gene_sce$Treatment, test.type = "t", pval.type = "any", full.stats = TRUE, sorted=FALSE)
diff_exp_df <- data.frame(matrix(data=NA, nrow=2000, ncol=4))
colnames(diff_exp_df) <- c("treat1.vs.rest.p.value", "treat1.vs.rest.logFC", "treat2.vs.rest.p.value", "treat2.vs.rest.logFC")
rownames(diff_exp_df) <- paste0("Gene_",1:2000)

diff_exp_df$treat1.vs.rest.p.value <- diff_exp$treat1$p.value
diff_exp_df$treat1.vs.rest.logFC <- diff_exp$treat1$summary.stats

diff_exp_df$treat2.vs.rest.p.value <- diff_exp$treat2$p.value
diff_exp_df$treat2.vs.rest.logFC <- diff_exp$treat2$summary.stats

rownames(regulon) <- regulon$target
combined_diff_exp_df <- cbind(regulon, diff_exp_df[regulon$target,])

test_that("addLogFC works correctly", {
  expect_identical(pruned.regulon, combined_diff_exp_df)
})


# test logFC of specific conditions

pruned.regulon2 <- addLogFC(expMatrix = gene_sce, 
                           clusters = gene_sce$Treatment,
                           regulon = regulon,
                           sig_type = "p.value", 
                           pval.type = "any", 
                           logFC_condition = "treat1",
                           logFC_ref = "treat2")

combined_diff_exp_df2 <- combined_diff_exp_df[,c("tf", "idxATAC", "target", "treat1.vs.rest.p.value", "treat1.vs.rest.logFC")]
colnames(combined_diff_exp_df2) <- colnames(pruned.regulon2)

test_that("addLogFC works correctly", {
  expect_identical(pruned.regulon2, combined_diff_exp_df2)
})
