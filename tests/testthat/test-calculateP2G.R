set.seed(1010)
# set up expression and peak matrices
n_cells <- 1000
n_genes <- 200
n_peaks <- 500
n_cell_states <- 20
cell_states <- sample(1:20, n_cells, replace = TRUE)
expression_means <- matrix(runif(n_genes*n_cell_states)*3, ncol=n_cell_states)
expression_means <- expression_means[,cell_states]
geneExpMatrix <- matrix(rnorm(length(expression_means), mean=expression_means), ncol = n_cells,
                        dimnames = list(NULL, paste0("Cell_", seq_len(n_cells))))
geneExpMatrix[geneExpMatrix < 0] <- 0
geneExpMatrix <- as(geneExpMatrix, "sparseMatrix")
rownames(geneExpMatrix) <- paste0("Gene_", seq_len(nrow(geneExpMatrix)))

gene.ranges <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr3", "chr4"), 50),
                   ranges = IRanges(start = seq(from = 1, length.out=n_genes, by = 1000),
                                    width = 100))
peak.ranges <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr3", "chr4"), c(90,160,200,50)),
                       ranges = IRanges(start = runif(n_peaks)*2e5, width = 100))

gene.start <- resize(gene.ranges, width=1)
peak.ranges <- peak.ranges
overlap <- S4Vectors::DataFrame(findOverlaps(resize(gene.start, 5000, "center"),
                                             peak.ranges))
colnames(overlap) <- c("RNA", "ATAC")


regulatory_links_idx <- sample(nrow(overlap), round(nrow(overlap)*0.3))
regulatory_pairs <- overlap[regulatory_links_idx,]

peakMatrix <- matrix(0, nrow = n_peaks,
                     ncol = n_cells,
                     dimnames = list(NULL, paste0("Cell_", seq_len(n_cells))))

for(peak_idx in unique(regulatory_pairs[,2])){
    target_genes <- unique(regulatory_pairs[,1][regulatory_pairs[,2]==peak_idx])
    target_expression <- Matrix::colSums(geneExpMatrix[target_genes,,drop=FALSE])
    target_expression_norm <- target_expression/max(target_expression)
    peakMatrix[peak_idx,] <- rbinom(length(target_expression_norm), 1, target_expression_norm)
}

# add sparsity
peakMatrix[sample(length(peakMatrix), round(length(peakMatrix)*0.3))] <- 0
geneExpMatrix[sample(length(geneExpMatrix), round(length(geneExpMatrix)*0.3))] <- 0

non.zero.genes <- which(Matrix::rowSums(geneExpMatrix) != 0)
non.zero.peaks <- which(Matrix::rowSums(peakMatrix) != 0)

new_gene_idx <- seq_along(non.zero.genes)
new_peak_idx <- seq_along(non.zero.peaks)

overlap <- overlap[(overlap[,1] %in% non.zero.genes) & (overlap[,2] %in% non.zero.peaks),]

overlap[,1] <- new_gene_idx[match(overlap[,1], sort(non.zero.genes))]
overlap[,2] <- new_peak_idx[match(overlap[,2], sort(non.zero.peaks))]
geneExpMatrix <- geneExpMatrix[non.zero.genes,]
peakMatrix <- peakMatrix[non.zero.peaks,]
gene.ranges <- gene.ranges[non.zero.genes]
peak.ranges <- peak.ranges[non.zero.peaks]

# calculate null distributions
null_correlations <- c()
for(peak_idx in unique(overlap[,2])){
    peak_chromosome = seqnames(peak.ranges[peak_idx])
    n_rep <- sum(overlap[,2]==peak_idx)
    distant_genes_idx <- which(as.logical(seqnames(gene.ranges)!=peak_chromosome))
    selected_genes <- sample(distant_genes_idx, 200*n_rep,replace=TRUE)
    for(j in seq_along(selected_genes)){
        null_correlations <- c(null_correlations, cor(peakMatrix[peak_idx,], geneExpMatrix[selected_genes[j],]))
    }
}

overlap$Correlation <- matrix(NA, nrow=nrow(overlap), ncol=1)
colnames(overlap$Correlation) <- "all"
for(i in seq_len(nrow(overlap))){
    overlap$Correlation[i,"all"] <- cor(peakMatrix[overlap[i,2],], geneExpMatrix[overlap[i,1],])
}

df <- overlap

overlap$p_val <- matrix(1, nrow=nrow(overlap), ncol=1)
colnames(overlap$p_val) <- "all"
non_neg_cor_idx <- which(overlap$Correlation[,"all"]>=0)
non_neg_corr_null <- null_correlations[null_correlations>=0]
for(i in non_neg_cor_idx){
    overlap$p_val[i,"all"] <- sum(non_neg_corr_null > overlap$Correlation[i,"all"])/length(non_neg_corr_null)
}

non_pos_cor_idx <- which(overlap$Correlation<=0)
non_pos_corr_null <- null_correlations[null_correlations<=0]
for(i in non_pos_cor_idx){
    overlap$p_val[i,"all"] <- sum(non_pos_corr_null < overlap$Correlation[i,"all"])/length(non_pos_corr_null)
}

overlap$FDR <- matrix(1, nrow=nrow(overlap), ncol=1)
colnames(overlap$FDR) <- "all"
overlap$FDR[, "all"] <- p.adjust(overlap$p_val[,"all"],method="BH")

stat_list <- .addFDR(overlap, geneStart = gene.ranges, peakSet = peak.ranges,
              geneExpr = geneExpMatrix, peakCounts = peakMatrix,
              n_random_conns = 1e5,
              cor_method = "pearson",
              batch_size=2e4,
              BPPARAM=BiocParallel::MulticoreParam())

df <- overlap
df$p_val <- matrix(stat_list$p_val, nrow=nrow(overlap), ncol=1)
colnames(df$p_val) <- "all"
df$FDR <- matrix(stat_list$FDR, nrow=nrow(overlap), ncol=1)
colnames(df$FDR) <- "all"



test_that(".addFDR works correctly", {
    expect_equal(df$p_val[,"all"], overlap$p_val[,"all"], tolerance = 2e-2)
    expect_true(cor(overlap$p_val[,"all"],df$p_val[,"all"])>0.9999)
    expect_true(cor(overlap$FDR[,"all"],df$FDR[,"all"])>0.999)
    expect_equal(df$Correlation[,"all"], overlap$Correlation[,"all"])
})
mcols(gene.ranges)$name <- rownames(geneExpMatrix)
peakMatrix_sce <- SingleCellExperiment(assay=list(counts=peakMatrix), rowRanges=peak.ranges)
geneExpMatrix_sce <- SingleCellExperiment(assay=list(counts=geneExpMatrix), rowRanges=gene.ranges)
cellNum <- optimizeMetacellNumber(peakMatrix_sce, geneExpMatrix_sce,
                                  reducedDim=Matrix::t(geneExpMatrix), exp_assay="counts",
                      peak_assay="counts", subsample_prop=0.1,
                      n_iter=1, cellNumMin=NULL,
                      cellNumMax=NULL, n_evaluation_points=4)

min_eval_point <- sqrt(min(20, round(ncol(peakMatrix)/10)))
max_eval_point <- sqrt(min(2000, round(ncol(peakMatrix)/10)))
test_that("optimizeMetacellNumber works correctly", {
    expect_s4_class(cellNum, "CellNumSol")
    expect_equal(max(cellNum@evaluation_points), max_eval_point)
    expect_equal(min(cellNum@evaluation_points), min_eval_point)
    expect_length(cellNum@AUC, length(cellNum@evaluation_points))
})
