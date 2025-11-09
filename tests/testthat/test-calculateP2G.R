set.seed(1010)
# set up expression and peak matrices
n_cells <- 1000
n_genes <- 200
n_peaks <- 500
n_cell_states <- 20
cellNum = 10
cell_states <- sample(1:20, n_cells, replace = TRUE)
expression_means <- matrix(runif(n_genes*n_cell_states)*3, ncol=n_cell_states)
expression_means <- expression_means[,cell_states]
geneExpMatrix <- matrix(rnorm(length(expression_means), mean=expression_means), ncol = n_cells,
                        dimnames = list(NULL, paste0("Cell_", seq_len(n_cells))))
geneExpMatrix[geneExpMatrix < 0] <- 0
geneExpMatrix <- as(geneExpMatrix, "sparseMatrix")
rownames(geneExpMatrix) <- gene_names <- paste0("Gene_", seq_len(nrow(geneExpMatrix)))

gene.ranges <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr3", "chr4"), 50),
                   ranges = IRanges(start = seq(from = 1, length.out=n_genes, by = 1000),
                                    width = 100), gene_name = gene_names)
peak.ranges <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr3", "chr4"), c(90,160,200,50)),
                       ranges = IRanges(start = runif(n_peaks)*2e5, width = 100))

gene.start <- resize(gene.ranges, width=1)
overlap <- S4Vectors::DataFrame(findOverlaps(resize(gene.start, 5000, "center"),
                                             peak.ranges))
colnames(overlap) <- c("RNA", "ATAC")

regulatory_links_idx <- sample(nrow(overlap), round(nrow(overlap)*0.3))
regulatory_pairs <- overlap[regulatory_links_idx,]

peakMatrix <- matrix(0, nrow = n_peaks,
                     ncol = n_cells,
                     dimnames = list(NULL, paste0("Cell_", seq_len(n_cells))))

# relate peak accessibility to gene expression
for(peak_idx in unique(regulatory_pairs[,2])){
    target_genes <- unique(regulatory_pairs[,1][regulatory_pairs[,2]==peak_idx])
    target_expression <- Matrix::colSums(geneExpMatrix[target_genes,,drop=FALSE])
    # gene expressed close to the max value is more likely to be related with
    # the open chromatin region
    target_expression_norm <- target_expression/max(target_expression)
    peakMatrix[peak_idx,] <- rbinom(length(target_expression_norm), 1, target_expression_norm)
}

# add sparsity
peakMatrix[sample(length(peakMatrix), round(length(peakMatrix)*0.3))] <- 0
geneExpMatrix[sample(length(geneExpMatrix), round(length(geneExpMatrix)*0.3))] <- 0

non.zero.genes <- which(Matrix::rowSums(geneExpMatrix) != 0)
non.zero.peaks <- which(Matrix::rowSums(peakMatrix) != 0)


###### tests .create_metacells
k = round(n_cells/cellNum)
kclusters <- as.numeric(clusterKmeans(as.matrix(geneExpMatrix),k = k)$clusters)

aggregated_gene_expr <- matrix(NA, nrow=n_genes, ncol=k, dimnames = list(rownames(geneExpMatrix), 1:k))
aggregated_peak_counts <- matrix(NA, nrow=n_peaks, ncol=k, dimnames = list(rownames(peakMatrix), 1:k))
for(i in 1:k){
    aggregated_gene_expr[,i] <- Matrix::rowSums(geneExpMatrix[,kclusters==i,drop=FALSE])/sum(kclusters==i)
    aggregated_peak_counts[,i] <- Matrix::rowSums(peakMatrix[,kclusters==i,drop=FALSE])/sum(kclusters==i)
}

aggregated_gene_expr <- aggregated_gene_expr[non.zero.genes,]
aggregated_peak_counts <- aggregated_peak_counts[non.zero.peaks,]
gene.start_filtered <- gene.start[non.zero.genes]
names(gene.start_filtered) <- mcols(gene.start_filtered)$gene_name

res_list_expected <-
    list(geneExpr = aggregated_gene_expr,
         peakCounts = aggregated_peak_counts,
         geneStart = gene.start_filtered,
         peakSet = peak.ranges[non.zero.peaks],
         old.idxRNA = as.integer(non.zero.genes),
         old.idxATAC = as.integer(non.zero.peaks),
         clust = as.factor(kclusters))

peakMatrix_sce <- SingleCellExperiment(assay=list(counts=peakMatrix), rowRanges=peak.ranges)
geneExpMatrix_sce <- SingleCellExperiment(assay=list(counts=geneExpMatrix), rowRanges=gene.ranges)

res_list <- .create_metacells(expMatrix=geneExpMatrix_sce,
                              exp_assay="counts",
                              peakMatrix=peakMatrix_sce,
                              peak_assay="counts",
                              reducedDim=Matrix::t(geneExpMatrix),
                              gene_symbol="gene_name",
                              frac_RNA=0,
                              frac_ATAC=0,
                              kNum=k)

test_that(".create_metacells works correctly", {
    expect_identical(res_list,res_list_expected)
})


###### test correlations with user-specified clusters #################
# remove non-expressed features
new_gene_idx <- seq_along(non.zero.genes)
new_peak_idx <- seq_along(non.zero.peaks)

overlap_2 <- S4Vectors::DataFrame(findOverlaps(resize(gene.start, 30000, "center"),
                                             peak.ranges))
colnames(overlap_2) <- c("RNA", "ATAC")

overlap_2 <- overlap_2[(overlap_2[,1] %in% non.zero.genes) & (overlap_2[,2] %in% non.zero.peaks),]

overlap_2[,1] <- new_gene_idx[match(overlap_2[,1], sort(non.zero.genes))]
overlap_2[,2] <- new_peak_idx[match(overlap_2[,2], sort(non.zero.peaks))]
overlap_2 <- overlap_2[order(overlap_2$ATAC, overlap_2$RNA),]



clusters <- rep(NA, length(cell_states))
# 4 cell states is one cluster
for(i in seq(1,16,5)){
    clusters[cell_states>=i] <- letters[(16-i)/5+1]
}

unique_clusters <- unique(clusters)
metacell_clusters <- list()
# determine cluster composition in terms of the kclusters
for(cluster in unique_clusters){
    metacell_clusters[[cluster]] <- c()
    for(kcluster in unique(kclusters)){
        # assign kcluster to cluster if its cells are overrepresented in the cluster
        if((sum((kclusters==kcluster) * (clusters==cluster))/(sum(kclusters==kcluster)))>=(1/(length(unique_clusters)))){
            metacell_clusters[[cluster]] <- c(metacell_clusters[[cluster]], kcluster)
        }
    }
}

# initialize carrelation_matrix
corr_matrix <- matrix(NA, nrow=nrow(overlap_2), ncol=length(unique_clusters)+1)
colnames(corr_matrix) <- c("all", sort(unique_clusters))
for(i in seq_len(nrow(overlap_2))){
    corr_matrix[i,1] <- cor(aggregated_gene_expr[overlap_2$RNA[i],], aggregated_peak_counts[overlap_2$ATAC[i],])
}
for(cluster in unique_clusters){
    cluster_composition <- metacell_clusters[[cluster]]
    for(j in 1:nrow(overlap_2)){
        corr_matrix[j,cluster] <- cor(aggregated_gene_expr[overlap_2$RNA[j],cluster_composition],
                                      aggregated_peak_counts[overlap_2$ATAC[j],cluster_composition])

    }
}

peakMatrix_sce <- SingleCellExperiment(assay=list(counts=peakMatrix), rowRanges=peak.ranges)
geneExpMatrix_sce <- SingleCellExperiment(assay=list(counts=geneExpMatrix), rowRanges=gene.ranges)

p2g <- calculateP2G(peakMatrix = peakMatrix_sce,
                    expMatrix = geneExpMatrix_sce,
                    reducedDim = as.matrix(Matrix::t(geneExpMatrix)),
                    exp_assay = "counts",
                    peak_assay = "counts",
                    gene_symbol = "gene_name",
                    cellNum=cellNum,
                    frac_RNA = 0,
                    frac_ATAC = 0,
                    clusters = clusters,
                    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE),
                    maxDist = 3e4,
                    cutoff_sig = 2
)

test_that("calculateP2g calculates correlation correctly with cluster argument", {
    expect_equal(corr_matrix,p2g$Correlation, tolerance = 1e-11)
})

colnames(overlap) <- c("RNA", "ATAC")

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

peakMatrix_sce <- SingleCellExperiment(assay=list(counts=peakMatrix), rowRanges=peak.ranges)
geneExpMatrix_sce <- SingleCellExperiment(assay=list(counts=geneExpMatrix), rowRanges=gene.ranges)
suppressWarnings(cellNum <- optimizeMetacellNumber(peakMatrix_sce, geneExpMatrix_sce,
                                  reducedDim=t(as.matrix(geneExpMatrix)), exp_assay="counts",
                      peak_assay="counts", subsample_prop=0.1,
                      n_iter=2, cellNumMin=NULL,
                      cellNumMax=NULL, n_evaluation_points=4,
                      gene_symbol="gene_name"))

min_eval_point <- sqrt(min(20, round(ncol(peakMatrix)/10)))
max_eval_point <- sqrt(min(2000, round(ncol(peakMatrix)/10)))
test_that("optimizeMetacellNumber works correctly", {
    expect_s4_class(cellNum, "CellNumSol")
    expect_equal(length(cellNum@evaluation_points),7)
    expect_equal(cellNum@args$subsample_prop, 0.1)
    expect_length(cellNum@AUC, length(cellNum@evaluation_points))
})
