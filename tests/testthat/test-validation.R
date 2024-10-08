set.seed(1010)
# set up expression and peak matrices
n_cells <- 100
geneExpMatrix <- matrix(abs(rnorm(1e4)), ncol = n_cells,
                        dimnames = list(NULL, paste0("Cell_", seq_len(n_cells))))
geneExpMatrix[geneExpMatrix < 0] <- 0
geneExpMatrix <- as(geneExpMatrix, "CsparseMatrix")

gene.ranges <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr3", "chr4"), 25),
                       ranges = IRanges(start = seq(from = 1, length.out=100, by = 1000),
                                        width = 100))
peak.ranges <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr3", "chr4"), c(40,50,60,70)),
                       ranges = IRanges(start = runif(220)*1e5, width = 100))
peakMatrix <- matrix(sample(c(0,0,0,0,0,1,1,2,3), 220*n_cells, replace = TRUE),
                     ncol = n_cells, dimnames = list(NULL, paste0("Cell_", seq_len(n_cells))))
peak_sce <- SingleCellExperiment(assays = list(counts = peakMatrix))
rowRanges(peak_sce) <- peak.ranges
gene_sce <- SingleCellExperiment(assays = list(logcounts = geneExpMatrix),
                                 colData = data.frame(Clusters = sample(LETTERS[1:3], 100, replace = TRUE)))
gene.ranges$name <- paste0("Gene_", 1:100)
rowRanges(gene_sce) <- gene.ranges


clusters_missing <- gene_sce$Clusters
clusters_missing[sample(length(clusters_missing),1)] <- NA
test_that(".validate_clusters thorws the error when clusters are specified incorrectly", {
  expect_error(.validate_clusters(gene_sce$Clusters[(sample(seq_along(gene_sce$Clusters),
                                                                      length(gene_sce$Clusters)-5))],
                                  gene_sce),
               "'clusters' length should be equal to the number of cells")
    expect_error(.validate_clusters(clusters_missing, gene_sce),
                 "'clusters' object contains NA")
})

gene_sce_incomplete <- gene_sce
rowRanges(gene_sce_incomplete) <- NULL
peak_sce_incomplete <- gene_sce
rowRanges(peak_sce_incomplete) <- NULL

test_that(".validate_input_sce thorws the error when input data is incorrect", {
    expect_error(.validate_input_sce(expMatrix = gene_sce[,0], exp_assay = "logcounts", env = environment()),
                 "SingleCellExperiment with no data")
    expect_error(.validate_input_sce(expMatrix = gene_sce, exp_assay = "logcounts",
                                     peakMatrix = peak_sce[0,], peak_assay = "counts",
                                     env = environment()),
                 "peakMatrix with no data")
    expect_error(.validate_input_sce(expMatrix = gene_sce, exp_assay = "logcounts",
                                     peakMatrix = peak_sce, peak_assay = "no_such_assay",
                                     env = environment()))
    expect_error(.validate_input_sce(expMatrix = gene_sce, exp_assay = "no_such_assay",
                                     peakMatrix = peak_sce, peak_assay = "counts",
                                     env = environment()))
    expect_error(.validate_input_sce(expMatrix = gene_sce[,sample(ncol(gene_sce), ncol(gene_sce)-5)], exp_assay = "logcounts",
                                     peakMatrix = peak_sce, peak_assay = "counts",
                                     env = environment()))
    expect_error(.validate_input_sce(expMatrix = gene_sce_incomplete,
                                     exp_assay = "logcounts",
                                     peakMatrix = peak_sce, peak_assay = "counts",
                                     row.ranges = TRUE, env = environment()))
    expect_error(.validate_input_sce(expMatrix = gene_sce,
                                     exp_assay = "logcounts",
                                     peakMatrix = peak_sce_incomplete, peak_assay = "counts",
                                     row.ranges = TRUE, env = environment()))
})





