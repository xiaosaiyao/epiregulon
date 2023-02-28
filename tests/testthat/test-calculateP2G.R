set.seed(1010)
n_cells <- 100
geneExpMatrix <- matrix(abs(rnorm(1e4)), ncol = n_cells,
                        dimnames = list(NULL, paste0("Cell_", seq_len(n_cells))))
geneExpMatrix[geneExpMatrix < 0] <- 0
geneExpMatrix <- as(geneExpMatrix, "sparseMatrix")

gene.ranges <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr3", "chr4"), 25),
                   ranges = IRanges(start = seq(from = 1, length.out=100, by = 1000),
                                    width = 100))
peak.ranges <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr3", "chr4"), c(40,50,60,70)),
                       ranges = IRanges(start = runif(220)*1e5, width = 100))
peakMatrix <- matrix(sample(c(0,0,0,0,0,1,1,2,3), 220*n_cells, replace = TRUE),
                     ncol = n_cells, dimnames = list(NULL, paste0("Cell_", seq_len(n_cells))))
peakMatrix[c(80,90),] <- 0
geneExpMatrix[c(10,30),] <- 0
peak_sce <- SingleCellExperiment(assays = list(counts = peakMatrix))
rowRanges(peak_sce) <- peak.ranges
reducedDimMatrix <- matrix(runif(n_cells*40), nrow = n_cells)
gene_sce <- SingleCellExperiment(assays = list(logcounts = geneExpMatrix))
gene.ranges$name <- paste0("Gene_", 1:100)
rowRanges(gene_sce) <- gene.ranges

# set seed to assure reproducibility with scran::clusterCells
set.seed(1100)
clusters <- kmeans(reducedDimMatrix, 10)$cluster

geneExpMatrix.avg <- t(apply(geneExpMatrix, 1, function(x) tapply(x,clusters, mean)))
# remove genes that are equal to 0
non.zero.genes <- which(rowSums(geneExpMatrix.avg) != 0)
peakMatrix.avg <- t(apply(peakMatrix, 1, function(x) tapply(x,clusters, mean)))
gene.start <- promoters(gene.ranges[non.zero.genes,])
overlap <- S4Vectors::DataFrame(findOverlaps(resize(gene.start, 5000, "center"),
                                       peak.ranges))

# save row numbers to later set them as row names for compatibility with
# calculateP2G output
overlap$row_numb <- seq_len(nrow(overlap))

overlap$Correlation <- mapply(cor, asplit(geneExpMatrix.avg[non.zero.genes,][overlap[,1],],1),
                              asplit(peakMatrix.avg[overlap[,2],],1))

overlap$Correlation <- mapply(cor, as.data.frame(t(geneExpMatrix.avg[non.zero.genes,][overlap[,1],])),
                              as.data.frame(t(peakMatrix.avg[overlap[,2],])))

overlap$distance <- distance(gene.start[overlap[,1], ], peak.ranges[overlap[,2], ])

overlap$TStat <- (overlap$Correlation /
              sqrt((pmax(1 - overlap$Correlation ^ 2, 0.00000000000000001, na.rm = TRUE))
                  / (ncol(peakMatrix.avg) - 2))) #T-statistic P-value

overlap$Pval <- 2 * stats::pt(-abs(overlap$TStat), ncol(peakMatrix.avg) - 2)
overlap$FDR <- stats::p.adjust(overlap$Pval, method = "fdr")
overlap$target <-  gene.ranges[non.zero.genes,][overlap[,1],]$name
overlap$chr <- as.character(seqnames(peak.ranges[overlap[,2], ]))
overlap$start <- GenomicRanges::start(peak.ranges[overlap[,2], ])
overlap$end <- GenomicRanges::end(peak.ranges[overlap[,2], ])
overlap <- overlap[order(overlap[,2], overlap[,1]),]
overlap <- overlap[overlap$Correlation > 0.5,]
overlap <- as.data.frame(overlap)
attr(overlap, "row.names") <- overlap$row_numb

set.seed(1100)
test_that("calculateP2G works correctly", {
  P2G <- calculateP2G(peakMatrix = peak_sce,
                      peak_assay = "counts",
                      expMatrix = gene_sce,
                      exp_assay = "logcounts",
                      reducedDim = reducedDimMatrix,
                      cellNum = 10,
                      maxDist = 5000)
  expect_equal(P2G[,c("Correlation", "FDR")], overlap[,c("Correlation", "FDR")], tolerance = 1e-10)
  expect_equal(P2G[,c("distance","target","chr","start","end")], overlap[,c("distance","target","chr","start","end")])
})


