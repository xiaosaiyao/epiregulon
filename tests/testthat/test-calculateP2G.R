set.seed(1010)
# set up expression and peak matrices
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



### test pseudobulk

sce_combined <- combineSCE(gene_sce, "logcounts", peak_sce, "counts", reducedDimMatrix, "reducedDim")
sce_grouped <- applySCE(sce_combined,
                        scuttle::aggregateAcrossCells,
                        ids = clusters,
                        statistics = "mean")


# extract gene expression and peak matrix
expGroupMatrix <- assay(sce_grouped, "counts")
peakGroupMatrix <- assay(altExp(sce_grouped), "counts")


geneExpMatrix.avg <- t(apply(geneExpMatrix, 1, function(x) tapply(x,clusters, mean)))
peakMatrix.avg <- t(apply(peakMatrix, 1, function(x) tapply(x,clusters, mean)))


test_that("pseudobulk works correctly", {
  expect_equal(expGroupMatrix, geneExpMatrix.avg)
  expect_equal(peakGroupMatrix, peakMatrix.avg)
})


### test overlap


# remove genes that are equal to 0
non.zero.genes <- which(rowSums(geneExpMatrix.avg) != 0)
non.zero.peaks <- which(rowSums(peakMatrix.avg) != 0)

gene.start <- promoters(gene.ranges[non.zero.genes,])
peak.ranges <- peak.ranges[non.zero.peaks,]
overlap <- S4Vectors::DataFrame(findOverlaps(resize(gene.start, 5000, "center"),
                                       peak.ranges))

# remove genes and peaks that are equal to 0
sce_grouped <- sce_grouped[which(rowSums(assay(sce_grouped)) != 0),]
altExp(sce_grouped) <- altExp(sce_grouped)[which(rowSums(assay(altExp(sce_grouped), "counts")) != 0),]

# get gene information
geneSet <- rowRanges(sce_grouped)
geneStart <- promoters(geneSet)

# get peak range information
peakSet <- rowRanges(altExp(sce_grouped))

# find overlap after resizing
o <- S4Vectors::DataFrame(findOverlaps(resize(geneStart, 5000, "center"),
                                       peakSet,
                                       ignore.strand = TRUE))



test_that("overlap works correctly", {
  expect_equal(o, overlap)
})


### test calculateP2G
overlap$Correlation <- mapply(cor, asplit(geneExpMatrix.avg[non.zero.genes,][overlap[,1],],1),
                              asplit(peakMatrix.avg[non.zero.peaks,][overlap[,2],],1))

overlap$Correlation <- mapply(cor, as.data.frame(t(geneExpMatrix.avg[non.zero.genes,][overlap[,1],])),
                              as.data.frame(t(peakMatrix.avg[non.zero.peaks,][overlap[,2],])))

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


set.seed(1100)
P2G <- calculateP2G(peakMatrix = peak_sce,
                    peak_assay = "counts",
                    expMatrix = gene_sce,
                    exp_assay = "logcounts",
                    reducedDim = reducedDimMatrix,
                    cellNum = 10,
                    maxDist = 5000,
                    cor_cutoff = 0.5)

test_that("calculateP2G works correctly", {

  expect_equal(as.vector(P2G$Correlation), as.vector(overlap$Correlation), tolerance = 1e-10)
  expect_equal(data.frame(P2G[,c("distance","target","chr","start","end")]), overlap[,c("distance","target","chr","start","end")])
})


