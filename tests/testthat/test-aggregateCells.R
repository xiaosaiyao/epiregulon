library(SingleCellExperiment)

# construct sce
ngene <- 5
ncell <- 10
ncolData <- 2
counts <- matrix(floor(runif(ngene*ncell, min=0, max=10)), nrow=ngene, ncol=ncell)
rownames(counts) <- paste0("gene",1:ngene)
colnames(counts) <- paste0("cell",1:ncell)

colData <- matrix(NA, ncell, ncolData)
rownames(colData) <- paste0("cell", 1:ncell)
colnames(colData) <- c("cluster", "batch")
colData <- data.frame(colData)
colData$cluster <- paste0("celltype", sample(1:3, ncell, replace=TRUE))
colData$batch[colData$cluster %in% c("celltype1", "celltype2")] <- 1
colData$batch[colData$cluster %in% c("celltype3")] <- 2
sce <- SingleCellExperiment(list(counts=counts), colData=colData)

sce_bulk <- aggregateAcrossCellsFast(sce, clusters=colData$cluster, fun_name="sum", aggregateColData=TRUE)

# check counts
aggregated_counts <- matrix(data=NA, nrow=ngene, ncol=length(unique(colData$cluster)))
rownames(aggregated_counts) <- paste0("gene",1:ngene)
colnames(aggregated_counts) <- unique(colData$cluster)
for (cluster in unique(colData$cluster)){
  cell_select <- which(colData(sce)$cluster == cluster)
  aggregated_counts[, cluster] <- rowSums(counts[, cell_select, drop=FALSE])
}

test_that("count matrix of aggregateAcrossCellsFast works ", {
  expect_equal(counts(sce_bulk), aggregated_counts[rownames(sce_bulk), colnames(sce_bulk)])
})

# test colData
colData_aggregated <- matrix(data=NA, nrow=length(unique(colData$cluster)), ncol= ncolData+2)
rownames(colData_aggregated) <- unique(colData$cluster) 
colnames(colData_aggregated) <- c("idx", "ncells", "cluster","batch")
colData_aggregated <- DataFrame(colData_aggregated)
colData_aggregated$idx <- unique(colData$cluster)
colData_aggregated$cluster <- unique(colData$cluster)
colData_aggregated$batch[colData_aggregated$cluster %in% c("celltype1", "celltype2")] <- 1
colData_aggregated$batch[colData_aggregated$cluster %in% c("celltype3")] <- 2
freq <- table(sce$cluster)
colData_aggregated[names(freq), "ncells"] <- as.vector(freq)

test_that("colData of aggregateAcrossCellsFast works ", {
  expect_equal(colData(sce_bulk), colData_aggregated[rownames(colData(sce_bulk)), colnames(colData(sce_bulk))])
})

