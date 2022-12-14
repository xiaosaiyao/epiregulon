% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructRegulon.R
\name{addTFMotifInfo}
\alias{addTFMotifInfo}
\title{Add TF binding motif occupancy information to the peak2gene object}
\usage{
addTFMotifInfo(p2g, grl, peakMatrix = NULL, archR_project_path = NULL)
}
\arguments{
\item{p2g}{A Peak2Gene dataframe created by ArchR or getP2Glinks() function}

\item{grl}{GRangeList object containing reference TF binding information}

\item{peakMatrix}{A matrix of scATAC-seq peak regions with peak ids as rows}

\item{archR_project_path}{Path to an ArchR project that have performed LSI dimensionality reduction and scRNA-seq integration}
}
\value{
A dataframe containing overlapping ids of scATAC-seq peak regions and reference TF binding regions
}
\description{
Add TF binding motif occupancy information to the peak2gene object
}
\examples{
set.seed(1)
# create a mock peak-to-gene matrix
p2g <- data.frame(idxATAC = c(rep(1,5), rep(2,5)), Chrom = "chr1", idxRNA = 1:10,
Gene = paste0("Gene_",1:10),Correlation = runif(10, 0,1))

# create mock a GRanges list of TF binding sites
grl <- GRangesList("TF1" = GRanges(seqnames = "chr1",
ranges = IRanges(start = c(50,1050), width = 100)),
"TF2" = GRanges(seqnames = "chr1",
ranges = IRanges(start = c(1050), width = 100))
)

# create a mock singleCellExperiment object for peak matrix
peak_gr <- GRanges(seqnames = "chr1",
             ranges = IRanges(start = seq(from = 1, to = 10000, by = 1000),
             width = 100))
peak_counts <- matrix(sample(x = 0:4, size = 100*length(peak_gr), replace = TRUE),
nrow = length(peak_gr), ncol = 100)
peak_sce <- SingleCellExperiment(list(counts = peak_counts))
rowRanges(peak_sce) <- peak_gr
rownames(peak_sce) <- paste0("peak",1:10)

# create overlaps between p2g matrix, TF binding sites and peak matrix
overlap <- addTFMotifInfo(p2g, grl, peakMatrix = peak_sce)
utils::head(overlap)
}
\author{
Xiaosai Yao, Shang-yang Chen
}
