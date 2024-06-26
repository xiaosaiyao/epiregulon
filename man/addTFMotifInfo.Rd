% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructRegulon.R
\name{addTFMotifInfo}
\alias{addTFMotifInfo}
\title{Add TF binding motif occupancy information to the peak2gene object}
\usage{
addTFMotifInfo(p2g, grl, peakMatrix = NULL)
}
\arguments{
\item{p2g}{A Peak2Gene dataframe created by ArchR or getP2Glinks() function}

\item{grl}{GRangeList object containing reference TF binding information. We recommend retrieving \code{grl} from \code{getTFMotifInfo} which
contains TF occupancy data derived from public and ENCODE ChIP-seq peaks. Alternatively, if the users would like to provide a GRangeList
of motif annotations. This can be derived using \code{motifmatchr::matchMotifs}. See details}

\item{peakMatrix}{A matrix of scATAC-seq peak regions with peak ids as rows}
}
\value{
A data frame containing overlapping ids of scATAC-seq peak regions and reference TF binding regions
}
\description{
Add TF binding motif occupancy information to the peak2gene object
}
\details{
This function annotates each regulatory element with possible transcription factors. We can either provide a GRangeList of known ChIP-seq
binding sites (TF occupancy) or positions of TF motifs (TF motifs). While public ChIP-seq data may not fully align with the ground truth TF occupancy in users' data
(due to technical challenges of ChIP-seq or cell type nature of TF occupancy), it does offer a few important advantages over TF motif information:
\enumerate{
\item{TF occupancy allows co-activators to be included. Co-activators are chromatin modifiers that do not directly bind to DNA but nonetheless play an important role
in gene regulation}
\item{TF occupancy can distinguish between members of the same class that may share similar motifs but that may have drastically different binding sites}
}
If multiple ChIP-seq are available for the same TF, we merge the ChIP-seq data to represent an universal set of possible binding sites. The predicted TF
occupancy is further refined by \code{\link{pruneRegulon}}.

If the users prefer to use TF motifs instead of TF occupancy, the users can create a GRangeList of motif annotation using \code{motifmatchr::matchMotifs}.
Here, we demonstrate how to annotate peaks with cisbp motif database

\if{html}{\out{<div class="sourceCode">}}\preformatted{library(motifmatchr)
library(chromVARmotifs)
data("human_pwms_v1")
peaks <- GRanges(seqnames = c("chr1","chr2","chr2"),
                               ranges = IRanges(start = c(76585873,42772928, 100183786),
                                                         width = 500))
eh <- AnnotationHub::query(ExperimentHub::ExperimentHub(),
pattern = c("scMultiome", "TF motifs", "human"))
pwms <- readRDS(eh[[eh$ah_id]]))
grl <- matchMotifs(pwms, peaks, genome = "hg38", out = "positions")
retain only TF symbols. TF symbols need to be consistent with gene names in regulon
names(grl) <- sapply(strsplit(names(grl), "_"), "[",3)
}\if{html}{\out{</div>}}
}
\examples{
set.seed(1)
# create a mock peak-to-gene matrix
p2g <- data.frame(idxATAC = c(rep(1,5), rep(2,5)), Chrom = 'chr1', idxRNA = 1:10,
Gene = paste0('Gene_',1:10),Correlation = runif(10, 0,1))

# create mock a GRanges list of TF binding sites
grl <- GRangesList('TF1' = GRanges(seqnames = 'chr1',
ranges = IRanges(start = c(50,1050), width = 100)),
'TF2' = GRanges(seqnames = 'chr1',
ranges = IRanges(start = c(1050), width = 100))
)

# create a mock singleCellExperiment object for peak matrix
peak_gr <- GRanges(seqnames = 'chr1',
             ranges = IRanges(start = seq(from = 1, to = 10000, by = 1000),
             width = 100))
peak_counts <- matrix(sample(x = 0:4, size = 100*length(peak_gr), replace = TRUE),
nrow = length(peak_gr), ncol = 100)
peak_sce <- SingleCellExperiment(list(counts = peak_counts))
rowRanges(peak_sce) <- peak_gr
rownames(peak_sce) <- paste0('peak',1:10)

# create overlaps between p2g matrix, TF binding sites and peak matrix
overlap <- addTFMotifInfo(p2g, grl, peakMatrix = peak_sce)
utils::head(overlap)
}
\author{
Xiaosai Yao, Shang-yang Chen
}
