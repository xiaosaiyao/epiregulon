% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructRegulon.R
\name{getP2Glinks}
\alias{getP2Glinks}
\title{Compute correlations between ATAC-seq peaks and RNA-seq genes}
\usage{
getP2Glinks(
  archr_path,
  cor_cutoff = 0.5,
  reducedDims = "IterativeLSI",
  useMatrix = "GeneIntegrationMatrix",
  ...
)
}
\arguments{
\item{archr_path}{String indicating the path to a ArchR project that has already performed LSI dimensionality reduction and scRNA-seq integration}

\item{cor_cutoff}{A numeric scalar to indicate the cutoff for correlations between ATAC-seq peaks and RNA-seq genes}

\item{reducedDims}{String specifying which dimensional reduction representation in the ArchR project to use}

\item{useMatrix}{String specifying which data matrix in the ArchR project to use}

\item{...}{other parameters to pass to addPeak2GeneLinks from ArchR package}
}
\value{
A Peak2Gene correlation data frame
}
\description{
Compute correlations between ATAC-seq peaks and RNA-seq genes
}
\author{
Shang-yang Chen
}
