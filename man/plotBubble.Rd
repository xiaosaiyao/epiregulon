% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plotBubble}
\alias{plotBubble}
\title{Generate bubble plots of relative activities across cluster/group labels for a list of TFs}
\usage{
plotBubble(
  activity_matrix,
  tf,
  clusters,
  bubblesize = c("FDR", "summary.logFC")
)
}
\arguments{
\item{activity_matrix}{A matrix of TF activities inferred from calculateActivity}

\item{tf}{A character vector indicating the names of the transcription factors to be plotted}

\item{clusters}{A character or integer vector of cluster or group labels for single cells}

\item{bubblesize}{String indicating the variable from findDifferentialActivity output to scale size of bubbles
by either \code{FDR} or \code{summary.logFC}. Default is \code{FDR}.}
}
\value{
A ggplot object
}
\description{
Generate bubble plots of relative activities across cluster/group labels for a list of TFs
}
\examples{
example_sce <- scuttle::mockSCE()
example_sce <- scuttle::logNormCounts(example_sce)
example_sce$cluster <- sample(LETTERS[1:5], ncol(example_sce), replace = TRUE)
plotBubble(activity_matrix = logcounts(example_sce),
tf = c("Gene_0001","Gene_0002"),  clusters = example_sce$cluster)
}
\author{
Shang-yang Chen
}
