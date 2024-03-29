% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DownstreamAnalysisFunctions.R
\name{clusterLabels}
\alias{clusterLabels}
\title{Cluster labels}
\usage{
clusterLabels(cellWalk, cellTypes, distMethod = "euclidean", plot = FALSE)
}
\arguments{
\item{cellWalk}{a cellWalk object}

\item{cellTypes}{character, vector of labels to use, all labels used by default}

\item{distMethod}{character, method used to compute distance}

\item{plot}{boolean, plot output matrix}
}
\value{
cellWalk object with label clustering stored in "cluster"
}
\description{
\code{clusterLabels()} Computes hierarchical clustering of labels
}
\examples{
data("SampleCellWalkRData")
cellWalk <- clusterLabels(SampleCellWalkRData$cellWalk)

}
