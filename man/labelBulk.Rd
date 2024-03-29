% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DownstreamAnalysisFunctions.R
\name{labelBulk}
\alias{labelBulk}
\title{Label bulk data}
\usage{
labelBulk(
  cellWalk,
  bulkPeaks,
  ATACMat,
  peaks,
  extendRegion,
  extendDistance,
  cellTypes,
  allScores = FALSE,
  parallel = FALSE,
  numCores = 1
)
}
\arguments{
\item{cellWalk}{a cellWalk object}

\item{bulkPeaks}{GRanges of peaks in bulk data or GRangesList of sets of peaks}

\item{ATACMat}{cell-by-peak matrix}

\item{peaks}{GRanges of peaks in ATACMat}

\item{extendRegion}{GRanges defining where to extend mapping to consider peaks in a larger region (e.g. in LD) with bulk data}

\item{extendDistance}{numeric maximum distance to extend region by (if region is missing, just distance is used, if distance is missing whole region is used)}

\item{cellTypes}{character, vector of labels to use, all labels used by default}

\item{allScores}{return full table of scores}

\item{parallel}{execute in parallel}

\item{numCores}{number of cores to use for parallel execution}
}
\value{
labels for each region in bulk data
}
\description{
\code{labelBulk()} Determines labels for bulk data by mapping them via
the calculated information matrix
}
\examples{
data("SampleCellWalkRData")
labelBulk(SampleCellWalkRData$cellWalk,
          SampleCellWalkRData$sampleEnhancers,
          SampleCellWalkRData$ATACMat,
          SampleCellWalkRData$peaks)

}
