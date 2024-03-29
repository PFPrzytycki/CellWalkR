% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataParsingFunctions.R
\name{mapCiceroToGenes}
\alias{mapCiceroToGenes}
\title{Map peaks to genes based on Cicero coaccessibility}
\usage{
mapCiceroToGenes(labelGenes, cicero_gene_activities)
}
\arguments{
\item{labelGenes}{data.frame with genes of interest in first column and corresponding
labels in second column}

\item{cicero_gene_activities}{gene activity matrix}
}
\value{
list of genes and corresponding peaks by label
}
\description{
\code{mapCiceroToGenes} Generates a mapping from genes to peaks per label
}
