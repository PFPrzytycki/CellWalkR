#' Cell-by-peak matrix
#'
#' A cell-by-peak matrix generated from a scATAC-seq experiment
#'
#' @format A matrix with 1723 rows and 100000 columns
"ATACMat"

#' Peaks for ATACMat
#'
#' A GRanges object of peaks in ATACMat
#'
#' @format A GRanges object with 100000 ranges
"peaks"

#' Label genes
#'
#' A data.frame with labeling data
#'
#' @format A data frame with 12410 rows and 3 columns
"labelGenes"

#' Label genes B
#'
#' A data.frame with labeling data
#'
#' @format A data frame with 1838 rows and 3 columns
"labelGenesB"

#' Peaks for filter
#'
#' A GRanges object of peaks for filter
#'
#' @format A GRanges object with 21331 ranges
"filter"

#' Enhancer regions
#'
#' A GRanges object of enhancer regions
#'
#' @format A GRanges object with 20634 ranges
"sampleEnhancers"
