CellWalkR TensorFlow Vignette
================
Pawel F. Przytycki
2022-01-03

Getting Started
---------------

CellWalkR has built in functions to run on a GPU rather than CPU which allows the code to run more than 15 times faster. This requires having TensorFlow installed along with a high performance GPU. For a simple setup, we recommend using a [prebuilt AWS image](https://www.louisaslett.com/RStudio_AMI/) with TensorFlow already installed.

An instance with GPUs will be necessary for faster analysis. We have tested P2.xlarge for up to 15,000 cells, and P3.2xlarge for up to 20,000 cells. Larger numbers of cells require GPUs with more memory.

TensorFlow
----------

``` r
install.packages("tensorflow")
library("tensorflow")
install_tensorflow()
```

Make sure a GPU available with:

``` r
tf$config$list_physical_devices()
```

With TensorFlow installed and a GPU available, GPU enabled versions of CellWalkR functions can be used. CellWalkR is mostly run in exactly the same way as in the primary [CellWalkR vignette](CellWalkR_Vignette.md) with a few minor changes.

As before, CellWalkR is installed from devtools and we will load a cell-by-peak matrix and label data:

``` r
#install
install.packages("devtools")
devtools::install_github("PFPrzytycki/CellWalkR")

#load library
library(CellWalkR)

#load scATAC-seq data
pathToMat <- system.file("extdata", "SamplePeakMat.mtx", package = "CellWalkR")
ATACMat <- Matrix::readMM(pathToMat)
pathToPeaks <- system.file("extdata", "SamplePeaks.txt", package = "CellWalkR")
peaks <- as(data.table::fread(pathToPeaks, header = FALSE)$V1, "GRanges")

#load labels
pathToLabels <- system.file("extdata", "SampleMarkers1.txt", package = "CellWalkR")
labelGenes <- data.table::fread(pathToLabels)

#compute cell-label edges
regions <- getRegions(geneBody = TRUE, genome = "hg38", names = "Entrez")
ATACGenePeak <- mapPeaksToGenes(labelGenes, ATACMat, peaks, regions)
labelEdges <- computeLabelEdges(labelGenes, ATACMat, ATACGenePeak)
```

TensorFlow can be used for building the network and computing the walk:

#### Computing Cell-Cell Edges

A TensorFlow version of the Jaccard Similarity computation can be used:

``` r
cellEdges <- computeCellSim(ATACMat, method=tensorJaccard)
```

#### Tuning Label Edges

Edge weights and filter weights can be quickly tuned with TensorFlow:

``` r
edgeWeights <- tuneEdgeWeights(cellEdges, 
                              labelEdgesList, 
                              labelEdgeOpts = 10^seq(1,7,1), 
                              sampleDepth = 1000,
                              tensorflow = TRUE)
```

#### Making a cellWalk Object

The CellWalk itself can be computed using TensorFlow:

``` r
cellWalk <- walkCells(cellEdges, 
                     labelEdgesList, 
                     labelEdgeWeights = 1e+07,
                     tensorflow = TRUE)
```

At this point, once the cellWalk object has been created. Downstream analyses are run the same way as in the primary [CellWalkR vignette](CellWalkR_Vignette.md#downstream-analysis).
