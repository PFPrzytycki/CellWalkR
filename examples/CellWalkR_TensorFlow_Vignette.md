CellWalkR TensorFlow Vignette
================
Pawel F. Przytycki
2021-02-09

Getting Started
---------------

CelLWalkR has built in functions to run on a GPU rather than CPU which allows the code to run more than 15 times faster. This requires having TensorFlow installed along with a high performance GPU. For a simple setup, we recommend using a [prebuilt AWS image](https://www.louisaslett.com/RStudio_AMI/) with TensorFlow already installed.

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

1.  A TensorFlow version of the Jaccard Similarity computation can be used:

    ``` r
    cellEdges <- computeCellSim(ATACMat, method=tensorJaccard)
    ```

2.  Edge weights and filter weights:

    ``` r
    edgeWeights <- tuneEdgeWeights(cellEdges, 
                              labelEdgesList, 
                              labelEdgeOpts = 10^seq(1,7,1), 
                              sampleDepth = 1000,
                              tensorflow = TRUE)
    ```

3.  Walking cells:

    ``` r
    cellWalk <- walkCells(cellEdges, 
                     labelEdgesList, 
                     labelEdgeWeights = 1e+07,
                     tensorflow = TRUE)
    ```
