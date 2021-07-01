Interacting with your cellWalk object
================
Pawel F. Przytycki
2021-07-01

Once a cellWalk object has been created, downstream analysis can be performed using an interactive interace. Several additional R packages are required for the interface to function:

``` r
install.packages("shiny")
install.packages("plotly")
install.packages("ggplot2")
install.packages("reshape2")
```

Once required packages are installed, an interface can be launched as follows:

``` r
launchViz(cellWalk)
```

To be able to compute label scores for bulk data, the cellWalk object will need to raw ATAC data to be associated with it:

``` r
cellWalk <- storeMat(cellWalk, ATACMat, peaks)
```

Alternatively, if label scores for bulk data have already been computed (as desribed in the [main vignette](CellWalkR_Vignette.md#Downstream%20Analysis)), they can be added to the cellWalk object before launching the vizualization:

``` r
cellWalk <- storeBulk(cellWalk, bulkPeaks, labelScores)
```
