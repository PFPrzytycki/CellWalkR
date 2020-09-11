CellWalkR Vignette
================
Pawel F. Przytycki
2020-09-11

## Getting Started

CellWalkR is an R package implementing the CellWalker method for
combining scATAC-seq data with labels and other epigenetic data. This
vignette shows an example of running CellWalkR on a small set of
scATAC-seq data to generate a cellWalk object which can then be used to
assign labels to cells as well as cell-type specific labels to bulk
data.

``` r
library(CellWalkR)
```

First, load scATAC-seq data in the form of a cell-by-peak matrix and
load the corresponding peaks into a GRanges object. If working with a
SnapATAC object, a cell-by-peak matrix should be stored in @pmat.
**Please note:** all data needed to recreate this vignette locally is
already loaded so read commands can be skipped. They are included as an
example for how new data could be loaded.

``` r
ATACMat = Matrix::readMM("extdata/SamplePeakMat.mtx")
peaks = as(data.table::fread("extdata/SamplePeaks.txt", header = FALSE)$V1, "GRanges")
```

Next, we compute cell-to-cell similarity in order to build edges in the
cell-to-cell portion of the graph. If working with a SnapATAC object,
Jaccard similarity is stored in @jmat.

``` r
cellEdges = computeCellSim(ATACMat, method="Jaccard")
```

cellEdges is a cell-by-cell matrix of cell-to-cell similarity. Any
matrix of cell similarity can be used.

``` r
head(cellEdges, n=c(5,5))
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>                                                            
#> [1,] 1.00000000 0.11338151 0.15001705 0.09244314 0.09813385
#> [2,] 0.11338151 1.00000000 0.14950372 0.08130564 0.09035017
#> [3,] 0.15001705 0.14950372 1.00000000 0.08960442 0.11350499
#> [4,] 0.09244314 0.08130564 0.08960442 1.00000000 0.06237177
#> [5,] 0.09813385 0.09035017 0.11350499 0.06237177 1.00000000
```

In order to generate label-to-cell edges, we need to define which
genomic regions correspond to which genes. These could be promoters,
gene bodies, or any other definition. In this example we will use full
gene bodies.

``` r
regions = getRegions(geneBody = TRUE, genome = "hg38", names = "Entrez")
```

regions is a GRanges object with a gene\_id field. This gene\_id field
needs to match the genes (or other identifiers) in the labeling data.

``` r
head(regions)
#> GRanges object with 6 ranges and 1 metadata column:
#>             seqnames            ranges strand |     gene_id
#>                <Rle>         <IRanges>  <Rle> | <character>
#>           1    chr19 58362552-58364751      - |           1
#>          10     chr8 18389282-18391481      + |          10
#>         100    chr20 44652034-44654233      - |         100
#>        1000    chr18 28176931-28179130      - |        1000
#>   100009613    chr11 70075234-70077433      - |   100009613
#>   100009667    chr10 68010663-68012862      - |   100009667
#>   -------
#>   seqinfo: 595 sequences (1 circular) from hg38 genome
```

Now we can load our first set of labeling data. If no labeling data is
available, can run findMarkers() on a set of scRNA-seq data.

``` r
labelGenes = data.table::fread("extdata/SampleMarkers1.txt")
```

The labeling data should consist of at least two columns, gene names (or
other identifiers that match regions) and associated labels, with an
optional third column for log-fold change in expression of that gene for
that label.

``` r
head(labelGenes)
#>   entrez  cluster   avg_diff
#> 1  10299 RG-early -1.2297890
#> 2   6167 RG-early  0.2546596
#> 3  11168 RG-early  0.2570446
#> 4   8760 RG-early -0.2578798
#> 5   8503 RG-early  0.2613031
#> 6  10208 RG-early  0.2618003
```

We then need to map between this data and the peaks in the scATAC-seq
data.

``` r
ATACGenePeak = mapPeaksToGenes(labelGenes, ATACMat, peaks, regions)
```

With this mapping we can compute label-to-cell edges, a matrix where the
number of rows is the number of labels and the number of columns is the
number of cells.

``` r
labelEdges = computeLabelEdges(labelGenes, ATACMat, ATACGenePeak)
head(labelEdges)
#>      RG-early        oRG         tRG vRG     RG-div1    RG-div2
#> [1,]        0 0.02591189 0.025962569   0 0.040257391 0.01894717
#> [2,]        0 0.01897995 0.017158561   0 0.022360435 0.02047899
#> [3,]        0 0.02585844 0.015839112   0 0.019537577 0.02276763
#> [4,]        0 0.01583931 0.003441275   0 0.002901203 0.01882845
#> [5,]        0 0.01466549 0.020814822   0 0.019853461 0.03259348
#> [6,]        0 0.02243765 0.011575348   0 0.006914287 0.01302038
```

## Tuning Label Edges

Although we now have cell-to-cell edges and label-to-cell edges, we
don’t know how to correctly weight the two relative to each other. The
tuneEdgeWeights method will run CellWalker across a range of possible
parameters and compute cell homogeneity for each. We make a list of
labelEdges because there can be many of them, and sample down to 1000
cells for faster computation.

``` r
labelEdgesList = list(labelEdges)
edgeWeights = tuneEdgeWeights(cellEdges, 
                              labelEdgesList, 
                              labelEdgeOpts=10^seq(1,7,1), 
                              sampleDepth=1000)
```

We can see which parameter had the highest cell homogeneity:

``` r
head(edgeWeights[order(edgeWeights$cellHomogeneity, decreasing = TRUE),])
#>    Var1 cellHomogeneity
#> 7 1e+07      0.08080729
#> 6 1e+06      0.05416757
#> 5 1e+05      0.05411512
#> 4 1e+04      0.04437068
#> 3 1e+03     -0.01519351
#> 2 1e+02     -0.23014493
```

And can generate a cellWalk object with this parameter. This object
stores the final influence matrix and can be used for downstream
analysis.

``` r
cellWalk = walkCells(cellEdges, 
                     labelEdgesList, 
                     labelEdgeWeights = 1e+07)
```

## Adding Filters

We may have some bulk epigenetic data that can help filter down which
peaks are relevant to our analysis. We can tune weights on each filter
to determine how signficant it is to our data. For our example we have
H3K4me3 data which indicates active promoters. Thus we apply this filter
permissively (setting filterOut=FALSE) and at the whole gene level
(filterGene=TRUE) rather than just to overlaping peaks.

``` r
filter = data.table::fread("extdata/SampleFilter.bed")
filter = GRanges(filter$V1, IRanges(filter$V2, filter$V3))
```

``` r
filters = list(filter)
labelGenesList = list(labelGenes)
filterWeights = tuneFilterWeights(cellEdges, 
                                  labelGenesList, 
                                  labelEdgesList, 
                                  labelEdgeWeights = 1e+07,
                                  ATACMat,
                                  ATACGenePeak,
                                  filters = filters,
                                  filterOut = c(FALSE),
                                  filterGene = c(TRUE),
                                  regions=regions, 
                                  sampleDepth=1000)
filterWeights
#>   Var1 cellHomogeneity
#> 1    0      0.06987248
#> 2    1      0.08754723
```

We see that adding this filter improves performance. We can make a new
cellWalk object using this filter:

``` r
labelEdges = computeLabelEdges(labelGenes, 
                               ATACMat, 
                               ATACGenePeak,
                               filters = filters, 
                               filterWeights = c(1),
                               filterOut = c(FALSE),
                               filterGene = c(TRUE),
                               regions = regions)
labelEdgesList = list(labelEdges)
cellWalk = walkCells(cellEdges, 
                     labelEdgesList, 
                     labelEdgeWeights = 1e+07)
```

## Downstream Analysis

Once we have created a cellWalk object, we can use it for downstream
analysis. Most directly, we can look at what labels are the most
strongly linked to each cell. This is based on the amount of
label-to-cell influence in the cell walk. Cell labeling can be used for
numerous further downstream analyses such as cell-type specific peak
calling.

``` r
head(cellWalk$cellLabels)
#> [1] "oRG"     "tRG"     "tRG"     "oRG"     "tRG"     "RG-div1"
```

We can dig further into cell labeling by examining how often labels are
confused for each other.

``` r
cellWalk = findUncertainLabels(cellWalk, plot=TRUE)
```

![](CellWalkR_Vignette_files/figure-gfm/downstream-uncertainty-1.png)<!-- -->

We can also directly examine label similarity by considering
label-to-label influence.

``` r
cellWalk = clusterLabels(cellWalk, plot=TRUE)
```

![](CellWalkR_Vignette_files/figure-gfm/downstream-labelLabels-1.png)<!-- -->

    #> NULL

A very powerful use for the cell walk is mapping data to labels via
cell-to-label influence. For example, we can map enhancers to cell
types.

``` r
sampleEnhancers = data.table::fread("extdata/sampleEnhancers.bed")
sampleEnhancers = GRanges(sampleEnhancers$V1, IRanges(sampleEnhancers$V2, sampleEnhancers$V3))
```

``` r
mappedLabel = labelBulk(cellWalk, 
                        sampleEnhancers[1:10], 
                        ATACMat, 
                        peaks)
mappedLabel
#>  [1] "vRG"      "RG-early" "RG-early" "RG-early" "RG-early" "oRG"     
#>  [7] "RG-early" "RG-early" "RG-early" "RG-early"
```

## Adding a Second Set of Labels

CellWalkR can be run on an arbitrary number of sets of labels and
filters, each with it’s own weight. Filters can selectively be applied
to some sets of labels and not others. Here for example, we will add a
second set of labels to which the above filter does not apply.

``` r
labelGenesB = data.table::fread("extdata/SampleMarkers2.txt")
```

``` r
ATACGenePeakB = mapPeaksToGenes(labelGenesB, ATACMat, peaks, regions)
labelEdgesB = computeLabelEdges(labelGenesB, ATACMat, ATACGenePeakB)
```

Now simply tune edge weights as before with a list of all label edges.

``` r
labelEdgesListB = list(labelEdges, labelEdgesB)
edgeWeightsB = tuneEdgeWeights(cellEdges, 
                               labelEdgesListB, 
                               labelEdgeOpts = 10^seq(4,7,1),
                               sampleDepth = 1000)
head(edgeWeightsB[order(edgeWeightsB$cellHomogeneity, decreasing = TRUE),])
#>     Var1  Var2 cellHomogeneity
#> 9  1e+04 1e+06       0.5011993
#> 13 1e+04 1e+07       0.4894711
#> 14 1e+05 1e+07       0.4518762
#> 10 1e+05 1e+06       0.4264239
#> 15 1e+06 1e+07       0.4203451
#> 16 1e+07 1e+07       0.3927239
```

We can then compute a new cell walk using the list of edges and a vector
of optimal weights.

``` r
cellWalkB = walkCells(cellEdges, 
                      labelEdgesListB, 
                      labelEdgeWeights = c(1e+04, 1e+06))
```
