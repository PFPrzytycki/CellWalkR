Load Other Data
================
Pawel F. Przytycki
2021-01-26

Cell edges and label edges can be generated directly from a SnapATAC object, an ArchR project, or from Cicero data.

SnapATAC
--------

If working with SnapATAC, the generated SnapATAC object can be used to generate cell and label edges. The user can specify whether to use the "bmat" or "gmat" for analysis.

``` r
cellEdges <- computeCellSim(snap, method="Jaccard")
ATACGenePeak <- mapSnapATACToGenes(labelGenes, snap, whichMat = "bmat", regions)
labelEdges <- computeLabelEdges(labelGenes, snap, ATACGenePeak)
```

These cell and label edges can then be used as in the [CellWalkR vignette](CellWalkR_Vignette.md).

ArchR
-----

If working with ArchR, the generated ArchR project can be used to generate cell and label edges. The user can specify whether to use the "TileMatrix" or "GeneScoreMatrix" for analysis.

``` r
cellEdges <- computeCellSim(ArchRproj, method="Jaccard")
ATACGenePeak <- mapArchRToGenes(labelGenes, ArchRproj, whichMat = "TileMatrix", regions)
labelEdges <- computeLabelEdges(labelGenes, ArchRproj, ATACGenePeak)
```

These cell and label edges can then be used as in the [CellWalkR vignette](CellWalkR_Vignette.md).

Cicero
------

If working with Cicero, two different objects are needed, both of which are part of standard Cicero analysis. First, "cicero\_data" which is a data.frame with peak ids in the first column, cell barcodes in the second column, and read counts in the third, is used fo compute cell edges. Second, the Cicero outputted "cicero\_gene\_activities" Matrix is used to compute label edges.

``` r
cellEdges <- computeCellSim(cicero_data, method="Jaccard")
ATACGenePeak <- mapCiceroToGenes(labelGenes, cicero_gene_activities,)
labelEdges <- computeLabelEdges(labelGenes, cicero_gene_activities, ATACGenePeak)
```

These cell and label edges can then be used as in the [CellWalkR vignette](CellWalkR_Vignette.md).
