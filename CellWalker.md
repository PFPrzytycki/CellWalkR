<img src="examples/CellWalkR_Vignette_files/figure-markdown_github/cellwalker_icon.png" id="id" class="class" width="50" height="50" /> CellWalkR
================

## About

CellWalkR is an R package that integrates single-cell open chromatin
(scATAC-seq) data with cell type labels and bulk epigenetic data to
identify cell type-specific regulatory regions. A GPU implementation and
downsampling strategies enable thousands of cells to be processed in
seconds. CellWalkR’s user-friendly interface provides interactive
analysis and visualization of cell labels and regulatory region
mappings.

## Installation

Install CellWalkR for R using devtools as follows:

``` r
$ R
> install.packages("devtools")
> devtools::install_github("PFPrzytycki/CellWalkR")
```

## Usage

For a guide to using CellWalkR, see the provided
[vignette](examples/CellWalkR_Vignette.md), which covers the following:

1.  [Data
    Pre-processing](examples/CellWalkR_Vignette.md#data-pre-processing)
2.  [Getting Started with
    CellWalkR](examples/CellWalkR_Vignette.md#getting-started-with-cellwalkr)
    1.  [Loading scATAC-seq
        Data](examples/CellWalkR_Vignette.md#loading-scatac-seq-data)
    2.  [Defining Label
        Nodes](examples/CellWalkR_Vignette.md#defining-label-nodes)
3.  [Building a
    Network](examples/CellWalkR_Vignette.md#building-a-network)
    1.  [Computing Cell-Cell
        Edges](examples/CellWalkR_Vignette.md#computing-cell-cell-edges)
    2.  [Computing Label-Cell
        Edges](examples/CellWalkR_Vignette.md#computing-label-cell-edges)
4.  [Tuning Label
    Edges](examples/CellWalkR_Vignette.md#tuning-label-edges)
5.  [Making a cellWalk
    Object](examples/CellWalkR_Vignette.md#making-a-cellwalk-object)
6.  [Adding Filters](examples/CellWalkR_Vignette.md#adding-filters)
7.  [Downstream
    Analysis](examples/CellWalkR_Vignette.md#downstream-analysis)
    1.  [Cell Labels](examples/CellWalkR_Vignette.md#cell-labels)
    2.  [Confusion
        Matrix](examples/CellWalkR_Vignette.md#confusion-matrix)
    3.  [Hierarchical Clustering of
        Labels](examples/CellWalkR_Vignette.md#hierarchical-clustering-of-labels)
    4.  [Plotting Cells](examples/CellWalkR_Vignette.md#plotting-cells)
    5.  [Bulk Data
        Mapping](examples/CellWalkR_Vignette.md#bulk-data-mapping)
8.  [Interactive
    Visualizaiton](examples/CellWalkR_Vignette.md#interactive-visualzation)
9.  [Adding a Second Set of
    Labels](examples/CellWalkR_Vignette.md#adding-a-second-set-of-labels)
10. [Detecting
    Doublets](examples/CellWalkR_Vignette.md#detecting-doublets)

If you use CellWalkR please cite:

1.  Przytycki, P.F., Pollard, K.S. “CellWalkR: An R Package for
    integrating and visualizing single-cell and bulk data to resolve
    regulatory elements.” *Bioinformatics* (2022).
    <https://doi.org/10.1093/bioinformatics/btac150>

2.  Przytycki, P.F., Pollard, K.S. “CellWalker integrates single-cell
    and bulk data to resolve regulatory elements across cell types in
    complex tissues.” *Genome Biology* (2021).
    <https://doi.org/10.1186/s13059-021-02279-1>

## AWS + TensorFlow

CellWalkR can also be run on AWS which vastly simplifies the process of
running on GPUs using TensorFlow. Using GPUs allows the code to run more
than 15 times faster. For a guide to running CellWalkR on AWS using GPUs
see this [vignette](examples/CellWalkR_TensorFlow_Vignette.md).
