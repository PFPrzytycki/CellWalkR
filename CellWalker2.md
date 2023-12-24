<img src="examples/CellWalkR_Vignette_files/figure-markdown_github/cellwalker2_icon.png" id="id" class="class" width="100" height="50" /> CellWalker2
================

## About

CellWalker2 is integrated into CellWalkerR package. 
CellWalker2 integrates different modalities of single-cell data, scRNA-Seq, scATAC-Seq or Multiomic. 
CellWalker2 annotates cells and compares cell type labels by scRNA-Seq, 
and identifities cell type-specific regulatory regions or other bulk-derived annotations by multiomic data. 
CellWalker2 can assign cells or bulk genomic annotations (e.g. TF motif or regulatory regions) to cell type hierarchy and compare different cell type hierarchies from different datasets. 
It also provides statistical signifiance of the association by permutation and visualization of cell labels and regulatory region mappings on a cell type hierarchy. 

## Installation

Install CellWalkR for R using devtools as follows:

``` r
$ R
> install.packages("devtools")
> devtools::install_github("PFPrzytycki/CellWalkR@dev")
```

## Usage

For a guide to using CellWalker2, see the provided
[vignette](examples/CellWalker2_Vignette.md), which covers the following:

1.  [Getting Started with CellWalker2](examples/CellWalker2_Vignette.md#getting-started-with-cellwalker2)
2.  [Use CellWalker2 for scRNA-Seq data](examples/CellWalker2_Vignette.md#use-cellwalker2-for-scRNA-Seq-data)
    1.  [Computing Cell-Cell
        Edges](examples/CellWalker2_Vignette.md#computing-cell-cell-edges)
    2.  [Computing Label-Cell
        Edges](examples/CellWalker2_Vignette.md#computing-label-cell-edges)
4.  [Tuning Label
    Edges](examples/CellWalker2_Vignette.md#tuning-label-edges)
5.  [Making a cellWalk
    Object](examples/CellWalker2_Vignette.md#making-a-cellwalk-object)
6.  [Adding Filters](examples/CellWalker2_Vignette.md#adding-filters)
7.  [Downstream
    Analysis](examples/CellWalker2_Vignette.md#downstream-analysis)
    1.  [Cell Labels](examples/CellWalker2_Vignette.md#cell-labels)
    2.  [Confusion
        Matrix](examples/CellWalker2_Vignette.md#confusion-matrix)
    3.  [Hierarchical Clustering of
        Labels](examples/CellWalker2_Vignette.md#hierarchical-clustering-of-labels)
    4.  [Plotting Cells](examples/CellWalker2_Vignette.md#plotting-cells)
    5.  [Bulk Data
        Mapping](examples/CellWalker2_Vignette.md#bulk-data-mapping)
8.  [Interactive
    Visualizaiton](examples/CellWalker2_Vignette.md#interactive-visualzation)
9.  [Adding a Second Set of
    Labels](examples/CellWalker2_Vignette.md#adding-a-second-set-of-labels)
10. [Detecting
    Doublets](examples/CellWalker2_Vignette.md#detecting-doublets)

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
