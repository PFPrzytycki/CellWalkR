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
    1.  [load scRNA-Seq data](examples/CellWalker2_Vignette.md#load-scRNA-Seq-data)
    2.  [Cell type annotation](examples/CellWalker2_Vignette.md#cell-type-annotation)
    3.  [Cell types (trees) mapping](examples/CellWalker2_Vignette.md#Cell-types-trees-mapping)
3.  [Use CellWalker2 for multiomic data](examples/CellWalker2_Vignette.md#use-cellwalker2-for-multiomic-data)
    1.  [Load data](examples/CellWalker2_Vignette.md#load-data)
    2.  [Process RNASeq data](examples/CellWalker2_Vignette.md#process-rnaseq-data)
    3.  [Compute cell-to-cell type edges](examples/CellWalker2_Vignette.md#compute-cell-to-cell-type-edges)
    4.  [Construct cell-to-cell graph](examples/CellWalker2_Vignette.md#construct-cell-to-cell-graph)
    5.  [Map region-specific pREs to cell types](examples/CellWalker2_Vignette.md#map-region-specific-pREs-to-cell-types)
    6.  [Identify cell type specific TFs](examples/CellWalker2_Vignette.md#identify-cell-type-specific-tfs)
    

If you use CellWalkR please cite:

1.  Przytycki, P.F., Pollard, K.S. “CellWalkR: An R Package for
    integrating and visualizing single-cell and bulk data to resolve
    regulatory elements.” *Bioinformatics* (2022).
    <https://doi.org/10.1093/bioinformatics/btac150>

2.  Przytycki, P.F., Pollard, K.S. “CellWalker integrates single-cell
    and bulk data to resolve regulatory elements across cell types in
    complex tissues.” *Genome Biology* (2021).
    <https://doi.org/10.1186/s13059-021-02279-1>

