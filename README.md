<img src="examples/CellWalkR_Vignette_files/figure-markdown_github/cellwalker_icon.png" id="id" class="class" width="50" height="50" /> CellWalkR
================

## About

CellWalkR is updated to Version 2! It combines CellWalker that integrates single-cell open chromatin
(scATAC-seq) data with cell type labels and bulk epigenetic data to
and CellWalker2 that integrate different modalities of single-cell data. 
CellWalkR can annotate cells and compare cell type labels by scATACSeq (CellWalker) and by scRNASeq (CellWalker2), 
and assign bulk identify cell type-specific regulatory regions by scATACSeq (CellWalker) or multiomic data (CellWalker2).
In addition, CellWalker2 can assign cells or bulk genomic annotations (e.g. TF motif or regulatory regions) to cell type hierarchy and compare different cell type hierarchies from different datasets. 
It also provides statistical signifiance of the association by permutation.  

## Installation

Install CellWalkR for R using devtools as follows:

``` r
$ R
> install.packages("devtools")
> devtools::install_github("PFPrzytycki/CellWalkR@dev")
```

## Usage

For a guide to using CellWalkR, see the provided
[readme](CellWalker.md) and [vignette](examples/CellWalkR_Vignette.md). 

For a guide to using CellWalkR2, see the provided
[readme](CellWalker2.md) and [vignette](examples/CellWalker2_Vignette.md). 

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
