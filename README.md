<img src="examples/CellWalkR_Vignette_files/figure-markdown_github/cellwalker2_icon.png" id="id" class="class" width="100" height="50" /> CellWalkR
================

## About

CellWalkR is updated to Version 2! This version combines the functionality of the **CellWalker** model with new features implemented in 
**CellWalker2**. 
As inputes, CellWalkerR takes cell type labels and count data, either gene by cell matrix for scRNASeq data or peak by cell matrix for scATACSeq data. 
Count matrices can be computed using a variety of upstream single-cell quantification software tools. 
CellWalkR can annotate cells, compare cell type labels, and assign cell type-specificity to genomic regions. 

The original **CellWalker** model integrates single-cell open chromatin (scATAC-seq) data with cell type labels and (optional) 
bulk epigenetic data to annotate cells, compare cell type labels, and probabilitically assign cell type labels to genomic coordinates (e.g., 
bulk-derived regulatory elements, sequence motifs, genetic variants). 

**CellWalker2** integrates different modalities of single-cell data (scATAC-seq, scRNA-seq, multiomic), 
including cells measured in different experiments. It also incorporates hierarchical relationships between cell type labels and enables comparisons of
cell type ontologies across contexts (conditions, species, datasets). 
Finally, **CellWalker2** adds permutation methods that provide statistical significance for cell-to-label, region-to-label, and label-to-label mappings. 

## Installation

Install CellWalkR for R using devtools as follows:

``` r
$ R
> install.packages("devtools")
> devtools::install_github("PFPrzytycki/CellWalkR@dev")
```

## Usage

For a guide to using CellWalker, see the provided
[readme](CellWalker.md) and [vignette](examples/CellWalkR_Vignette.md). 

For a guide to using CellWalker2, see the provided
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
