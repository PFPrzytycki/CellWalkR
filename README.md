<img src="examples/CellWalkR_Vignette_files/figure-markdown_github/cellwalker2_icon.png" id="id" class="class" width="100" height="50" /> CellWalkR
================

## About

CellWalkR is updated to Version 2! This version combines the functionality of the **CellWalker** model with new features implemented in 
**CellWalker2**. 
As inputs, CellWalkerR takes cell type labels (lineages, states, etc; defined by their marker genes) and count data: 
a gene by cell matrix from scRNASeq data and/or a peak by cell matrix from scATACSeq data. 
Count matrices can be computed using a variety of upstream single-cell quantification software tools. 
Optionally, users may provide genome coordinates for sets of annotations they wish to map to 
cell types. Examples of annotations are bulk-derived regulatory elements, sequence motifs, genetic variants, 
or gene sets. CellWalkR builds a graph where the nodes are cells, cell types, and (if provided) annotations.
Then, CellWalkR uses a graph diffusion method to annotate cells, compare cell type labels, 
and assign cell type-specificity to annotations. 

The original **CellWalker** model integrates single-cell open chromatin (scATAC-seq) data with cell type labels and (optionally) 
bulk epigenetic data to annotate cells, compare cell type labels, and probabilitically assign cell type labels to annotations. 
**CellWalker2** extends this functionality with many new features, including:
+ unordered list
+ integrates different modalities of single-cell data (scATAC-seq, scRNA-seq, multiomic), 
+ enables integrative modeling of cells measured in different experiments,
+ incorporates hierarchical relationships between cell type labels,
+ compares cell type ontologies across contexts (conditions, species, datasets),
+ provide permutation-based measures of statistical significance for cell-to-label, region-to-label, and label-to-label mappings. 

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
