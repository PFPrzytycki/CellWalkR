CellWalker2 ntroduction and installation
================
Zhirui Hu
2023-12-14


## Introduction

CellWalker2 extends the CellWalker method for combining scRNA-seq data with cell type ontologies to annotate cells and compare (hierarchically related) cell type labels, and combinding multiomic, scATAC-seq and scRNA-seq data with cell type labels and other bulk epigenetic data to interpret bulk level annotations at single cell level (see [paper](https://doi.org/10.1186/s13059-021-02279-1) for algorithmic details). This vignette shows examples of running CellWalker2 on 1) a subset of two scRNA-seq datasets from human peripheral blood mononuclear cells (PBMCs)  and 2) a small subset of multiomic (scATAC-seq and scRNA-seq) data from human developing cortex (see [paper](https://doi.org/10.1016/j.cell.2021.07.039)) to generate a cellWalk2 object which can then be used to assign cell type labels to cells, map cell type trees as well as assign cell-type specific labels to annotations from bulk assays.

CellWalker2 uses a Gene-by-Cell count matrix as input for scRNA-seq data and a Cell-by-Peak matrix as input for scATAC-seq data. scATAC-Seq data can be pre-processed by many pipelines that output a Cell-by-Peak matrix (e.g., CellRanger, ArchR, SnapATAC, Cicero). Users define cell types via marker genes, and they can provide a hierarchical tree relating cell tyeps to each other. Additionally, genome coordinates for sets of annotations can be provided to CellWalker2.

## Getting Started with CellWalker2

Currently, CellWalker2 must be installed using devtools:


```r
install.packages("devtools")
devtools::install_github("PFPrzytycki/CellWalk@dev")
```

