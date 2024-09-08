CellWalker2 introduction and installation
================
Zhirui Hu
2023-12-14


## Introduction

CellWalker2 extends the CellWalker method for combining scRNA-seq data with cell type ontologies to annotate cells and compare (hierarchically related) cell type labels, and combinding multiomic, scATAC-seq and scRNA-seq data with cell type labels and other bulk epigenetic data to interpret bulk level annotations at single cell level (see [paper](https://doi.org/10.1186/s13059-021-02279-1) for algorithmic details). This vignette shows examples of running CellWalker2 for different tasks: 
- Cell annotation and mapping cell types using scRNA-Seq data: a subset of two scRNA-seq datasets from human peripheral blood mononuclear cells (PBMCs)
- Map bulk-derived annotations to cell types: a subset of multiomic (scATAC-seq,  scRNA-seq and both) data from human developing cortex (see [paper](https://doi.org/10.1016/j.cell.2021.07.039)) to generate a cellWalk2 object. Then examples show how to assign bulk-derived brain region-specific predicted regulatory regions (pREs) to cell types and identify cell type-specific transcription factors in cortex.

## Prepare your data

### Assay data 

CellWalker2 uses a Gene-by-Cell count matrix as input for scRNA-seq data and a Cell-by-Peak count matrix as input for scATAC-seq data. scATAC-Seq data can be pre-processed by many pipelines that output a Cell-by-Peak matrix (e.g., CellRanger, ArchR, SnapATAC, Cicero). CellWalker2 provide optionally batch effect removel for multiple scRNA-Seq datasets using Seurat piplines but don't remove batch effect for scATAC-Seq data.  

### Cell type labels
Users can provide  marker genes for each cell type, and they can also provide a hierarchical tree relating cell tyeps to each other. If the marker genes are not provided, CellWalker2 will utilize Seurat pipeline to find differentially expressed genes per cell type.  

### Bulk derived-annotations
Additionally, genome coordinates can be provided to CellWalker2. Genome coordinates can be groups of regulatory regions (e.g. regulatory regions in different brain regions) or groups of regions with TF motifs or ChIPSeq peaks (for identify cell type-specific TFs).

## Getting Started with CellWalker2

Currently, CellWalker2 must be installed using devtools:


```r
install.packages("devtools")
devtools::install_github("PFPrzytycki/CellWalk@dev")
```

