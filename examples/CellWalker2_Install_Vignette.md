CellWalker2 introduction and installation
================
Zhirui Hu
2023-12-14


## Introduction

CellWalker2 extends the CellWalker method for combining scRNA-seq data with cell type ontologies to annotate cells and compare (hierarchically related) cell type labels, and combinding multiomic, scATAC-seq and scRNA-seq data with cell type labels and other bulk epigenetic data to interpret bulk level annotations at single cell level (see [paper](https://doi.org/10.1186/s13059-021-02279-1) for algorithmic details). This vignette shows examples of running CellWalker2 for different tasks: 
- **Cell annotation and mapping cell types using scRNA-Seq data**: a subset of two scRNA-seq datasets from human peripheral blood mononuclear cells (PBMCs)
- **Map bulk-derived annotations to cell types**: a subset of multiomic (scATAC-seq,  scRNA-seq and both) data from human developing cortex (see [paper](https://doi.org/10.1016/j.cell.2021.07.039)) to generate a cellWalk2 object. Then examples show how to assign bulk-derived brain region-specific predicted regulatory regions (pREs) to cell types and identify cell type-specific transcription factors in cortex.

CellWalker2 will output Z-scores matrix for associations between labels or annotations. It also generates heatmaps or cell type tree plots of Z-scores.  

## Prepare your data

### Assay data 

CellWalker2 uses a Gene-by-Cell count matrix as input for scRNA-seq data and a Cell-by-Peak count matrix as input for scATAC-seq data. scATAC-Seq data can be pre-processed by many pipelines that output a Cell-by-Peak matrix (e.g., CellRanger, ArchR, SnapATAC, Cicero). CellWalker2 also uses multimodal data taking both RNA-Seq and ATAC-Seq matrices of the same cell barcodes as input. CellWalker2 provides optional batch effect removal for multiple scRNA-Seq datasets using Seurat piplines but doesn't remove batch effect for scATAC-Seq data.  

### Cell type labels
Users can provide  marker genes for each cell type, and can also provide a hierarchical tree relating cell types to each other. If the marker genes are not provided, CellWalker2 will utilize Seurat pipeline to find differentially expressed genes per cell type. If cell types are not provided, CellWalker2 will utilize Seurat pipeline to identify cell clusters using RNA-Seq data.  

### Bulk derived-annotations
Additionally, groups of genome coordinates can be provided to CellWalker2. Genome coordinates can be a group of regulatory regions (e.g. regulatory regions in a brain region), or regions with motifs or ChIPSeq peaks of a TF (for identifying cell type-specific TFs).

## Installation

Before installing CellWalker2, you might need to install some system libraries. If you use conda, you may create a new conda environment and install R, devtools and other system libraries. 
```bash
conda create -n cellWalker2 -c conda-forge r-base r-devtools
conda activate cellWalker2
# install system libraries
conda install -c conda-forge zlib openssl curl libxml2 bzip2 xz pcre2 gsl gmp glpk

```
Currently, CellWalker2 must be installed using devtools. Also CellWalker2 is only compatible with Seurat v4 and you might need to install it using remotes.


```r
$ R
if(!require(devtools))
{
  install.packages("devtools")
}
if (!require("BiocManager", quietly = TRUE))
{
  install.packages("BiocManager") # some dependencies are from Bioconductor
}

# install Seurat version 4 if you don’t have it
if(!require(remotes)){
  install.packages("remotes") # If not already installed
}
remotes::install_version("Seurat", version = "4.3.0.1")

devtools::install_github("PFPrzytycki/CellWalkR@cellwalker2") # # do not select updating Seurat package
```

