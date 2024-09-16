Use CellWalker2 for multiomic data
================
Zhirui Hu
2023-12-14



In this vignette, we will use CellWalker2 to 1) assign cell type labels to predicted regulatory regions (pREs) identified from basal ganglia vs cortex (see [paper](https://doi.org/10.1016/j.cell.2020.06.002)) in human developing brain samples and 2) identify enriched TF motifs within pREs that tend to be more open in specific cell types (candidate cell type specific pREs).

#### Load data

We load subset of multiomic, scRNA-Seq and scATAC-Seq data from human developing cortex ([GSE162170](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162170)). For ATAC-Seq, we need to input a peak-by-cell count matrix accompanied with a peak coordinate matrix with the same rows. We can also load cell type tree as a phylo object and cell type markers. We rename some of the columns to make them readable by CellWalker2.


```r
# load single cell data
# 'data/SampleCellTypeTree.rda' should be loaded upon loading the package
devtools::load_all() # load SampleCortexSingleCellData
#> ℹ Loading CellWalkR

ATAC_Mat = SampleCortexSingleCellData$ATAC_Mat
peaks = SampleCortexSingleCellData$peaks
ATAC_Mat0 = SampleCortexSingleCellData$ATAC_Mat0
peaks0 = SampleCortexSingleCellData$peaks0
counts = SampleCortexSingleCellData$counts
counts2 = SampleCortexSingleCellData$counts2
RNA_markers = SampleCortexSingleCellData$RNA_markers
tr = SampleCellTypeTree

colnames(RNA_markers)[5] = 'p_val_adj' # will be used for select markers
colnames(peaks)[1:3] = c('seqnames', 'start', 'end')
```

#### Process RNASeq data

We first normalize and standarize scRNA-Seq data for both RNASeq part of multiomic data and unpaired RNASeq data. We don't need to compute cell type markers and build cell tree as they are already provided. We also don't need to compute cell-to-cell similarities as it will be computed by integrating all sources of data in the following step.


```r
dataset1 = processRNASeq(counts, do.findMarkers = F, computeKNN = F, computeSimilarity = F, 
                         buildTree = F) # RNASeq part of multiomic data, normalize data
dataset2 = processRNASeq(counts2, do.findMarkers = F, computeKNN = F, computeSimilarity = F, 
                         buildTree = F) # unpaired RNASeq, normalize data
```

#### Compute cell-to-cell type edges

Then we can compute cell to cell type edges based on standarized marker gene expression. Users can adjust the p-values and log2 fold change cutoff to select markers for computing edge weights. Here we just use the default values.


```r
# option 1: connect all cells with RNA data to cell types
labelEdges1 = computeTypeEdges(dataset1$expr_norm, RNA_markers)
labelEdges2 = computeTypeEdges(dataset2$expr_norm, RNA_markers)
labelEdges = rbind(labelEdges1, labelEdges2)
```

Users can also choose which cells to be connected to cell types.


```r
# option 2: connect only cells with unpaired RNA data to cell types
labelEdges = computeTypeEdges(dataset2$expr_norm, RNA_markers)
```

#### Construct cell-to-cell graph

Next, we construct cell-to-cell graph by integrating all the single cell data. In the current version of CellWalker2, multiomic data is required for this functionality, and unpaired scRNA-Seq or scATAC-Seq data can be added to the graph. CellWalker2 will compute cell-to-cell similarity based on gene expression and chromatin accessibility, construct a KNN graph, and output the shared K nearest neighbors between each pair of cells.

For gene expression, users can adjust how many principal components (PCs) to use to compute cell-to-cell distance (`ndim`). For chromatin accessibility, users can choose the similarity metric among Cosine, Jaccard, or LSI. Additionally, `ATAC_weight` specifies the weight of ATACSeq similarity for multiomic data (the weight for RNASeq is 1-`ATAC_weight`).


```r
cellgraph = constructCellGraph(counts, ATAC_Mat0, peaks0, counts2, ATAC_Mat, peaks) # with ATACSeq
```

We can also construct cell graph without unpaired scATAC-Seq:


```r
cellgraph = constructCellGraph(counts,ATAC_Mat0, peaks0, counts2) # without ATACSeq
```

#### Map region-specific pREs to cell types

We further input basal ganglia or cortex specific pREs and compute cell-to-label edges. In this case, the labels are *basal ganglia* and *cortex*. The edge weight is based on the chromatin accessibility of basal ganglia or cortex pREs in the cell.


```r
pRE = read.csv(system.file("extdata", "pRE_region_bg_cortex.csv", package = "CellWalkR"))
colnames(pRE)[c(1,11)] = c('seqnames', 'cluster') #rename some column names to be readable by CellWalker2

labelEdges1 = computeBulkEdges(pRE, peaks0, ATAC_Mat0)
labelEdges2 = computeBulkEdges(pRE, peaks, ATAC_Mat) 
labelEdges2 = rbind(labelEdges1, labelEdges2) # connect all ATAC cells to bulk annotations
```

Alternatively, we can only connect unpaired ATAC cells to bulk annotations:


```r
labelEdges2 = computeBulkEdges(pRE, peaks, ATAC_Mat) 
```

Finally, we assign cell type labels to basal ganglia or cortex specific pREs. By default, `compute.Zscore = T`, CellWalker2 will compute Z-scores by comparing the observed influence score with its null distribution generated by permuting cell-to-cell type edges. Otherwise, CellWalker2 will output influence score matrix from bulk labels to cell types only. `nround` is the number of randomization to compute Z-scores.

We set `labelEdgeWeight = NULL`, so CellWalker2 will tune the weight ratio between cell-to-label edges and cell-to-cell edges for both cell type labels and bulk annotations to minimize the entropy mapping each bulk annotation to cell types.




```r
cellWalk2 = annotateBulkRegion(cellgraph, labelEdges, labelEdges2, tr1 = tr, wtree = c(1, 0.1),
                               labelEdgeWeights = NULL, sampleDepth = 2000,  parallel = T, 
                               numCores = 8) # with tuning edgeWeights
#> tunning labelEdgeWeights...
#> labelHomogeneity at optimal edgeWeight:
#>    Var1 Var2 cellHomogeneity
#> 53 1000   10       0.1007592
#> run CellWalker:
```

If you don't want to tune `leblEdgeWeight`, you can use the following command using the default values or input other values.


```r
cellWalk2 = annotateBulkRegion(cellgraph, labelEdges, labelEdges2, tr1 = tr, wtree = c(1, 0.1))
```

We plot Z-score on the cell type tree to compare the cell type specificity of different bulk annotations.


```r
tr$node.label = colnames(cellWalk2$zscore)[-1:-(tr$Nnode+1)]
p1 = plotZscoreTree(tr, cellWalk2$zscore, cutoff = 3)
p1
```

![plot of chunk plotTree](CellWalker2_Vignette_files/figure-markdown_github/plotTree-1.png)

#### Identify cell type specific transcription factors (TFs)

First, we input the genomic coordinates of pREs as a bed file and identify TF motifs from JASPAR2020 within each pRE using Signac. For faster computation in this demo, we also filter for TF motifs that occur more frequently. The result file contains the genomic coordinates of each pRE and the motifs within it. Currently only hg38 is supported in `findMotifs` but users can input custom motifs file.


```r
pRE = read.table(system.file("extdata", "pRE-hg38.bed", package = "CellWalkR"))
colnames(pRE) = c('seqnames', 'start', 'end')
motifs = findMotifs(pRE)
#> Warning: Keys should be one or more alphanumeric characters followed by an underscore, setting key from peaks to peaks_
#> Building motif matrix
#> Finding motif positions
#> Creating Motif object

motifs = data.table::as.data.table(motifs)
motifs[, count:= .N, by = cluster] # select TFs appear in more pREs
motifs = as.data.frame(motifs[count > 1000])
```

Then, we compute cell-to-label edges. We treat each TF as a bulk label and each TF connects to each cell through the chromatin accessibilty of pREs containing the TF motif. We connect both cells from multiomic and unpaired ATACSeq to TFs.


```r
labelEdges1 = computeBulkEdges(motifs, peaks0, ATAC_Mat0)
labelEdges2 = computeBulkEdges(motifs, peaks, ATAC_Mat)
labelEdges2 = rbind(labelEdges1, labelEdges2)
```

Finally, we identify TFs whose motifs are enriched in pREs that are active in specific cell types. We convert `motifs` to a binary matrix showing the present/absence of motif in each pRE that can be used by `annotateBulkRegion`.

Different from the previous section, we randomly assign pREs to motifs to generate the null distribution for computing Z-scores. We set `groups1` to all zero as no permutation between cell-to cell types edges and `group2` to all one as permutations occur between all motifs and regions. We need to input count matrix and peaks for ATACSeq data as we need to recompute `labelEdges` for each randomization.


```r
regionMat = convertToMatrix(motifs) # a data.table with sequence name of pRE as the first column and clusters (TFs) as the following columns
groups1 =  rep(0, nrow(labelEdges)) # no permutation between cell-to cell types edges
groups2 = rep(1, nrow(regionMat)) # permutation between all motifs and regions
cellWalk2 = annotateBulkRegion(cellgraph, labelEdges, labelEdges2, groups1, groups2,
                               regionMat, list(ATAC_Mat0, ATAC_Mat), list(peaks0, peaks), 
                               tr1 = tr, wtree = c(1, 0.1), labelEdgeWeights = NULL, sampleDepth = 2000,
                               parallel = T, numCores = 8)
#> tunning labelEdgeWeights...
#> labelHomogeneity at optimal edgeWeight:
#>     Var1  Var2 cellHomogeneity
#> 18 10000 0.001      0.09224497
#> run CellWalker:
```

We can visualize Z-scores using dotplot. Each row is a cell type (including internal nodes on the tree), and each column is a TF. The internal node is names by two cell types among its desendents and the depth of the node.


```r
Zscore = cellWalk2$zscore[, tr$edge[,2]] # reorder cell types so that cell types closer in the tree will be closer on the heatmap
p1 = plotZscoreDotplot(Zscore, th = 5)
p1
```

![plot of chunk heatmap](CellWalker2_Vignette_files/figure-markdown_github/heatmap-1.png)

## Session Information


```r
sessionInfo()
#> R version 4.2.1 (2022-06-23)
#> Platform: aarch64-apple-darwin20 (64-bit)
#> Running under: macOS Monterey 12.6
#> 
#> Matrix products: default
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] CellWalkR_1.0.0             doParallel_1.0.17           iterators_1.0.14            foreach_1.5.2               splatter_1.21.1             scater_1.26.1               ggplot2_3.4.2              
#>  [8] scuttle_1.8.4               SingleCellExperiment_1.20.1 SummarizedExperiment_1.28.0 Biobase_2.58.0              MatrixGenerics_1.10.0       matrixStats_1.0.0           GenomicRanges_1.50.2       
#> [15] GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.2            BiocGenerics_0.44.0        
#> 
#> loaded via a namespace (and not attached):
#>   [1] rtracklayer_1.58.0                scattermore_1.2                   R.methodsS3_1.8.2                 SeuratObject_4.1.3                tidyr_1.3.0                      
#>   [6] JASPAR2020_0.99.10                bit64_4.0.5                       knitr_1.43                        irlba_2.3.5.1                     DelayedArray_0.24.0              
#>  [11] R.utils_2.12.2                    data.table_1.14.8                 KEGGREST_1.38.0                   TFBSTools_1.36.0                  RCurl_1.98-1.12                  
#>  [16] generics_0.1.3                    ScaledMatrix_1.6.0                callr_3.7.3                       cowplot_1.1.1                     usethis_2.2.1                    
#>  [21] RSQLite_2.3.1                     RANN_2.6.1                        future_1.32.0                     bit_4.0.5                         tzdb_0.4.0                       
#>  [26] spatstat.data_3.0-1               httpuv_1.6.11                     DirichletMultinomial_1.40.0       viridis_0.6.3                     xfun_0.39                        
#>  [31] hms_1.1.3                         evaluate_0.21                     promises_1.2.0.1                  fansi_1.0.4                       restfulr_0.0.15                  
#>  [36] caTools_1.18.2                    igraph_1.5.0                      DBI_1.1.3                         htmlwidgets_1.6.2                 spatstat.geom_3.2-1              
#>  [41] purrr_1.0.1                       ellipsis_0.3.2                    ggnewscale_0.4.9                  backports_1.4.1                   dplyr_1.1.2                      
#>  [46] annotate_1.76.0                   deldir_1.0-9                      sparseMatrixStats_1.10.0          vctrs_0.6.3                       remotes_2.4.2                    
#>  [51] ROCR_1.0-11                       abind_1.4-5                       cachem_1.0.8                      withr_2.5.0                       BSgenome.Hsapiens.UCSC.hg38_1.4.5
#>  [56] BSgenome_1.66.3                   progressr_0.13.0                  checkmate_2.2.0                   sctransform_0.3.5                 GenomicAlignments_1.34.1         
#>  [61] treeio_1.22.0                     prettyunits_1.1.1                 goftest_1.2-3                     cluster_2.1.4                     ape_5.7-1                        
#>  [66] lazyeval_0.2.2                    seqLogo_1.64.0                    crayon_1.5.2                      spatstat.explore_3.2-1            labeling_0.4.2                   
#>  [71] pkgconfig_2.0.3                   nlme_3.1-162                      vipor_0.4.5                       pkgload_1.3.2                     devtools_2.4.5                   
#>  [76] rlang_1.1.1                       globals_0.16.2                    lifecycle_1.0.3                   miniUI_0.1.1.1                    rsvd_1.0.5                       
#>  [81] rprojroot_2.0.3                   polyclip_1.10-4                   lmtest_0.9-40                     Matrix_1.5-4.1                    aplot_0.1.10                     
#>  [86] zoo_1.8-12                        beeswarm_0.4.0                    ggridges_0.5.4                    processx_3.8.1                    png_0.1-8                        
#>  [91] viridisLite_0.4.2                 rjson_0.2.21                      bitops_1.0-7                      R.oo_1.25.0                       KernSmooth_2.23-21               
#>  [96] Biostrings_2.66.0                 blob_1.2.4                        DelayedMatrixStats_1.20.0         stringr_1.5.0                     parallelly_1.36.0                
#> [101] spatstat.random_3.1-5             readr_2.1.4                       gridGraphics_0.5-1                CNEr_1.34.0                       beachmat_2.14.2                  
#> [106] scales_1.2.1                      memoise_2.0.1                     magrittr_2.0.3                    plyr_1.8.8                        ica_1.0-3                        
#> [111] zlibbioc_1.44.0                   compiler_4.2.1                    BiocIO_1.8.0                      RColorBrewer_1.1-3                fitdistrplus_1.1-11              
#> [116] Rsamtools_2.14.0                  cli_3.6.1                         XVector_0.38.0                    urlchecker_1.0.1                  listenv_0.9.0                    
#> [121] patchwork_1.1.2                   pbapply_1.7-0                     ps_1.7.5                          MASS_7.3-60                       tidyselect_1.2.0                 
#> [126] stringi_1.7.12                    RcppHungarian_0.2                 highr_0.10                        yaml_2.3.7                        locfit_1.5-9.8                   
#> [131] BiocSingular_1.14.0               ggrepel_0.9.3                     grid_4.2.1                        fastmatch_1.1-3                   tools_4.2.1                      
#> [136] future.apply_1.11.0               rstudioapi_0.14                   TFMPvalue_0.0.9                   gridExtra_2.3                     farver_2.1.1                     
#> [141] Rtsne_0.16                        digest_0.6.31                     shiny_1.7.4                       pracma_2.4.2                      motifmatchr_1.20.0               
#> [146] Rcpp_1.0.10                       later_1.3.1                       RcppAnnoy_0.0.20                  httr_1.4.6                        AnnotationDbi_1.60.2             
#> [151] colorspace_2.1-0                  brio_1.1.3                        XML_3.99-0.14                     fs_1.6.2                          tensor_1.5                       
#> [156] reticulate_1.34.0                 splines_4.2.1                     uwot_0.1.15                       yulab.utils_0.0.6                 RcppRoll_0.3.0                   
#> [161] tidytree_0.4.2                    spatstat.utils_3.0-3              sp_2.0-0                          ArchR_1.0.2                       ggplotify_0.1.0                  
#> [166] plotly_4.10.2                     sessioninfo_1.2.2                 xtable_1.8-4                      jsonlite_1.8.5                    ggtree_3.6.2                     
#> [171] poweRlaw_0.70.6                   testthat_3.1.9                    ggfun_0.1.1                       R6_2.5.1                          profvis_0.3.8                    
#> [176] pillar_1.9.0                      htmltools_0.5.5                   mime_0.12                         glue_1.6.2                        fastmap_1.1.1                    
#> [181] BiocParallel_1.32.6               BiocNeighbors_1.16.0              codetools_0.2-19                  pkgbuild_1.4.2                    Signac_1.10.0                    
#> [186] utf8_1.2.3                        lattice_0.21-8                    spatstat.sparse_3.0-2             tibble_3.2.1                      ggbeeswarm_0.7.2                 
#> [191] leiden_0.4.3                      gtools_3.9.4                      GO.db_3.16.0                      limma_3.54.2                      survival_3.5-5                   
#> [196] rmarkdown_2.22                    desc_1.4.2                        munsell_0.5.0                     GenomeInfoDbData_1.2.9            reshape2_1.4.4                   
#> [201] gtable_0.3.3                      Seurat_4.3.0.1
```
