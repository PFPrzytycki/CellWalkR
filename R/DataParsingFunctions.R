#' Compute Cell Similarity
#'
#' \code{computeCellSim} computes a cell-to-cell similarity matrix using the given method.
#' Currently only Jaccard similarity is implemented.
#'
#' @param ATACMat either a cell-by-peak matrix or a SnapATAC object
#' @param method method used to compute cell similarity
#' @return matrix of cell-to-cell similarity
#' @export
computeCellSim = function(ATACMat, method="Jaccard"){
  if(missing(ATACMat)){
    stop("Must provide either a cell-by-peak matrix or a SnapATAC object")
  }
  if(method!="Jaccard"){
    stop("Currently only Jaccard similarity is implemented")
  }

  if(is(ATACMat,"snap")){
    if(!requireNamespace("SnapATAC", quietly = TRUE)){
      stop("Must install SnapATAC")
    }
    ATACMat = SnapATAC::makeBinary(ATACMat, mat="bmat")
    ATACMat = SnapATAC::runJaccard(obj=ATACMat, mat="bmat", tmp.folder=tempdir(), max.var = dim(ATACMat)[1])
    ATACMat@jmat@jmat
  }
  else if(is(ATACMat,"Matrix")|is(ATACMat,"matrix")){
    if(min(dim(ATACMat))==0){
      stop("Not enough cells or peaks in matrix")
    }
    if(dim(ATACMat)[1]>dim(ATACMat)[2]){
      warning("ATACMat has more cells than peaks, does it need to be transposed?")
    }
    sparseJaccard(ATACMat>0)
  }
  else{
    stop("Must provide either a cell-by-peak matrix or a SnapATAC object")
  }
}

#' Get genomic regions
#'
#' \code{getRegions} returns a GRanges object of how gene name identifiers map to
#' genomic regions either based on gene promoters or the full gene body
#'
#' @param geneBody boolean, use gene body or just promoter
#' @param genome which genome build (hg38, mm10)
#' @param names which gene names to use (Entrez, Ensembl)
#' @return GRanges list of regions with gene_id
#' @export
getRegions = function(geneBody=TRUE, genome="hg38", names="Entrez"){
  if(!requireNamespace("GenomicFeatures", quietly = TRUE)){
    stop("Please install GenomicFeatures to use getRegions")
  }
  if(genome=="hg38"){
    if(names=="Ensembl"){
      stop("Not an implemented gene name identifier")
    }
    else if(names=="Entrez"){
      if(!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)){
        stop("Must install TxDb.Hsapiens.UCSC.hg38.knownGene")
      }
      TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    }
    else{
      stop("Not a valid gene name identifier")
    }
  }
  else if(genome=="mm10"){
    if(names=="Ensembl"){
      if(!requireNamespace("TxDb.Mmusculus.UCSC.mm10.ensGene", quietly = TRUE)){
        stop("Must install TxDb.Mmusculus.UCSC.mm10.ensGene")
      }
      TxDb = TxDb.Mmusculus.UCSC.mm10.ensGene::TxDb.Mmusculus.UCSC.mm10.ensGene
    }
    else if(names=="Entrez"){
      if(!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE)){
        stop("Must install TxDb.Mmusculus.UCSC.mm10.knownGene")
      }
      TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    }
    else{
      stop("Not a valid gene name identifier")
    }
  }
  else{
    stop("Not an implemented genome build")
  }

  regions = GenomicRanges::promoters(GenomicFeatures::genes(TxDb))
  if(geneBody){
    regions = c(regions, GenomicFeatures::genes(TxDb))
  }

  regions
}

#' Map peaks to genes
#'
#' \code{mapPeaksToGenes} Generates a mapping from genes to peaks per label
#'
#' @param labelGenes data.frame with genes of interest in first column and corresponding
#' labels in second column
#' @param ATACMat cell-by-peak matrix
#' @param peaks GRanges object of peaks (equal to number of rows in ATACMat)
#' @param regions GRanges object of genomic regions associated with genes
#' @return list of genes and corresponding peaks by label
#' @export
mapPeaksToGenes = function(labelGenes, ATACMat, peaks, regions){
  if(missing(labelGenes)){
    stop("Must provide a table with genes of interest in first column and corresponding labels in second column")
  }
  if(is.null(dim(labelGenes)) || dim(labelGenes)[2]<2){
    stop("Must provide a table with genes of interest in first column and corresponding labels in second column")
  }
  labels = unique(labelGenes[,2])
  if(length(labels)==1){
    warning("Only one label exists")
  }
  if(missing(ATACMat)){
    stop("Must provide either a cell-by-peak matrix or a SnapATAC object")
  }
  if(is(ATACMat,"snap")){
    if(!requireNamespace("SnapATAC", quietly = TRUE)){
      stop("Must install SnapATAC")
    }
    if(dim(ATACMat@gmat)[1]==0){
      stop("No gmat, call SnapATAC function createGmatFromMat first")
    }
    whichGenes = match(unique(labelGenes[,1]),colnames(ATACMat@gmat))
    if(length(whichGenes)<10){
       warning("Fewer than 10 gene names match")
    }
    stop("No need to map to genes, use gmat")
  }
  else{
    if(missing(regions)){
      stop("Must provide a GRanges object of genomic regions associatated with genes")
    }
    if(!is(regions,"GRanges")){
      stop("regions must be a GRanges object")
    }
    if(length(regions)==0){
      stop("regions must contain at least one genomic range")
    }
    if(length(regions$gene_id)!=length(regions)){
      stop("Each region must have an associated gene_id")
    }
    whichGenes = match(unique(labelGenes[,1]),regions$gene_id)
    if(length(whichGenes)<10){
      warning("Fewer than 10 gene names match")
    }
    if(missing(peaks)){
      stop("Must provide a GRanges object of peaks")
    }
    if(!is(peaks,"GRanges")){
      stop("peaks must be a GRanges object")
    }
    if(!is(ATACMat,"Matrix") & !is(ATACMat,"matrix")){
      stop("Must provide either a cell-by-peak matrix or a SnapATAC object")
    }
    if(min(dim(ATACMat))==0){
      stop("Not enough cells or peaks in matrix")
    }
    if(dim(ATACMat)[1]>dim(ATACMat)[2]){
      warning("ATACMat has more cells than peaks, does it need to be transposed?")
    }
    if(!length(peaks)==dim(ATACMat)[2]){
      stop("Number of peaks must be the same as the number of columns in ATACMat")
    }
    cellTypePeakMap = sapply(labels, function(x) {
      markerGenes = labelGenes[labelGenes[,2]==x,1]
      markerGenes = markerGenes[markerGenes %in% regions$gene_id]
      markerLocs = regions[match(markerGenes, regions$gene_id)]
      markerOverlaps = GenomicRanges::findOverlaps(peaks,markerLocs)
      list(peak=markerOverlaps@from,gene=markerGenes[(markerOverlaps@to-1)%%length(markerGenes)+1])
    })
  }

}

#' Compute Label Edges
#'
#' \code{computeLabelEdges} generates a matrix of edges from each label to each cell.
#' Edges are scaled by the passed weights and filters. Each filter can either apply
#' just to specific peaks or to whole genes (filterGene, TRUE by default) and can be permissive or
#' filtering out (filterOut, FALSE by default). Regions mapping peaks to genes must be
#' provided to use filters. If filters are applied to specific peaks, peaks for the ATAC
#' data must also be provided.
#'
#' @param labelGenes data.frame with genes of interest in first column and corresponding
#' labels in second column, optionally log-fold change in expression values in the third
#' @param ATACMat either a cell-by-peak matrix or a SnapATAC object
#' @param ATACGenePeak per label mapping of peaks to genes returned by mapPeaksToGenes()
#' @param method scaling method, either "Expression","Correlation", or "None"
#' @param filters list of GRanges lists for locations that pass filters
#' @param filterWeights numeric vector of weights assigned to each filter
#' @param filterOut boolean for each filter of whether those regions are permitted or filtered out
#' @param filterGene boolean for each filter of whether it applies to genes as a whole or just overlapping peaks (note that weights do not apply to peaks)
#' @param regions regions map peaks to genes
#' @param peaks GRanges of peaks in ATACMat
#' @return matrix of weights from each label to each cell
#' @export
computeLabelEdges = function(
  labelGenes,
  ATACMat,
  ATACGenePeak,
  method="Expression",
  filters,
  filterWeights,
  filterOut,
  filterGene,
  regions,
  peaks){
  if(missing(labelGenes)){
    stop("Must provide a table with genes of interest in first column and corresponding labels in second column")
  }
  if(is.null(dim(labelGenes)) || dim(labelGenes)[2]<2){
    stop("Must provide a table with genes of interest in first column and corresponding labels in second column")
  }
  labels = unique(labelGenes[,2])
  if(length(labels)==1){
    warning("Only one label exists")
  }

  if(missing(ATACMat)){
    stop("Must provide either a cell-by-peak matrix or a SnapATAC object")
  }
  if(is(ATACMat,"snap")){
    stop("Not yet implemented")
  }
  else{
    if(!is(ATACMat,"Matrix") & !is(ATACMat,"matrix")){
      stop("Must provide either a cell-by-peak matrix or a SnapATAC object")
    }
    if(min(dim(ATACMat))==0){
      stop("Not enough cells or peaks in matrix")
    }
    if(dim(ATACMat)[1]>dim(ATACMat)[2]){
      warning("ATACMat has more cells than peaks, does it need to be transposed?")
    }
    if(missing(ATACGenePeak)){
      stop("Must provide a per label mapping of peaks to genes returned by mapPeaksToGenes()")
    }

    if(!missing(filters)){
      if(missing(regions)){
        stop("Must provide a GRanges object of genomic regions associatated with genes to apply filters")
      }
      if(!is(regions,"GRanges")){
        stop("regions must be a GRanges object")
      }
      if(length(regions)==0){
        stop("regions must contain at least one genomic range")
      }
      if(length(regions$gene_id)!=length(regions)){
        stop("Each region must have an associated gene_id")
      }
      if(missing(filterWeights) || length(filters)!=length(filterWeights)){
        warning("Assuming all filters have full weight")
        filterWeights = rep(1, length(filters))
      }
      if(missing(filterOut) || length(filters)!=length(filterOut)){
        warning("Assuming all filters are permissive")
        filterOut = rep(FALSE, length(filters))
      }
      if(missing(filterGene) || length(filters)!=length(filterGene)){
        warning("Assuming all filters apply at the gene level")
        filterGene = rep(TRUE, length(filters))
      }
      geneFilter = rep(1, length(labelGenes[,1]))
      for(filterIndex in which(filterGene)){
        filterRanges = filters[[filterIndex]]
        regionFilterOverlap = countOverlaps(regions, filterRanges)
        if(filterOut[filterIndex]){
          geneFilterHit = labelGenes[,1] %in% names(which(regionFilterOverlap>0))
          geneFilter[geneFilterHit] = geneFilter[geneFilterHit]*(1-filterWeights[filterIndex])
        }
        else{
          geneFilterHit = labelGenes[,1] %in% names(which(regionFilterOverlap>0))
          geneFilter[!geneFilterHit] = geneFilter[!geneFilterHit]*(1-filterWeights[filterIndex])
        }
      }
    }
    else{
      geneFilter = rep(1, length(labelGenes[,1]))
    }

    cellPeakCounts = Matrix::rowSums(ATACMat)
    cellsInMarkers = c()
    for(label in labels){
      if(!label %in% dimnames(ATACGenePeak)[[2]]){
        stop("A label is missing from the peak to gene mapping. Is it the right mapping?")
      }
      peaksInMarkers = ATACGenePeak[["peak",label]]
      genesInMarkers = ATACGenePeak[["gene",label]]
      if(max(peaksInMarkers)>dim(ATACMat)[2]){
        stop("Indecies dont match in peak to gene mapping. Is it the right mapping?")
      }
      if(length(which(genesInMarkers %in% labelGenes[,1]))!=length(genesInMarkers)){
        stop("Genes dont match in peak to gene mapping. Is it the right mapping?")
      }
      if(method=="None"){
        geneExp = rep(1, length(genesInMarkers))
      }
      else if(method=="Correlation"){
        stop("Correlation not yet implemented")
      }
      else if(method=="Expression"){
        if(dim(labelGenes)[2]<3){
          stop("Must provide a table with genes of interest in first column and corresponding labels in
               second column and log fold gene expression in the third")
        }
        geneExp = labelGenes[match(genesInMarkers,labelGenes[labelGenes[,2]==label,1]),3]
      }
      else{
        stop("Not a recognized scaling method")
      }

      geneFilterScale = geneFilter[match(genesInMarkers,labelGenes[labelGenes[,2]==label,1])]

      if(!missing(filters) && sum(!filterGene)>0){
        if(missing(peaks) || !is(peaks,"GRanges")){
          stop("Must provide a GRanges object of peaks to apply filters to peaks")
        }
        if(!length(peaks)==dim(ATACMat)[2]){
          stop("Number of peaks must be the same as the number of columns in ATACMat")
        }

        for(filterIndex in which(!filterGene)){
          filterRanges = filters[[filterIndex]]
          peakFilterOverlap = countOverlaps(peaks, filterRanges)
          if(filterOut[filterIndex]){
            peaksInMarkers = peaksInMarkers[peakFilterOverlap==0]
            geneExp = geneExp[peakFilterOverlap==0]
            geneFilterScale = geneFilterScale[peakFilterOverlap==0]
          }
          else{
            peaksInMarkers = peaksInMarkers[peakFilterOverlap>0]
            geneExp = geneExp[peakFilterOverlap>0]
            geneFilterScale = geneFilterScale[peakFilterOverlap>0]
          }
        }
      }
      cellsInLabelMarkers = Matrix::colSums(Matrix::t(ATACMat[,peaksInMarkers])*geneExp*geneFilterScale)/cellPeakCounts
      cellsInLabelMarkers[cellsInLabelMarkers<0] = 0
      cellsInMarkers = cbind(cellsInMarkers,cellsInLabelMarkers)
    }
    colnames(cellsInMarkers) = labels
    cellsInMarkers
  }
}

#' Run CellWalker Method
#'
#' \code{walkCells} applies the CellWalker method to the given cell-to-cell and
#' cell-to-label edges.
#'
#' @param cellEdges matrix of cell-to-cell similarity (c-by-c)
#' @param labelEdgesList list of cell-to-label similarity matrices each with dimensions c-by-l
#' @param labelEdgeWeights vector of weights to test for label edges
#' @param sampleDepth depth to subsample cells to for faster calculation
#' @return cellWalk object with influence matrix and labels
#' @export
walkCells = function(cellEdges, labelEdgesList, labelEdgeWeights, sampleDepth){
  if(missing(cellEdges)){
    stop("Must provide a matrix of cell-to-cell similarity")
  }
  if(!is(cellEdges, "matrix") & !is(cellEdges,"Matrix")){
    stop("Cell-to-cell similarity must be a matrix")
  }
  if(dim(cellEdges)[1] != dim(cellEdges)[2]){
    stop("Cell-to-cell similarity matrix must have dimensions c-by-c")
  }

  c = dim(cellEdges)[1]
  if(missing(sampleDepth)){
    sampleDepth = c
  }
  else if(sampleDepth>c){
    warning("Sample depth is greater than number of cells, will use all cells")
    sampleDepth = c
  }

  if(missing(labelEdgesList) || !is(labelEdgesList, "list")){
    stop("Must provide a list of cell-to-label similarity matrices")
  }
  if(sum(sapply(labelEdgesList, function(x) dim(x)[1])!=c)!=0){
    stop("Each cell-to-label similarity matrix must have dimensions c-by-l, number of rows equal to the number of cells")
  }

  if(missing(labelEdgeWeights) || !length(labelEdgeWeights)==length(labelEdgesList)){
    stop("Must provide a weight for each set of cell-to-label edges")
  }

  l = sum(sapply(labelEdgesList, function(x) dim(x)[2]))
  combinedGraph = combineMultiLabelGraph(labelEdgesList, cellEdges, labelEdgeWeights)

  if(! is.null(colnames(cellEdges))){
    colnames(combinedGraph) = c( c(unlist(sapply(labelEdgesList, function(x) colnames(x)))), colnames(cellEdges))
    rownames(combinedGraph) = c( c(unlist(sapply(labelEdgesList, function(x) colnames(x)))), colnames(cellEdges))
  }

  if(!sampleDepth == c){
    selectCells = sample(dim(cellEdges)[1], sampleDepth)+l
    combinedGraph = combinedGraph[c(1:l,selectCells),c(1:l,selectCells)]
  }

  combinedGraph = combinedGraph[Matrix::colSums(combinedGraph)!=0, Matrix::colSums(combinedGraph)!=0]
  infMat = randomWalk(combinedGraph)

  if(! is.null(colnames(cellEdges))){
    rownames(infMat) = rownames(combinedGraph)
    colnames(infMat) = colnames(combinedGraph)
  }

  normMat = normalizeInfluence(infMat[-(1:l),(1:l)])
  colnames(normMat) = c(unlist(sapply(labelEdgesList, function(x) colnames(x))))
  if(max(table(colnames(normMat)))>1){
    warning("Multiple labels share the same name")
  }
  cellLabels = apply(normMat, 1, function(x) colnames(normMat)[order(x, decreasing = TRUE)][1])

  cellWalk = list(infMat=infMat, normMat=normMat, cellLabels=cellLabels)
  class(cellWalk) = "cellWalk"
  cellWalk
}

#' Tune edge weights
#'
#' \code{tuneEdgeWeights} generates cell homogeneity scores for possible edge
#' weight settings across all labelEdge matrices. To speed up calculations cells
#' can be subsampled to a specified sampleDepth (default 10000).
#'
#' @param cellEdges matrix of cell-to-cell similarity (c-by-c)
#' @param labelEdgesList list of cell-to-label similarity matrices each with dimensions c-by-l
#' @param labelEdgeOpts vector of weights to test for label edges
#' @param sampleDepth integer depth to subsample cells to for faster calculation
#' @param cellTypes list of names of cell types, if not provided unique cell labels are used
#' @param parallel boolean, execute in parallel
#' @param numCores integer, number of cores to use for parallel execution
#' @param trackProgress boolean, print percent run completed
#' @return data frame of weights with corresponding cell homogeneity
#' @export
tuneEdgeWeights = function(
  cellEdges,
  labelEdgesList,
  labelEdgeOpts = 10^seq(-4,4,1),
  sampleDepth = 10000,
  cellTypes,
  parallel = FALSE,
  numCores = 1,
  trackProgress = FALSE ){

  if(missing(cellEdges)){
    stop("Must provide a matrix of cell-to-cell similarity")
  }
  if(!is(cellEdges, "matrix") & !is(cellEdges,"Matrix")){
    stop("Cell-to-cell similarity must be a matrix")
  }
  if(dim(cellEdges)[1] != dim(cellEdges)[2]){
    stop("Cell-to-cell similarity matrix must have dimensions c-by-c")
  }
  c = dim(cellEdges)[1]
  if(sampleDepth>c){
    warning("Sample depth is greater than number of cells, will use all cells")
    sampleDepth = c
  }

  if(missing(labelEdgesList) || !is(labelEdgesList, "list")){
    stop("Must provide a list of cell-to-label similarity matrices")
  }
  if(sum(sapply(labelEdgesList, function(x) dim(x)[1])!=c)!=0){
    stop("Each cell-to-label similarity matrix must have dimensions c-by-l, number of rows equal to the number of cells")
  }

  numLabelLists = length(labelEdgesList)
  l = sum(sapply(labelEdgesList, function(x) dim(x)[2]))

  parameterCombos = expand.grid(rep(list(labelEdgeOpts),numLabelLists))

  cellTypesMissing = FALSE
  if(missing(cellTypes)){
    cellTypesMissing = TRUE
  }

  if(parallel){
    if(!requireNamespace("parallel", quietly = TRUE)){
      stop("Must install parallel")
    }
    if(numCores >= parallel::detectCores()) {
      numCores = parallel::detectCores() - 1
      warning("numCores greater than avaialble cores, using available cores")
    }
    if(trackProgress){
      warning("Cant track progress in parallel")
    }

    cellHomogeneity = parallel::mcmapply(function(param){
      theseWeights = as.numeric(parameterCombos[param,])

      cellWalk = walkCells(cellEdges, labelEdgesList, theseWeights, sampleDepth)

      if(cellTypesMissing){
        cellTypes = unique(cellWalk[["cellLabels"]])
      }

      computeCellHomogeneity(cellWalk, cellTypes)

    }, param=1:dim(parameterCombos)[1], mc.cores = numCores)
  }
  else{
    cellHomogeneity = c()
    for(param in 1:dim(parameterCombos)[1]){
      theseWeights = as.numeric(parameterCombos[param,])

      cellWalk = walkCells(cellEdges, labelEdgesList, theseWeights, sampleDepth)

      if(cellTypesMissing){
        cellTypes = unique(cellWalk[["cellLabels"]])
      }

      cellHomogeneity = c(cellHomogeneity, computeCellHomogeneity(cellWalk, cellTypes))

      if(trackProgress & sum(sapply(seq(.1,.9,.1), function(x)
        (param)/dim(parameterCombos)[1] >= x &
        (param-1)/dim(parameterCombos)[1] <= x))>=1){
        print(paste0(as.integer(param/dim(parameterCombos)[1]*100), "% done"))
      }

    }
  }

  data.frame(parameterCombos, cellHomogeneity)

}

#' Tune filter weights
#'
#' \code{tuneFilterWeights} generates cell homogeneity scores for filter weight settings (by default
#' 0 and 1 for on and off)
#' across all filters. Each filter can either apply just to specific peaks or to whole genes
#' (filterGene, TRUE by default) and can be permissive or
#' filtering out (filterOut, FALSE by default). Regions mapping peaks to genes must be
#' provided to use filters. If filters are applied to specific peaks, peaks for the ATAC
#' data must also be provided. To speed up calculations cells
#' can be subsampled to a specified sampleDepth (default 5000).
#'
#' @param cellEdges matrix of cell-to-cell similarity (c-by-c)
#' @param labelGenesList list of tables with genes of interest in first column and corresponding
#' labels in second column, optionally log-fold change in expression values in the third
#' @param labelEdgesList list of cell-to-label similarity matrices each with dimensions c-by-l
#' @param labelEdgeWeights vector of edge weights corresponding to list of cell-to-label similarity matrices
#' @param ATACMat either a cell-by-peak matrix or a SnapATAC object
#' @param ATACGenePeak per label mapping of peaks to genes returned by mapPeaksToGenes()
#' @param method scaling method, either "Expression","Correlation", or "None"
#' @param filters list of GRanges for locations that pass filters
#' @param filterWeightOpts vector of parameters to test for weights assigned to each filter
#' @param filterOut boolean for each filter of whether those regions are permitted or filtered out
#' @param filterGene boolean for each filter of whether it applies to genes as a whole or just overlapping peaks (note that weights do not apply to peaks)
#' @param regions regions map peaks to genes
#' @param peaks GRanges of peaks in ATACMat
#' @param sampleDepth integer, depth to subsample cells to for faster calculation
#' @param cellTypes list of names of cell types, if not provided unique cell labels are used
#' @param parallel boolean, execute in parallel
#' @param numCores integer, number of cores to use for parallel execution
#' @param trackProgress boolean, print percent run completed
#' @return data frame of parameters with corresponding cell homogeneity
#' @export
tuneFilterWeights = function(
  cellEdges,
  labelGenesList,
  labelEdgesList,
  labelEdgeWeights,
  ATACMat,
  ATACGenePeak,
  method = "Expression",
  filters,
  filterWeightOpts=c(0,1),
  filterOut,
  filterGene,
  regions,
  peaks,
  sampleDepth=5000,
  cellTypes,
  parallel = FALSE,
  numCores = 1,
  trackProgress = FALSE){

  if(missing(cellEdges)){
    stop("Must provide a matrix of cell-to-cell similarity")
  }
  if(!is(cellEdges, "matrix") & !is(cellEdges,"Matrix")){
    stop("Cell-to-cell similarity must be a matrix")
  }
  if(dim(cellEdges)[1] != dim(cellEdges)[2]){
    stop("Cell-to-cell similarity matrix must have dimensions c-by-c")
  }
  c = dim(cellEdges)[1]
  if(sampleDepth>c){
    warning("Sample depth is greater than number of cells, will use all cells")
    sampleDepth = c
  }

  if(missing(labelGenesList) || !is(labelGenesList, "list")){
    stop("Must provide a  list of tables with genes of interest in first column and corresponding
         labels in second column, optionally log-fold change in expression values in the third")
  }

  if(missing(labelEdgesList) || !is(labelEdgesList, "list")){
    stop("Must provide a list of cell-to-label similarity matrices")
  }
  if(sum(sapply(labelEdgesList, function(x) dim(x)[1])!=c)!=0){
    stop("Each cell-to-label similarity matrix must have dimensions c-by-l, number of rows equal to the number of cells")
  }

  if(missing(labelEdgeWeights) || !length(labelEdgeWeights)==length(labelEdgesList)){
    stop("Must provide a weight for each set of cell-to-label edges")
  }

  if(missing(ATACMat)){
    stop("Must provide either a cell-by-peak matrix or a SnapATAC object")
  }
  if(is(ATACMat,"snap")){
    stop("Not yet implemented")
  }
  if(missing(ATACGenePeak)){
    stop("Must provide a per label mapping of peaks to genes returned by mapPeaksToGenes()")
  }

  if(missing(filters) || !is(filters, "list")){
    stop("Must provide list of GRanges for locations that pass filters")
  }
  if(missing(regions)){
      stop("Must provide a GRanges object of genomic regions associatated with genes to apply filters")
  }

  numFilters = length(filters)
  parameterCombos = expand.grid(rep(list(filterWeightOpts),numFilters))

  cellTypesMissing = FALSE
  if(missing(cellTypes)){
    cellTypesMissing = TRUE
  }

  if(parallel){
    if(!requireNamespace("parallel", quietly = TRUE)){
      stop("Must install parallel")
    }
    if(numCores >= parallel::detectCores()) {
      numCores = parallel::detectCores() - 1
      warning("numCores greater than avaialble cores, using available cores")
    }
    if(trackProgress){
      warning("Cant track progress in parallel")
    }
    if(missing(filterOut)){
      filterOut = rep(FALSE, length(filters))
    }
    if(missing(filterGene)){
      filterGene = rep(TRUE, length(filters))
    }
    if(missing(peaks)){
      if(sum(!filterGene)>0 && missing(peaks)){
        stop("Must provide a GRanges object of peaks to apply filters to peaks")
      }
      else{
        peaks = NULL
      }
    }
    cellHomogeneity = parallel::mcmapply(function(param){
      theseWeights = as.numeric(parameterCombos[param,])
      labelEdgesList = list()
      for(labelGenes in labelGenesList){
        newEdges = computeLabelEdges(labelGenes, ATACMat, ATACGenePeak, method,
                                     filters, theseWeights,filterOut,
                                     filterGene, regions,peaks)
        labelEdgesList = c(labelEdgesList,list(newEdges))
      }
      l = sum(sapply(labelEdgesList, function(x) dim(x)[2]))

      cellWalk = walkCells(cellEdges, labelEdgesList, labelEdgeWeights, sampleDepth)

      if(cellTypesMissing){
        cellTypes = unique(cellWalk[["cellLabels"]])
      }

      computeCellHomogeneity(cellWalk, cellTypes)

    }, param=1:dim(parameterCombos)[1], mc.cores = numCores)
  }
  else{
    cellHomogeneity = c()
    for(param in 1:dim(parameterCombos)[1]){
      #add option for filter to only apply to some cells? to some labels?
      theseWeights = as.numeric(parameterCombos[param,])
      labelEdgesList = list()
      for(labelGenes in labelGenesList){
        newEdges = computeLabelEdges(labelGenes, ATACMat, ATACGenePeak, method,
                                     filters, theseWeights,filterOut,
                                     filterGene, regions,peaks)
        labelEdgesList = c(labelEdgesList,list(newEdges))
      }
      l = sum(sapply(labelEdgesList, function(x) dim(x)[2]))

      cellWalk = walkCells(cellEdges, labelEdgesList, labelEdgeWeights, sampleDepth)
      infMat = cellWalk[["infMat"]]
      normMat = cellWalk[["normMat"]]
      cellLabels = cellWalk[["cellLabels"]]

      if(cellTypesMissing){
        cellTypes = unique(cellWalk[["cellLabels"]])
      }

      cellHomogeneity = c(cellHomogeneity, computeCellHomogeneity(cellWalk, cellTypes))

      if(trackProgress & sum(sapply(seq(.1,.9,.1), function(x)
        (param)/dim(parameterCombos)[1] >= x &
        (param-1)/dim(parameterCombos)[1] <= x))>=1){
        print(paste0(as.integer(param/dim(parameterCombos)[1]*100), "% done"))
      }

    }
  }

  data.frame(parameterCombos, cellHomogeneity)
}

#' Find markers in scRNA-seq data
#'
#' \code{findMarkers} Uses Seurat to find markers in scRNA-seq data and returns
#' a table of markers and labels that can be used by CellWalker. This implementation
#' uses default parameters and is designed for a simple first pass for creating
#' labels.
#'
#' @param RNAMat gene-by-barcode matrix of scRNA data
#' @param genes character vector of gene names
#' @param barcodes character vector of barcodes
#' @param dims integer vector of PCA dimensions to use
#' @param resolution numeric resolution to use in Louvain clustering
#' @return table of marker genes, labels, and logFC in expression
#' @export
findMarkers = function(RNAMat, genes, barcodes, dims=1:10, resolution=0.5){
  if(!requireNamespace("Seurat", quietly = TRUE)){
    stop("Must install Seurat")
  }
  if(missing(RNAMat) || !is(RNAMat, "matrix")){
    stop("Must provide a matrix of RNA data")
  }
  if(missing(genes) || length(genes)!=dim(RNAMat)[1]){
    stop("Must provide a vector of gene names equal to the number of rows in the RNA matrix")
  }
  if(missing(barcodes) || length(barcodes)!=dim(RNAMat)[2]){
    stop("Must provide a vector of barcodes equal to the number of columns in the RNA matrix")
  }
  rownames(RNAMat) = genes
  colnames(RNAMat) = barcodes
  RNASeurat = Seurat::CreateSeuratObject(counts = RNAMat)
  RNASeurat = Seurat::NormalizeData(RNASeurat)
  RNASeurat = Seurat::FindVariableFeatures(RNASeurat)
  RNASeurat = Seurat::ScaleData(RNASeurat)
  RNASeurat = Seurat::RunPCA(RNASeurat)
  RNASeurat = Seurat::FindNeighbors(RNASeurat, dims = dims)
  RNASeurat = Seurat::FindClusters(RNASeurat, resolution = resolution)
  RNASeuratMarkers = Seurat::FindAllMarkers(RNASeurat)
  data.frame(gene=RNASeuratMarkers$gene,
             cluster=RNASeuratMarkers$cluster,
             logFC=RNASeuratMarkers$avg_logFC)
}
