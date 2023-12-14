#' Compute Cell Similarity by ATACSeq
#'
#' \code{computeCellSim} computes a cell-to-cell similarity matrix using the given method.
#' Currently Jaccard, Cosine and LSI similarity is implemented.
#'
#' @param ATACMat either a cell-by-peak matrix, a SnapATAC object, an ArchR project, or a Cicero count table
#' @param method function name, the method used to compute cell similarity, should operate on rows of a matrix. The default
#' is sparseJaccard, other options are sparseCosine and LSI which computes the euclidean distance in LSI space. LSI only works for input cell-by-peak matrix.
#' @param ... additional parameters used by method
#' @return matrix of cell-to-cell similarity
#' @export
#' @examples
#' data("SampleCellWalkRData")
#' computeCellSim(SampleCellWalkRData$ATACMat)
#'
computeCellSim = function(ATACMat, method=c("sparseJaccard", "sparseCosine", "LSI"), ...){
  if(missing(ATACMat)){
    stop("Must provide either a cell-by-peak matrix or a SnapATAC object")
  }
  # if(!is(method, "function")){
  #   stop("Must provide a similarity function")
  # }

  if(is(ATACMat,"snap")){
    if(!requireNamespace("SnapATAC", quietly = TRUE)){
      stop("Must install SnapATAC")
    }
    if(dim(ATACMat@bmat)[1]==0){
      stop("No bmat")
    }
    stopifnot("Only sparseJaccard and sparseCosine is allowed:" = method %in% c('sparseJaccard', 'sparseCosine'))
    ATACMat = SnapATAC::makeBinary(ATACMat, mat="bmat")
    get(method)(ATACMat@bmat)
  }
  else if(is(ATACMat,"ArchRProject")){
    if(!requireNamespace("ArchR", quietly = TRUE)){
      stop("Must install ArchR")
    }
    if(!"TileMatrix" %in% ArchR::getAvailableMatrices(ATACMat)){
      stop("TileMatrix missing from ArchR project")
    }
    stopifnot("Only sparseJaccard and sparseCosine is allowed:" = method %in% c('sparseJaccard', 'sparseCosine'))
    ATACData = ArchR::getMatrixFromProject(ATACMat, "TileMatrix", verbose=FALSE, binarize = TRUE)
    ATACMat = Matrix::t(assay(ATACData))
    get(method)(ATACMat)
  }
  else if(is(ATACMat,"Matrix")|is(ATACMat,"matrix")){
    if(min(dim(ATACMat))==0){
      stop("Not enough cells or peaks in matrix")
    }
    if(dim(ATACMat)[1]>dim(ATACMat)[2]){
      warning("ATACMat has more cells than peaks, does it need to be transposed?")
    }
    if(method %in% c('sparseJaccard', 'sparseCosine'))
    {
      get(method)(ATACMat>0)
    }else if(method == 'LSI'){
      get(method)(t(ATACMat))
    }else{
      stop('method not implemented yet')
    }
  }
  else if(is(ATACMat, "data.frame")){
      stopifnot("Only sparseJaccard and sparseCosine is allowed:" = method %in% c('sparseJaccard', 'sparseCosine'))
      peaks = unique(ATACMat[,1])
      cells = unique(ATACMat[,2])
      rowIndex = match(ATACMat[,2], cells)
      colIndex = match(ATACMat[,1], peaks)
      ATACMat = Matrix::sparseMatrix(rowIndex, colIndex, x=ATACMat[,3])
      simMat = method(ATACMat>0)
      dimnames(simMat) = list(cells,cells)
      simMat
  }else{
    stop("Must provide either a cell-by-peak matrix, SnapATAC object, ArchR project, or Cicero count table")
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
#' @examples
#' getRegions()
#'
getRegions = function(geneBody=TRUE, genome="hg38", names="Entrez"){
  if(genome=="hg38"){
    if(names=="Ensembl"){
      stop("Not an implemented gene name identifier")
    }
    else if(names=="Entrez"){
      pathToGenes <- system.file("extdata", "Hg38_knownGene.rds", package = "CellWalkR")
      genes <- readRDS(pathToGenes)
    }
    else{
      stop("Not a valid gene name identifier")
    }
  }
  else if(genome=="mm10"){
    if(names=="Ensembl"){
      pathToGenes <- system.file("extdata", "Mm10_ensGene.rds", package = "CellWalkR")
      genes <- readRDS(pathToGenes)
    }
    else if(names=="Entrez"){
      pathToGenes <- system.file("extdata", "Mm10_knownGene.rds", package = "CellWalkR")
      genes <- readRDS(pathToGenes)
    }
    else{
      stop("Not a valid gene name identifier")
    }
  }
  else{
    stop("Not an implemented genome build")
  }

  regions = GenomicRanges::promoters(genes)
  if(geneBody){
    regions = c(regions, genes)
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
#' @examples
#' data("SampleCellWalkRData")
#' regions <- getRegions()
#' mapPeaksToGenes(SampleCellWalkRData$labelGenes,
#'                 SampleCellWalkRData$ATACMat,
#'                 SampleCellWalkRData$peaks,
#'                 regions)
#'
mapPeaksToGenes = function(labelGenes, ATACMat, peaks, regions){
  if(missing(labelGenes)){
    stop("Must provide a table with genes of interest in first column and corresponding labels in second column")
  }
  if(is.null(dim(labelGenes)) || dim(labelGenes)[2]<2){
    stop("Must provide a table with genes of interest in first column and corresponding labels in second column")
  }
  labels = as.character(unique(labelGenes[,2]))
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


#' Map SnapATAC peaks to genes
#'
#' \code{mapSnapATACToGenes} Generates a mapping from genes to peaks per label
#'
#' @param labelGenes data.frame with genes of interest in first column and corresponding
#' labels in second column
#' @param snap a snap object
#' @param whichMat string determing whether to use "bmat" or "gmat" from snap object
#' @param regions GRanges object of genomic regions associated with genes (only needed for "bmat")
#' @return list of genes and corresponding peaks by label
#' @export
#'
mapSnapATACToGenes = function(labelGenes, snap, whichMat = "bmat", regions){
  if(missing(snap)){
    stop("Must provide a SnapATAC object")
  }
  if(!is(snap,"snap")){
    stop("Must provide a SnapATAC object")
  }
  if(!requireNamespace("SnapATAC", quietly = TRUE)){
    stop("Must install SnapATAC")
  }
  if(whichMat=="bmat"){
    mapPeaksToGenes(labelGenes, ATACMat = snap@bmat, peaks=snap@feature, regions)
  }
  else if(whichMat=="gmat"){
    whichGenes = match(unique(labelGenes[,1]),colnames(snap@gmat))
    if(length(which(!is.na(whichGenes)))<10){
      warning("Fewer than 10 gene names match")
    }
    sapply(as.character(unique(labelGenes[,2])), function(x) {
      markerGenes = labelGenes[labelGenes[,2]==x,1]
      markerGenes = markerGenes[markerGenes %in% colnames(snap@gmat)]
      list(peak=match(markerGenes,colnames(snap@gmat)),
           gene=markerGenes)
    })
  }
  else{
    stop("Must use bmat or gmat")
  }
}

#' Map ArchR peaks to genes
#'
#' \code{mapArchRToGenes} Generates a mapping from genes to peaks per label
#'
#' @param labelGenes data.frame with genes of interest in first column and corresponding
#' labels in second column
#' @param ArchRproj an ArchR project
#' @param whichMat string determing whether to use "TileMatrix" or "GeneScoreMatrix" from ArchR project
#' @param regions GRanges object of genomic regions associated with genes (only needed for "TileMatrix")
#' @return list of genes and corresponding peaks by label
#' @export
#'
mapArchRToGenes = function(labelGenes, ArchRproj, whichMat = "TileMatrix", regions){
  if(missing(ArchRproj)){
    stop("Must provide an ArchR project")
  }
  if(!is(ArchRproj,"ArchRProject")){
    stop("Must provide an ArchR project")
  }
  if(!requireNamespace("ArchR", quietly = TRUE)){
    stop("Must install ArchR")
  }
  if(whichMat=="TileMatrix"){
    if(!"TileMatrix" %in% ArchR::getAvailableMatrices(ArchRproj)){
      stop("TileMatrix missing from ArchR project")
    }
    ATACData = ArchR::getMatrixFromProject(ArchRproj, "TileMatrix", verbose=FALSE, binarize = TRUE)
    peaks = GenomicRanges::GRanges(ATACData@elementMetadata$seqnames, IRanges(start = ATACData@elementMetadata$start, width = 500))
    ATACMat = Matrix::t(assay(ATACData))
    mapPeaksToGenes(labelGenes, ATACMat = ATACMat, peaks=peaks, regions)
  }
  else if(whichMat=="GeneScoreMatrix"){
    if(!"GeneScoreMatrix" %in% ArchR::getAvailableMatrices(ArchRproj)){
      stop("GeneScoreMatrix missing from ArchR project")
    }
    ATACData = ArchR::getMatrixFromProject(ArchRproj, "GeneScoreMatrix", verbose=FALSE)
    whichGenes = match(unique(labelGenes[,1]),rowData(ATACData)$name)
    if(length(which(!is.na(whichGenes)))<10){
      warning("Fewer than 10 gene names match")
    }
    sapply(as.character(unique(labelGenes[,2])), function(x) {
      markerGenes = labelGenes[labelGenes[,2]==x,1]
      markerGenes = markerGenes[markerGenes %in% rowData(ATACData)$name]
      list(peak=match(markerGenes,rowData(ATACData)$name),
           gene=markerGenes)
    })
  }
  else{
    stop("Must use TileMatrix or GeneScoreMatrix")
  }
}

#' Map peaks to genes based on Cicero coaccessibility
#'
#' \code{mapCiceroToGenes} Generates a mapping from genes to peaks per label
#'
#' @param labelGenes data.frame with genes of interest in first column and corresponding
#' labels in second column
#' @param cicero_gene_activities gene activity matrix
#' @return list of genes and corresponding peaks by label
#' @export
#'
mapCiceroToGenes = function(labelGenes, cicero_gene_activities){
  if(missing(cicero_gene_activities)){
    stop("Must provide a gene activity matrix")
  }
  if(!is(cicero_gene_activities,"Matrix")){
    stop("Must provide a gene activity matrix")
  }
  whichGenes = match(unique(labelGenes[,1]),rownames(cicero_gene_activities))
  if(length(which(!is.na(whichGenes)))<10){
    warning("Fewer than 10 gene names match")
  }
  sapply(as.character(unique(labelGenes[,2])), function(x) {
    markerGenes = labelGenes[labelGenes[,2]==x,1]
    markerGenes = markerGenes[markerGenes %in% rownames(cicero_gene_activities)]
    list(peak=match(markerGenes,rownames(cicero_gene_activities)),
         gene=markerGenes)
  })
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
#' @param ATACMat either a cell-by-peak matrix, a SnapATAC object, or an ArchR project
#' @param ATACGenePeak per label mapping of peaks to genes returned by mapPeaksToGenes()
#' @param method scaling method, either "Expression","Correlation", or "None"
#' @param filters list of GRanges lists for locations that pass filters
#' @param filterWeights numeric vector of weights assigned to each filter
#' @param filterOut boolean for each filter of whether those regions are permitted or filtered out
#' @param filterGene boolean for each filter of whether it applies to genes as a whole or just overlapping peaks (note that weights do not apply to peaks)
#' @param regions regions map peaks to genes
#' @param peaks GRanges of peaks in ATACMat
#' @param whichMat string determing whether to use "bmat" or "gmat" from snap object or "TileMatrix" or "GeneScoreMatrix" from ArchR project
#' @return matrix of weights from each label to each cell
#' @export
#' @examples
#' data("SampleCellWalkRData")
#' computeLabelEdges(SampleCellWalkRData$labelGenes,
#'                   SampleCellWalkRData$ATACMat,
#'                   SampleCellWalkRData$ATACGenePeak)
#'
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
  peaks,
  whichMat="bmat"){
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
    stop("Must provide either a cell-by-peak matrix, a SnapATAC object, or an ArchR project")
  }
  else{
    if(is(ATACMat,"snap")){
      if(whichMat == "bmat"){
        peaks = ATACMat@feature
        ATACMat = ATACMat@bmat
      }
      else if(whichMat == "gmat"){
        ATACMat = ATACMat@gmat
      }
      else{
        stop("Must use bmat or gmat")
      }
    }
    if(is(ATACMat,"ArchRProject")){
      if(whichMat == "TileMatrix"){
        ATACData = ArchR::getMatrixFromProject(ATACMat, "TileMatrix", verbose=FALSE, binarize = TRUE)
        peaks = GenomicRanges::GRanges(ATACData@elementMetadata$seqnames, IRanges(start = ATACData@elementMetadata$start, width = 500))
        ATACMat = Matrix::t(assay(ATACData))
      }
      else if(whichMat == "GeneScoreMatrix"){
        ATACData = ArchR::getMatrixFromProject(ATACMat, "GeneScoreMatrix", verbose=FALSE)
        ATACMat = Matrix::t(assay(ATACData))
      }
      else{
        stop("Must use TileMatrix or GeneScoreMatrix")
      }
    }
    if(!is(ATACMat,"Matrix") & !is(ATACMat,"matrix")){
      stop("Must provide either a cell-by-peak matrix")
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
        regionFilterOverlap = GenomicRanges::countOverlaps(regions, filterRanges)
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
        # geneExp = labelGenes[match(genesInMarkers,labelGenes[labelGenes[,2]==label,1]),3]
        geneExp = labelGenes[labelGenes[,2]==label,3][match(genesInMarkers,labelGenes[labelGenes[,2]==label,1])]
      }
      else{
        stop("Not a recognized scaling method")
      }

      # geneFilterScale = geneFilter[match(genesInMarkers,labelGenes[labelGenes[,2]==label,1])]
      geneFilterScale = geneFilter[labelGenes[,2]==label][match(genesInMarkers,labelGenes[labelGenes[,2]==label,1])]

      if(!missing(filters) && sum(!filterGene)>0){
        if(missing(peaks) || !is(peaks,"GRanges")){
          stop("Must provide a GRanges object of peaks to apply filters to peaks")
        }
        if(!length(peaks)==dim(ATACMat)[2]){
          stop("Number of peaks must be the same as the number of columns in ATACMat")
        }

        for(filterIndex in which(!filterGene)){
          filterRanges = filters[[filterIndex]]
          peakFilterOverlap = GenomicRanges::countOverlaps(peaks, filterRanges)
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
    cellsInMarkers = cellsInMarkers[,colSums(cellsInMarkers)!=0]
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
#' @param steps integer indicating number of steps to take if walk should not be run to convergence
#' @param tensorflow boolean to indicate whether to compute on GPU
#' @return cellWalk object with influence matrix and labels
#' @export
#' @examples
#' data("SampleCellWalkRData")
#' labelEdgesList <- list(SampleCellWalkRData$labelEdges)
#' cellWalk <- walkCells(SampleCellWalkRData$cellEdges, labelEdgesList, 1)
#'
walkCells = function(cellEdges, labelEdgesList, labelEdgeWeights, sampleDepth, steps = Inf, tensorflow=FALSE){
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
  infMat = randomWalk(combinedGraph, tensorflow=tensorflow, steps=steps)

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
#' @param steps integer indicating number of steps to take if walk should not be run to convergence
#' @param cellTypes list of names of cell types, if not provided unique cell labels are used
#' @param parallel boolean, execute in parallel
#' @param numCores integer, number of cores to use for parallel execution
#' @param trackProgress boolean, print percent run completed
#' @param tensorflow boolean to indicate whether to compute on GPU
#' @param method criterion to choose the optimal parameters. If annotate cells or map cell types, use computeCellHomogeneity otherwise use  computeLabelHomogeneity
#' @return data frame of weights with corresponding cell homogeneity
#' @export
#' @examples
#' data("SampleCellWalkRData")
#' labelEdgesList <- list(SampleCellWalkRData$labelEdges)
#' tuneEdgeWeights(SampleCellWalkRData$cellEdges, labelEdgesList, 10^seq(1,7,1))
#'
tuneEdgeWeights = function(
  cellEdges,
  labelEdgesList,
  labelEdgeOpts = 10^seq(-4,4,1),
  sampleDepth = 10000,
  steps = Inf,
  cellTypes,
  parallel = FALSE,
  numCores = 1,
  trackProgress = FALSE,
  tensorflow=FALSE,
  method="computeCellHomogeneity"){

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
    if(tensorflow){
      warning("Cant run both on GPU and in parallel, running in parallel")
    }

    cellHomogeneity = parallel::mcmapply(function(param){
      theseWeights = as.numeric(parameterCombos[param,])
      cellWalk = walkCells(cellEdges, labelEdgesList, theseWeights, sampleDepth, steps)

      if(cellTypesMissing){
        cellTypes = unique(cellWalk[["cellLabels"]])
      }
      if(method == 'computeCellHomogeneity')
      {
        get(method)(cellWalk, cellTypes)
      }else if(method == 'computeLabelHomogeneity'){
        get(method)(cellWalk, colnames(labelEdgesList[[1]]), colnames(labelEdgesList[[2]]))
      }else{
        stop(method, 'not implemented')
      }

    }, param=1:dim(parameterCombos)[1], mc.cores = numCores)
  }
  else{
    cellHomogeneity = c()
    for(param in 1:dim(parameterCombos)[1]){
      theseWeights = as.numeric(parameterCombos[param,])

      cellWalk = walkCells(cellEdges, labelEdgesList, theseWeights, sampleDepth, steps, tensorflow=tensorflow)

      if(cellTypesMissing){
        cellTypes = unique(cellWalk[["cellLabels"]])
      }

      if(method == 'computeCellHomogeneity')
      {
        score = get(method)(cellWalk, cellTypes)
      }else if(method == 'computeLabelHomogeneity'){
        score = get(method)(cellWalk, colnames(labelEdgesList[[1]]), colnames(labelEdgesList[[2]]))
      }else{
        stop(method, 'not implemented')
      }

      cellHomogeneity = c(cellHomogeneity, score)

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
#' @param steps integer indicating number of steps to take if walk should not be run to convergence
#' @param cellTypes list of names of cell types, if not provided unique cell labels are used
#' @param parallel boolean, execute in parallel
#' @param numCores integer, number of cores to use for parallel execution
#' @param trackProgress boolean, print percent run completed
#' @param tensorflow boolean to indicate whether to compute on GPU
#' @return data frame of parameters with corresponding cell homogeneity
#' @export
#' @examples
#' data("SampleCellWalkRData")
#' regions <- getRegions()
#' filters <- list(SampleCellWalkRData$filter)
#' labelGenesList <- list(SampleCellWalkRData$labelGenes)
#' labelEdgesList <- list(SampleCellWalkRData$labelEdges)
#' \dontrun{tuneFilterWeights(SampleCellWalkRData$cellEdges,
#'                   labelGenesList,
#'                   labelEdgesList,
#'                   1,
#'                   SampleCellWalkRData$ATACMat,
#'                   SampleCellWalkRData$ATACGenePeak,
#'                   filters=filters,
#'                   regions=regions)}
#'
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
  steps,
  cellTypes,
  parallel = FALSE,
  numCores = 1,
  trackProgress = FALSE,
  tensorflow = FALSE){

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
    if(tensorflow){
      warning("Cant run both on GPU and in parallel, running in parallel")
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

      cellWalk = walkCells(cellEdges, labelEdgesList, labelEdgeWeights, sampleDepth, steps)

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

      cellWalk = walkCells(cellEdges, labelEdgesList, labelEdgeWeights, sampleDepth, steps, tensorflow=tensorflow)
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
#' @param RNAMat gene-by-barcode raw count matrix of scRNASeq data
#' @param genes character vector of gene names
#' @param barcodes character vector of cell barcodes
#' @param dims integer vector of PCA dimensions to use
#' @param resolution numeric resolution to use in Louvain clustering
#' @param only.pos only return positive markers
#' @param min.cells used by Seurat, include features detected in at least this many cells
#' @param scale.factor used by Seurat, scale factor for cell-level normalization
#' @param nfeatures used by Seurat, Number of features to select as top variable features
#' @return table of marker genes, labels, and logFC in expression
#' @export
#'
findMarkers = function(RNAMat, genes, barcodes, dims=1:20, resolution=0.5, only.pos = FALSE, min.cells = 3, scale.factor = 1e4, nfeatures = 3000){
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
  RNASeurat = Seurat::CreateSeuratObject(counts = RNAMat, min.cells = min.cells)
  RNASeurat = Seurat::NormalizeData(RNASeurat, scale.factor = scale.factor)
  RNASeurat = Seurat::FindVariableFeatures(RNASeurat, nfeatures = nfeatures)
  RNASeurat = Seurat::ScaleData(RNASeurat, verbose = F)
  RNASeurat = Seurat::RunPCA(RNASeurat)
  RNASeurat = Seurat::FindNeighbors(RNASeurat, dims = dims)
  RNASeurat = Seurat::FindClusters(RNASeurat, resolution = resolution)
  RNASeuratMarkers = Seurat::FindAllMarkers(RNASeurat, only.pos = only.pos)
  data.frame(gene=RNASeuratMarkers$gene,
             cluster=as.character(RNASeuratMarkers$cluster),
             avg_logFC=RNASeuratMarkers$avg_logFC)
}

#' preprocess scRNAseq data
#'
#' \code{processRNASeq} Uses Seurat to normalize raw count matrix, generate cell-to-cell similarity graph and find markers in scRNA-seq data and returns
#' a table of markers and labels, cell-cell graph that can be used by CellWalker.
#'
#' @param RNAMat gene-by-barcode raw count matrix of scRNASeq data, either a matrix or a data.frame
#' @param meta.data Additional cell-level metadata. Should be a data.frame where the rows are cell names matched with the column names of the counts matrix or \code{barcodes}
#' and the columns are additional metadata fields.
#' @param group.col the column name indicates the groups/cell clusters. If provided, will find markers for each group; otherwise, will use Seurat to identify clusters first.
#' @param genes character vector of gene names. If NULL, use rownames of \code{RNAMat} as genes.
#' @param barcodes character vector of cell barcodes
#' @param do.findMarkers Whether to find cell type markers. Default: True.
#' @param computeKNN Whether to compute cell-cell graph using KNN. Default: True.
#' @param computeSimilarity Whether to compute cell-cell similarity matrix based on PCs. Default: False.
#' @param knn Defines k for the k-nearest neighbor algorithm
#' @param buildTree Whether to construct a cell type tree. Default: True
#' @param dims integer vector of PCA dimensions to use
#' @param resolution numeric resolution to use in Louvain clustering
#' @param only.pos only return positive markers
#' @param min.cells used by Seurat, include features detected in at least this many cells.
#' @param scale.factor used by Seurat, scale factor for cell-level normalization
#' @param nfeatures used by Seurat, Number of features to select as top variable features
#' @param ... other parameters can be passed to \code{CreateSeuratObject}
#' @return normalized gene expression matrix,  table of marker genes, labels, and logFC in expression,
#' a matrix of cell-to-cell similarity matrix, a cell type tree and meta.data
#' @export
#'
processRNASeq = function(RNAMat, meta.data = NULL, group.col = NULL, do.findMarkers = TRUE, computeKNN = TRUE, computeSimilarity = FALSE, knn = 20, buildTree = TRUE,
                         genes=NULL, barcodes=NULL, dims=1:20, resolution=0.5, only.pos = FALSE, min.cells = 3, scale.factor = 1e4,
                         nfeatures = 3000, ...){
  if(!requireNamespace("Seurat", quietly = TRUE)){
    stop("Must install Seurat")
  }
  if(missing(RNAMat) || (!is(RNAMat, "matrix") & !is(RNAMat, "Matrix") & !is(RNAMat, "data.frame"))){
    stop("Must provide a matrix or a dataframe of RNA data")
  }
  if(is.null(genes) & is.null(rownames(RNAMat))){
    stop('Must provide a vector of gene names')
  }
  if(is.null(barcodes) & is.null(colnames(RNAMat))){
    stop('Must provide a vector of barcodes or cell names')
  }
  if(!is.null(genes) & length(genes)!=dim(RNAMat)[1]){
    stop("Must provide a vector of gene names equal to the number of rows in the RNA matrix")
  }
  if(!is.null(barcodes) & length(barcodes)!=dim(RNAMat)[2]){
    stop("Must provide a vector of barcodes equal to the number of columns in the RNA matrix")
  }
  if(!is.null(genes)) rownames(RNAMat) = genes
  if(!is.null(barcodes)) colnames(RNAMat) = barcodes

  RNASeuratMarkers  = cellgraph = tr = NULL

  RNASeurat = Seurat::CreateSeuratObject(counts = RNAMat, meta.data = meta.data, min.cells = min.cells, ...)
  RNASeurat = Seurat::NormalizeData(RNASeurat, scale.factor = scale.factor)
  RNASeurat = Seurat::FindVariableFeatures(RNASeurat, nfeatures = nfeatures)
  RNASeurat = Seurat::ScaleData(RNASeurat, features = rownames(RNAMat))
  RNASeurat = Seurat::RunPCA(RNASeurat, verbose = F)

  if(computeKNN) {
    RNASeurat = Seurat::FindNeighbors(RNASeurat, dims = dims, k.param = knn)
    cellgraph = RNASeurat@graphs$RNA_snn
  }else if(computeSimilarity){
    distance <- dist(RNASeurat@reductions$pca@cell.embeddings[,dims]) #include the weight of each PC
    distance = distance/max(distance)
    cellgraph = as.matrix(1 - distance) # similarity
    diag(cellgraph) = 1
  }

  if(do.findMarkers | buildTree) # clustering if needed
  {
    if(is.null(group.col) | is.null(meta.data))
    {
      if(!computeKNN) RNASeurat = Seurat::FindNeighbors(RNASeurat, dims = dims, k.param = knn)
      RNASeurat = Seurat::FindClusters(RNASeurat, resolution = resolution)
    }else{
      stopifnot("meta.data must contain group.col as column" = group.col %in% colnames(meta.data))
      Seurat::Idents(RNASeurat) = group.col
    }
  }

  if(do.findMarkers) {
    RNASeuratMarkers = Seurat::FindAllMarkers(RNASeurat, only.pos = only.pos)
  }

  if(buildTree){
    RNASeurat =  Seurat::BuildClusterTree(RNASeurat)
    tr =  Seurat::Tool(object = RNASeurat, slot = 'Seurat::BuildClusterTree')
  }

  return(list("expr_norm" = RNASeurat@assays$RNA@scale.data, "markers" = RNASeuratMarkers, "cellGraph" = cellgraph, 'tr' = tr, 'meta' = RNASeurat@meta.data))
}

#' merge scRNAseq data
#'
#' \code{mergeRNASeq} Uses Seurat to normalize raw count matrix of each data set, integrate and
#' generate cell-to-cell similarity graph that can be used by CellWalker.
#'
#' @param RNAMatList a list of gene-by-barcode raw count matrix of scRNASeq data, either a matrix or a data.frame.
#' Each of them has cell barcodes as colnames and genes as rownames.
#' @param integrate Whether to use Seurat to remove batch effects and integrate all datasets provided. Default: True
#' @param min.cells used by Seurat, include features detected in at least this many cells.
#' @param scale.factor used by Seurat, scale factor for cell-level normalization
#' @param nfeatures used by Seurat, Number of features to select as top variable features and
#' integration feautres (if \code{integrate} is TRUE)
#' @param ndim integer vector of CCA dimensions to use for integration or PCA dimensions
#' @param computeKNN Whether to compute cell-cell graph using KNN. Default: True. If false, compute cell-cell similarity matrix based on PCs.
#' @param knn Defines k for the k-nearest neighbor algorithm
#' @param ... other parameters can be passed to \code{CreateSeuratObject}
#' @return normalized gene expression matrix and a matrix of cell-to-cell similarity matrix.
#' @export
#'
mergeRNASeq = function(RNAMatList, integrate = TRUE, min.cells = 3, scale.factor = 1e4, nfeatures = 3000, ndim = 30, computeKNN = TRUE, knn = 20, ...)
{
  if(missing(RNAMatList) || !is(RNAMatList, "list")){
    stop("Must provide a list of count matrices")
  }
  # normalize and identify variable features for each dataset independently
  ifnb.list <- lapply(X = RNAMatList, FUN = function(x) {
    if(!is(x, "matrix") & !is(x, "Matrix") & !is(x, "data.frame")){
      stop("Must provide a matrix or a dataframe for each RNA data")
    }
    if(is.null(rownames(x)) | is.null(colnames(x))){
      stop("Each RNA data must have gene names as rownames and cell names as colnames")
    }
    x= Seurat::CreateSeuratObject(counts = x, min.cells = min.cells, ...)
    x = Seurat::NormalizeData(x, scale.factor = scale.factor)
    x = Seurat::FindVariableFeatures(x, nfeatures = nfeatures)
  })

  if(integrate)
  {
    features <- Seurat::SelectIntegrationFeatures(object.list = ifnb.list,nfeatures = nfeatures)
    ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
      x <- Seurat::ScaleData(x, features = features, verbose = FALSE)
      x <- Seurat::RunPCA(x, features = features, verbose = FALSE)
    })
    #
    ## find anchors and integrate
    anchors <- Seurat::FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "LogNormalize",
                                             anchor.features = features, dims = 1:ndim)
    RNASeurat <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "LogNormalize", dims = 1:ndim)
    Seurat::DefaultAssay(RNASeurat) <- "integrated"
  }else{
    RNASeurat <- merge(ifnb.list[[1]], ifnb.list[[-1]]) #add.cell.ids = paste0("_", 1:length(ifnb.list))
    RNASeurat <- Seurat::FindVariableFeatures(RNASeurat, selection.method = "vst", nfeatures = nfeatures)
  }
  RNASeurat <- Seurat::ScaleData(RNASeurat, verbose = FALSE)
  RNASeurat <- Seurat::RunPCA(RNASeurat, verbose = F)
  if(computeKNN)
  {
    RNASeurat <- Seurat::FindNeighbors(RNASeurat, dims = 1:ndim, verbose = F, k.param = knn)
    if(integrate) cellgraph = RNASeurat@graphs$integrated_snn else cellgraph = RNASeurat@graphs$RNA_snn
  }else{
    distance <- dist(RNASeurat@reductions$pca@cell.embeddings[,1:ndim]) #include the weight of each PC
    distance = distance/max(distance)
    cellgraph = as.matrix(1 - distance) # similarity
    diag(cellgraph) = 1
  }
  if(integrate)
  {
    return(list("expr_norm" = RNASeurat@assays$integrated@scale.data, "cellGraph" = cellgraph))
  }else{
    return(list("expr_norm" = RNASeurat@assays$RNA@scale.data, "cellGraph" = cellgraph))
  }
}

#' Download Sample Data
#'
#' \code{downloadSampleData} Download sample data for use with vignette
#'
#' @param dest folder to download data to
#' @return destination folder
#' @export
#'
downloadSampleData = function(dest="inst/extdata"){
  if(!dir.exists(dest)){
    dir.create(dest, showWarnings = FALSE)
  }
  downloadFiles = download.file(url = "https://figshare.com/ndownloader/files/32630126", destfile = file.path(dest,"SamplePeakMat.mtx"))
  downloadFiles = download.file(url = "https://figshare.com/ndownloader/files/32670425", destfile = file.path(dest,"SamplePeaks.txt"))
  dest
}

#' Load Sample Data
#'
#' \code{loadSampleData} Download sample data for use with vignette, only used for testing.
#' To load sampel data, simply call \code{data("SampleCellWalkRData")}
#' @export
#' @return sample data
#' @examples
#' loadSampleData()
#'
loadSampleData = function(){
  # data("SampleCellWalkRData")
  # SampleCellWalkRData
}
