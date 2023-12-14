#' Convert bulk annotation dataframe to a binary matrix
#'
#' \code{convertToMatrix} Reshape bulk annotation dataframe to a region by cluster matrix for permutation
#' @param bulkRegions a dataframe of genomic coordinates of regions. Must contain either four columns
#' with names 'seqnames', 'start', 'end' and 'cluster' or two columns: 'sequence_name' and 'cluster'
#' @return a data.table with sequence name as the first column and clusters as the following columns repressenting whether the sequence is in the cluster
#' @export

convertToMatrix = function(bulkRegions)
{
  if(missing(bulkRegions) | !is(bulkRegions, 'data.frame'))
  {
    stop('Must input genomic coordinates of regions as a dataframe')
  }
  if('Weight' %in% colnames(bulkRegions)) warning('will ignore region weight')

  stopifnot('cluster' %in% colnames(bulkRegions))

  bulkRegions = as.data.table(bulkRegions)
  if(!('sequence_name' %in% colnames(bulkRegions))){
    if(!all(c('seqnames', 'start', 'end') %in% colnames(bulkRegions))){
      stop('regions must have seqnames, start, end')
    }
    bulkRegions[, sequence_name :=  paste0(seqname, ':', start ,'-', end)]
  }
  bulkRegions[, count := 1]

  sum2 = function(x)
  {
    as.numeric(any(x > 0))
  }
  pRE_mat = data.table::dcast(bulkRegions, sequence_name~cluster, value.var = 'count', fun.aggregate = sum2, fill = 0) # max entry = 1 every region + cluster is unique
  pRE_mat = as.data.frame(pRE_mat)
  return(pRE_mat)
}



#' Compute Label Edges
#'
#' \code{computeBulkEdges} generates a matrix of edges from each bulk label (e.g. a group of genomic regions) to each cell by ATACSeq data.
#' The edge weight is proportion of read counts in each set of regions.
#'
#' @param regions a dataframe of genomic coordinates of regions. Must contain either four columns
#' with names 'seqnames', 'start', 'end' and 'cluster' or two columns: 'sequence_name' and 'cluster'. Can have an additional column named 'weight' to weight each region.
#' @param ATAC_Mat a peak by cell data.frame or matrix, colnames are cell barcodes
#' @param peaks a dataframe of genomic coordinates of ATACSeq data. Must contain at least three columns
#'  with names 'seqnames', 'start' and 'end'
#' @param drop whether to drop the labels with no overlapping with regions
#' @return a matrix of edges from each bulk label to each cell
#' @export

computeBulkEdges = function(regions, peaks, ATAC_Mat, drop = T)
{
  if(missing(regions) | !is(regions, 'data.frame'))
  {
    stop('Must input genomic coordinates of regions as a dataframe')
  }
  if(is.null(colnames(regions)) | !('cluster' %in% colnames(regions)))
  {
    stop('regions must have colnames and cluster as a column')
  }
  if(!all(c('seqnames', 'start', 'end') %in% colnames(regions))){
    if('sequence_name' %in% colnames(regions))
    {
      regions = convertSeqnames(regions)
    }else{
    stop('regions must have at least four columns with names seqnames, start, end and cluster')
    }
  }

  if(missing(peaks) | !is(peaks, 'data.frame'))
  {
    stop('Must input genomic coordinates of peaks for the ATACSeq data as a dataframe')
  }
  if(ncol(peaks) < 3 | is.null(colnames(peaks)) | any(colnames(peaks)[1:3]!=c('seqnames', 'start', 'end'))){
    stop('peaks must have at least 3 columns with names seqnames, start and end')
  }
  if(missing(ATAC_Mat) | !(is(ATAC_Mat, 'Matrix') | is(ATAC_Mat, 'matrix') | is(ATAC_Mat, 'data.frame')))
  {
    stop('Must input ATACSeq data as a matrix or dataframe')
  }
  if(is.null(colnames(ATAC_Mat))){
    stop('ATAC_Mat must have cell barcodes as colnames')
  }
  if('weight' %in% colnames(regions))
  {
    labelEnhancers = data.frame('id' = 1:nrow(regions), 'cluster' = regions$cluster, 'weight' = regions$weight)
    labelEnhancers$weight = labelEnhancers$weight/median(labelEnhancers$weight)
  }else{
    labelEnhancers = data.frame('id' = 1:nrow(regions), 'cluster' = regions$cluster)
  }

  regions = as(regions, 'GRanges')
  seqlevelsStyle(regions) <- 'UCSC'
  peaks = as(peaks, 'GRanges')
  seqlevelsStyle(peaks) <- 'UCSC'
  cellPeakCounts = colSums(ATAC_Mat)

  markerOverlaps = GenomicRanges::findOverlaps(peaks, regions)
  hits = unique(markerOverlaps@from)
  ATAC_Mat = as.matrix(ATAC_Mat[hits, ])
  peaks = peaks[hits, ]
  ATACGenePeak = sapply(unique(labelEnhancers$cluster),function(x) {
    ids = labelEnhancers[labelEnhancers$cluster == x, 'id']
    markerOverlaps = GenomicRanges::findOverlaps(peaks, regions[ids,])
    list(peak=markerOverlaps@from,'gene'= ids[markerOverlaps@to])
  })
  if(ncol(labelEnhancers) > 2)
  {
    labelEdges2 <- computeEdgeLabels(labelEnhancers, t(ATAC_Mat), ATACGenePeak, cellPeakCounts, method = 'Weight', drop = drop)
  }else{
    labelEdges2 <- computeEdgeLabels(labelEnhancers, t(ATAC_Mat), ATACGenePeak, cellPeakCounts, drop = drop)
  }

  return(labelEdges2)

}
#' Assign regions defined by bulk sequencing with cell types
#'
#' \code{annotateBulkRegion} Random walk on the cells, cell type labels and bulk labels graph and compute Z-scores from bulk labels to each cell type label.
#' Cell type labels can be organized in a hierarchical structure.
#' @param cellGraph cell-cell graph, a cell-to-cell similarity matrix, row and colnames are cell barcodes
#' @param labelEdges  a matrix or dataframe of edges from each cell type label to each cell. Each row is a cell,
#' each column is a cell type label
#' #' @param labelEdges2  a matrix or dataframe of edges from each bulk label to each cell. Each row is a cell,
#' each column is a bulk label
#' @param groups1 a numeric vector showing groups of cells for permutation between cells and cell type labels. If a single number, 0 means no permutation, non-zero means permutation for all cells.
#' If a vector, must have the same length as the rows of \code{labelEdge}. Only permute edges within the same group of the cells. For cells with entry 0, no permutation.
#' @param groups2 a numeric vector showing groups of regions for permutation between regions and bulk annotations .If a single number, 0 means no permutation, non-zero means permutation for all regions.
#' If a vector, must have the same length as the rows of \code{regionMat}. Only permute edges within the same group of the regions. For regions with entry 0, no permutation.
#' if both \code{groups1} and \code{groups2} are NULL, permute among all cells. If either of them is NULL, assume no permutation for the null object.
#' @param regionMat a data.table or dataframe with sequence name (chrX:XXX-XXX) as the first column and clusters as the following columns. The second to the last columns
#' are binary indicating whether a region belongs to a cluster.
#' @param ATACMatList a list of a peak by cell data.frame or matrix, colnames are cell barcodes. Must provided if permute bulk annotation to regions.
#' @param peaksList a list of a dataframe of genomic coordinates of ATACSeq data. Must contain at least three columns
#'  with names 'seqnames', 'start' and 'end'. Must provided if \code{ATACMatList} is not NULL.
#' @param labelEdgeWeights the edge weight ratio between cell-cell type edges and cell-cell edges and between cell-bulk label edges and cell-cell edges.
#' If it's NULL, will tune edge weights and get the optimal value of it; otherwise will use the input weight
#' @param tr1 if not null, incorporating the cell type hierarchy into the graph. input the cell type tree must be a phylo object
#' @param wtree the edge weights traversing up/down the cell type tree. default: 1
#' @param sampleDepth subsample cells to for faster calculation to tune edge weight ratio, default: 5000
#' @param compute.Zscore whether to do permutation and compute Z-score, If not, will return influence score between labels.
#' @param nround rounds of permutations to compute null distribution. default: 50
#' @param ... other arguments pass to \code{tuneEdgeWeights} and \code{walkCells}
#' @return A cellWalk2 object. if \code{compute.Zscore} is False a list including influence score from bulk labels to cell types and labelEdgeWeights.
#' otherwise a list including influence and Zscore from bulk labels to cell types, as well as labelEdgeWeights.
#' @export
#' @import doParallel
#' @import foreach

annotateBulkRegion = function(cellGraph, labelEdges, labelEdges2, groups1 = 1, groups2 = 0, regionMat = NULL, ATACMatList=NULL, peaksList=NULL, labelEdgeWeights = c(1, 1e-2), tr1 = NULL,
                              wtree = c(1,1), sampleDepth = 5000, compute.Zscore = T, nround = 50, ... )
{
  args = list(...)
  if('steps' %in% names(args)){
    steps = args[['steps']]
  }else{
    steps = Inf
  }
  if('tensorflow' %in% names(args)){
    tensorflow = args[['tensorflow']]
  }else{
    tensorflow = F
  }
  if(missing(labelEdges) || (!is(labelEdges, "data.frame") & (!is(labelEdges, "matrix") & (!is(labelEdges, "Matrix"))))){
    stop("Must provide a dataframe or matrix of cell-to-cell type label edges")
  }
  if(missing(labelEdges2)|| (!is(labelEdges2, "data.frame") & (!is(labelEdges2, "matrix") & (!is(labelEdges2, "Matrix"))))){
    stop("Must provide a dataframe or matrix of cell-to-bulk label edges")
  }
  if(missing(cellGraph) || (!is(cellGraph, "data.frame") & (!is(cellGraph, "matrix") & (!is(cellGraph, "Matrix"))))){
    stop("Must provide a dataframe or matrix of cell-to-cell similarity graph")
  }
  if(is.null(colnames(cellGraph)) || is.null(rownames(cellGraph))){
    stop("cellGraph must have same rownames and colnames as cell barcodes")
  }
  if(any(colnames(cellGraph) != rownames(cellGraph)))
  {
    stop('colname and rowname of cellGraph must be the same')
  }
  if(length(groups1) == 1)
  {
    groups1 = rep(groups1, nrow(labelEdges))
  }

  if(!is.null(tr1))
  {
    stopifnot('tree must have tip labels' = !is.null(tr1$tip.label))
    res = checkLabelEdges(labelEdges, groups1, cellGraph, tips = tr1$tip.label)
  }else{
    res = checkLabelEdges(labelEdges, groups1, cellGraph)
  }
  labelEdges = res[[1]]
  groups1 = res[[2]]

  res = checkLabelEdges(labelEdges2, rep(0, nrow(labelEdges2)), cellGraph)
  labelEdges2 = res[[1]]
  if(any(colnames(labelEdges2) %in% colnames(labelEdges))) stop('bulk labels must be distinct from cell type names')

  diag(cellGraph) = 0

  if(any(groups2!= 0) & (is.null(regionMat) | is.null(ATACMatList) | is.null(peaksList)))
  {
    stop('Must provide bulkRegions, ATACMat and peaks if permute region-label')
  }
  if(!is.null(regionMat) & (!is(regionMat, 'data.frame') & !is(regionMat, 'data.table')))
  {
    stop('Must input genomic coordinates of regions by cluster matrix as a data.frame or data.table')
  }
  if(!is.null(regionMat)){
    if(length(groups2)==1)
    {
      groups2 = rep(groups2, nrow(regionMat))
    }else if(length(groups2) != nrow(regionMat))
    {
      stop('length of groups2 must be the same as rows of regionMat')
    }
  }

  if(!is.null(ATACMatList) & (!is(ATACMatList, 'list') | is.null(peaksList) | !is(peaksList, 'list') | length(peaksList) != length(ATACMatList)))
  {
    stop('Must provide ATACMat and peaks as lists, peaksList must have the same length as ATACMatList')
  }

  if(is.null(labelEdgeWeights))
  {
    message('tunning labelEdgeWeights... ')
    edgeWeights <- tuneEdgeWeights(cellGraph,
                                   list(labelEdges, labelEdges2),
                                   sampleDepth = sampleDepth,
                                   method = "computeLabelHomogeneity",
                                   ...)
    opt_idx = which.max(edgeWeights$cellHomogeneity)
    message('labelHomogeneity at optimal edgeWeight:')
    print(edgeWeights[opt_idx,])
    labelEdgeWeights = unlist(edgeWeights[opt_idx, ])
  }else if(!is(labelEdgeWeights, 'numeric') || length(labelEdgeWeights)!=2){
    stop("Must provide a numeric value weight for each set of cell-to-label edges")
  }

  if(!is.null(tr1))
  {
    if(!is(wtree, "numeric") | length(wtree) !=2 )
    {
      stop('tree edge weights must be a numeric vector with length 2')
    }
    res = tree2Mat(tr1,wtree[1], wtree[2])
    allCellTypes = res[[2]]
    cellTypesM = res[[1]]
  }else{
    allCellTypes = colnames(labelEdges)
    cellTypesM = Matrix(0, length(allCellTypes), length(allCellTypes))
  }
  cellTypesM2 = matrix(0, ncol(labelEdges2), ncol(labelEdges2))
  allCellTypes2 = colnames(labelEdges2)
  labelMatrix = rbind(cbind(cellTypesM, matrix(0, nrow = nrow(cellTypesM), ncol = ncol(cellTypesM2))),
                      cbind(matrix(0, nrow = nrow(cellTypesM2), ncol = ncol(cellTypesM)),cellTypesM2))
  colnames(labelMatrix) = rownames(labelMatrix) = c(allCellTypes, allCellTypes2)
  l1 =  length(allCellTypes)
  l_all = length(allCellTypes) + length(allCellTypes2)

  # construct the whole graph permuted or original
  if(!compute.Zscore) nround = 0
  message('run CellWalker:')
  info1_rand = foreach(r = 0:nround) %dopar%
  {
    message('permutation round ', r)
    ## resample labelEdges
    labelEdges_rand = labelEdges
    if(r > 0 & any(groups1!=0))  # if all groups = 0, no permutation
    {
      for(g in unique(groups1)) # don't include group 0 (these cells don't connect to any of the labels in current set)
      {
        if(g == 0) next
        cells = which(groups1 == g)
        labelEdges_rand[cells, ] = simulate_rand0(labelEdges[cells, ])
      }
    }

    if(r > 0 & any(groups2!=0))  # if all groups = 0, no permutation
    {
      regionMat_rand = regionMat
      for(g in unique(groups2)) # don't include group 0 (some of these cells don't connect to any of the labels in current set)
      {
        if(g == 0) next
        regs = which(groups2 == g)
        regionMat_rand[regs, -1] = simulate_rand_binary0(regionMat[regs, -1]) ## permute
      }
      regionMat_rand[regionMat_rand==0] = NA
      regions = reshape2::melt(regionMat_rand, na.rm=T)
      colnames(regions)[1:2] = c('sequence_name', 'cluster')
      regions = convertSeqnames(regions)
      regions$cluster = as.character(regions$cluster)

      labelEdges2_rand = do.call('rbind', lapply(1:length(ATACMatList), function(i){
        computeBulkEdges(regions, peaksList[[i]], ATACMatList[[i]], drop = F)
      }))
      labelEdges2_rand = checkLabelEdges(labelEdges2_rand, rep(0, nrow(labelEdges2_rand)), cellGraph)[[1]]
      if(!setequal(rownames(labelEdges2), rownames(labelEdges2_rand)) | !setequal(colnames(labelEdges2), colnames(labelEdges2_rand))){
        stop("region to bulk permutation error: cells or bulk labels don't match")
      }
      labelEdges2_rand = labelEdges2_rand[rownames(labelEdges2), colnames(labelEdges2)]
    }else{
      labelEdges2_rand = labelEdges2
    }

    if(!is.null(tr1))
    {
      expandLabelEdges = labelEdgeWeights[1] * cbind(labelEdges_rand, matrix(0,dim(labelEdges_rand)[1], tr1$Nnode)) #add 0s to labelEdges to allow room for internal nodes
    }else{
      expandLabelEdges = labelEdges_rand * labelEdgeWeights[1]
    }
    labelEdges2_rand = labelEdges2_rand * labelEdgeWeights[2]
    cell2label = cbind(expandLabelEdges, labelEdges2_rand)
    combinedGraph = rbind(cbind(labelMatrix,t(cell2label)), cbind(cell2label,cellGraph))

    # compute CellWalk on expanded graph
    infMat <- randomWalk(combinedGraph, tensorflow = tensorflow, steps=steps) #5 steps is usually enough
    aa = infMat[(l1+1):l_all, 1:l1]
    colnames(aa) = allCellTypes
    rownames(aa) = allCellTypes2
    return(as.matrix(aa))
  }

  info1 = info1_rand[[1]]
  if(compute.Zscore) zscore = compute_zscore(info1, info1_rand[2:(nround+1)], nround) else zscore = NULL
  cellWalk = list(infMat=info1, zscore= zscore, labelEdgeWeights = labelEdgeWeights)
  class(cellWalk) = "cellWalk2"
  return(cellWalk)
}
