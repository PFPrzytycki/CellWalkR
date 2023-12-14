#' Compute Label Edges
#'
#' \code{computeTypeEdges} generates a matrix of edges from each cell type label to each cell from gene expression.
#' The edge weight is normalized gene expression of markers weighted by the log2FC.
#'
#' @param exprMat_norm a gene by cell data.frame or matrix, rownames are genes, colnames are cell barcodes
#' @param markers marker genes for each cell type, columns are gene, cluster, p_val_adj (optional) and avg_log2FC (optional).
#' @param pval.cutoff select markers with adjust pvalue (\code{p_val_adj}) < pval.cutoff. Default: 0.05
#' @param log2FC.cutoff select markers with abs(log2FC.cutoff)> 0.5 if only.pos = F or log2FC.cutoff>0.5 if only.pos = T. Default: 0.5
#' @param only.pos only include positive markers
#' @return a matrix of edges from each cell type label to each cell
#' @export
#' @import data.table

computeTypeEdges <- function(exprMat_norm, markers, pval.cutoff = 0.05, log2FC.cutoff = 0.5, only.pos = F)
{
  if(!requireNamespace("data.table", quietly = TRUE)){
    stop("Must install data.table")
  }
  if(missing(exprMat_norm) || (!is(exprMat_norm, "data.frame") & !is(exprMat_norm, "matrix")  & !is(exprMat_norm, "Matrix"))){
    stop("Must provide a dataframe or matrix of RNA data")
  }
  if(is.null(colnames(exprMat_norm)) | is.null(rownames(exprMat_norm)))
  {
    stop('exprMat_norm must have column and row names')
  }
  if(missing(markers) || !is(markers, "data.frame")){
    stop("Must provide a dataframe of markers")
  }
  if(is.null(markers$gene) || is.null(markers$cluster)){
    stop("markers must have 'gene' and 'cluster' columns")
  }


  if(is.null(markers$avg_log2FC)) markers$avg_log2FC = log2FC.cutoff
  if(is.null(markers$p_val_adj)) markers$p_val_adj =  pval.cutoff

  markers = data.table(markers)
  if(only.pos)
  {
    markers = markers[avg_log2FC >= log2FC.cutoff & p_val_adj <= pval.cutoff]
  }else{
    markers = markers[abs(avg_log2FC) >= log2FC.cutoff & p_val_adj <= pval.cutoff]
  }

  markers_inter = markers[markers$gene %in% rownames(exprMat_norm)]

  if(length(unique(markers_inter$cluster)) < length(unique(markers$cluster)))
  {
    stop("Some clusters don't have markers. Need to lower the pval/logFC thresholds or check if expression matrix contains the markers")
  }

  labelEdges = markers_inter[, list("score" = colSums(exprMat_norm[.SD[['gene']],]* .SD[['avg_log2FC']]) / sum(abs(.SD[['avg_log2FC']])), "cell" = colnames(exprMat_norm)),
                             by = cluster]
  labelEdges = reshape2::acast(labelEdges, cell~cluster, value.var = 'score')
  labelEdges[labelEdges <0] = 0

  labelEdges
}



#' Annotate cells by cell types
#' \code{annotateCells} Random walk on the cells and cell type labels graph and compute influence scores from each cell to each cell type labels.
#' Cell type labels can be organized in a hierarchical structure.
#' @param cellGraph cell-cell graph, a cell-to-cell similarity matrix, row and colnames are cell barcodes
#' @param labelEdges  a matrix or dataframe of edges from each cell type label to each cell. Each row is a cell,
#' each column is a cell type label
#' @param weight1 the edge weight ratio between cell-label edges and cell-cell edges.
#' If it's NULL, will tune edge weights and get the optimal value of it; otherwise will use the input weight
#' @param tr1 if not null, incorporating the cell type hierarchy into the graph. input the cell type tree as a phylo object
#' @param wtree the edge weight between cell type labels, default: 1
#' @param sampleDepth subsample cells to for faster calculation to tune edge weight ratio, default: 3000
#' @param ... other arguments pass to \code{tuneEdgeWeights} and \code{walkCells}
#' @return a list including a CellWalker object and weight1
#' @export

annotateCells <- function(cellGraph, labelEdges, weight1 = NULL, sampleDepth =3000, ... , tr1 = NULL, wtree = 1)
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
    stop("Must provide a dataframe or matrix of cell-to-label edges")
  }
  if(is.null(colnames(labelEdges)) || is.null(rownames(labelEdges))){
    stop("labelEdges must have cell barcodes as rownames and cell type labels as colnames")
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
  if(!setequal(colnames(cellGraph), rownames(labelEdges)))
  {
    stop('labelEdges and cellGraph needs to have the same set of cell barcodes')
  }
  if(!is.null(tr1) & !is(tr1, 'phylo'))
  {
    stop('tree must be a phylo object or NULL')
  }

  labelEdges = labelEdges[rownames(cellGraph),]
  diag(cellGraph) = 0

  if(is.null(weight1))
  {
    labelEdgesList <- list(labelEdges)
    edgeWeights <- tuneEdgeWeights(cellGraph,
                                   labelEdgesList,
                                   sampleDepth = sampleDepth,
                                   ...)
    opt_idx = which.max(edgeWeights$cellHomogeneity)
    message('cellHomogeneity at each edgeWeight:')
    print(edgeWeights)
    weight1 = edgeWeights[opt_idx, 1]
  }

  if(is.null(tr1))
  {
    cellWalk <- walkCells(cellGraph,
                          labelEdgesList,
                          labelEdgeWeights = weight1,
                          steps=steps, tensorflow = tensorflow)

    return(list(cellWalk, weight1))
  }


  res1 = tree2Mat(tr1)
  allCellTypes = res1[[2]]
  cellTypesM = res1[[1]]
  colnames(cellTypesM) = rownames(cellTypesM) = allCellTypes

  l = (length(allCellTypes) + 1)/2 #vector of cell type names
  l_all = length(allCellTypes) #vector with cell type names and internal node names
  labelEdges = labelEdges[, allCellTypes[1:l]]
  expandLabelEdges = cbind(labelEdges, matrix(0,dim(labelEdges)[1],l_all-l)) #add 0s to labelEdges to allow room for internal nodes
  labelMatrix = cellTypesM
  cell2label = weight1*expandLabelEdges
  combinedGraph = rbind(cbind(wtree*labelMatrix,t(cell2label)), #combined graph w/ internal nodes
                        cbind(cell2label,cellGraph))

  infMat <- randomWalk(combinedGraph, tensorflow = tensorflow, steps = steps) #5 steps is usually enough
  normMat <- normalizeInfluence(infMat[-(1:l_all),1:l])
  colnames(normMat) <- allCellTypes[1:l]
  rownames(normMat) <- rownames(cellGraph)
  cellLabels <- apply(normMat, 1, function(x) {
    if(max(x) < 0) return(NA)
    colnames(normMat)[order(x, decreasing = TRUE)][1]
  }) #top scoring label
  cellWalkH <- list(infMat=infMat, normMat=normMat, cellLabels=cellLabels) #make cellWalk object
  class(cellWalkH) <- "cellWalk"

  return(list(cellWalkH, weight1))
}

#' Mapping cell type labels
#' \code{mapCellTypes} Random walk on the cells and cell type labels graph and permuted graphs, compute Z-scores for each pair of labels.
#' Cell type labels can be organized in a hierarchical structure.
#' @param cellGraph cell-cell graph, a cell-to-cell similarity matrix, row and colnames are cell barcodes
#' @param labelEdgeList a list of cell-to-label edges. Each is a matrix or dataframe of edges from a cell to each cell type label. Each row is a cell,
#' each column is a cell type label. The cell names of each matrix can be partially overlapped with the cells in \code{cellGraph}.
#' @param labelEdgeWeights a vector of the edge weight ratios between cell-label edges and cell-cell edges. One for each label set.
#' If it is NULL, will tune edge weights and get the optimal value for each  label set;
#' otherwise must have the same length as  \code{labelEdgeList}.
#' @param wtrees a dataframe/matrix the edge weight between nodes on the tree for each set of cell type labels. Each row is a label set. It must have two columns,
#' the first and second column is the edge weights traversing up/down the tree. If NULL, will set all edge weights to 1;
#' otherwise must have the same row as  \code{labelEdgeList}.
#' @param treeList a list of cell type trees, each is a phylo object of the hierarchical structure of a set of cell type labels. Cell type labels in labelEdgeList much contain all the tips on the tree.
#' If NULL, will not include hierarchical relationship between labels in the graph; otherwise must have the same length as \code{labelEdgeList}.
#' @param compute.Zscore whether to do permutation and compute Z-score, If not, will return influence score between labels.
#' @param nround rounds of permutations to compute null distribution. default: 50
#' @param groupsList a list of vectors showing groups of cells for each label set. Only permute edges between cells and labels within the same group of the cells. default: NULL, permute among all cells with all labels.
#' If not NULL, must have the same length as \code{labelEdgeList}. Each vector must have the same length as the rows of each labelEdge matrix. No permutation if group == 0.
#' @param sampleDepth subsample cells to for faster calculation to tune edge weight ratios (\code{labelEdgeWeights}), default: 2000
#' @param ... other arguments pass to \code{tuneEdgeWeights} and \code{walkCells}
#' @return A cellWalk2 object. if \code{compute.Zscore} is False a list including influence score between each set of cell types and labelEdgeWeights.
#' otherwise  a list including influence and Zscore between each set of cell types, as well as labelEdgeWeights.
#' @export

mapCellTypes <- function(cellGraph, labelEdgesList, labelEdgeWeights = NULL, wtrees = NULL, treeList = NULL, compute.Zscore = TRUE,  nround = 50,
                        groupsList = NULL, sampleDepth =2000, ...) # groups, permutation within cell groups, cells same order as cellEdges
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
  diag(cellGraph) = 0

  if(missing(labelEdgesList) || !is(labelEdgesList, "list")){
    stop("Must provide a list of cell-to-label edges matrices")
  }
  if(!is.null(groupsList) & (!is(groupsList, "list") || length(groupsList) != length(labelEdgesList))){
    stop('groupsList must be a list of vectors and must have the same length as labelEdgesList')
  }
  if(is.null(groupsList))
  {
    groupsList = lapply(labelEdgesList, function(labelEdges) rep(1, nrow(labelEdges)))
  }
  if(!is.null(treeList) & length(treeList)!=length(labelEdgesList))
  {
    stop('treeList must have the same length as labelEdgesList')
  }

  for (i in seq_along(labelEdgesList)){
    labelEdges = labelEdgesList[[i]]
    groups = groupsList[[i]]
    if(!is.null(treeList))
    {
      stopifnot('tree must be provided as a phylo object'=is(treeList[[i]], 'phylo'))
      stopifnot('tree must have tip labels' = !is.null(treeList[[i]]$tip.label))
      res = checkLabelEdges(labelEdges, groups, cellGraph, tips = treeList[[i]]$tip.label, suffix = i) #reorder columns of labelEdges as the cell type matrix
    }else{
      res = checkLabelEdges(labelEdges, groups, cellGraph, suffix = i) # adding '_x' to cell type labels to avoid duplicate
    }

    labelEdgesList[[i]] = res[[1]]
    groupsList[[i]] = res[[2]]
  }

  if(!is.null(treeList))
  {
    if(!is.null(wtrees))
    {
      if(!(is(wtrees, "data.frame") | is(wtrees, "matrix")) | nrow(wtrees) != length(labelEdgesList) | ncol(wtrees) !=2 )
      {
         stop('tree edge weights must be a dataframe/matrix with two columns and same number of rows as labelEdgesList')
      }
    }else{
      wtrees = matrix(1, length(labelEdgesList), 2)
    }

    res = lapply(seq_along(treeList), function(i){
      tr = treeList[[i]]
      w = wtrees[i, ]
      tree2Mat(tr,w[1], w[2], i) # adding '_x' to tips to avoid duplicates
    })
    allCellTypes = c(sapply(res, function(x) x[[2]]))
    ncellTypes = sapply(res, function(x) length(x[[2]]))
    cellTypesM = Matrix::bdiag(sapply(res, function(x) x[1])) # becomes a sparse matrix
    ##cellTypesM = as.matrix(cellTypesM)
    colnames(cellTypesM) = rownames(cellTypesM) = allCellTypes

  }else{
    allCellTypes = c(sapply(labelEdgesList, colnames))
    ncellTypes = sapply(labelEdgesList, ncol)
    cellTypesM = matrix(0, length(allCellTypes), length(allCellTypes))
    colnames(cellTypesM) = rownames(cellTypesM) = allCellTypes
  }
  l_all = length(allCellTypes)

  if(is.null(labelEdgeWeights))
  {
    message('tunning labelEdgeWeights... ')
    edgeWeights <- tuneEdgeWeights(cellGraph,
                                   labelEdgesList,
                                   sampleDepth = sampleDepth,
                                   ...)
    opt_idx = which.max(edgeWeights$cellHomogeneity)
    message('cellHomogeneity at optimal edgeWeight:')
    print(edgeWeights[opt_idx,])
    labelEdgeWeights = unlist(edgeWeights[opt_idx, ])
  }else if(!is(labelEdgeWeights, 'numeric') || length(labelEdgeWeights)!=length(labelEdgesList)){
    stop("Must provide a numeric value weight for each set of cell-to-label edges")
  }

  # construct the whole graph permuted or original
  if(!compute.Zscore) nround = 0
  message('run CellWalker:')
  info1_rand = foreach(r = 0:nround) %dopar%
  {
    message('permutation round ', r)
    ## resample labelEdges
    expandLabelEdges = lapply(seq_along(labelEdgesList), function(i)
    {
      labelEdges_rand = labelEdges = labelEdgesList[[i]]
      groups = groupsList[[i]]
      if(r > 0 & any(groups!=0)) # if all groups = 0, no permutation
      {
        for(g in unique(groups)) # don't include group 0 (some of these cells don't connect to any of the labels in current set)
        {
          if(g == 0) next
          cells = which(groups == g)
          labelEdges_rand[cells, ] = simulate_rand0(labelEdges[cells, ])
        }
      }
      if(!is.null(treeList))
      {
        labelEdgeWeights[i] * cbind(labelEdges_rand, matrix(0,dim(labelEdges)[1], treeList[[i]]$Nnode)) #add 0s to labelEdges to allow room for internal nodes
      }else{
        labelEdges_rand * labelEdgeWeights[i]
      }
    })

    cell2label = do.call('cbind', expandLabelEdges)
    combinedGraph = rbind(cbind(cellTypesM,t(cell2label)), #combined graph w/ internal nodes
                          cbind(cell2label,cellGraph))

    # compute CellWalk on expanded graph
    infMat <- randomWalk(combinedGraph, tensorflow = tensorflow, steps=steps) #5 steps is usually enough
    aa = infMat[1:l_all, 1:l_all]
    colnames(aa) = rownames(aa) = allCellTypes
    return(as.matrix(aa))
  }

  info1 = info1_rand[[1]]
  if(compute.Zscore) zscore = compute_zscore(info1, info1_rand[2:(nround+1)], nround)
  idx = cumsum(ncellTypes)
  sts = c(1, idx[1:(length(idx)-1)]+1)
  params = expand.grid(1:length(idx), 1:length(idx))
  params = params[params[,1]!=params[,2,], ]
  info = Map(function(u,v) info1[sts[u]:idx[u],sts[v]:idx[v]], params[,2],  params[,1])
  if(compute.Zscore) zscore = Map(function(u,v) zscore[sts[u]:idx[u],sts[v]:idx[v]], params[,2],  params[,1]) else zscore = NULL
  cellWalk = list(infMat=info, zscore= zscore, labelEdgeWeights = labelEdgeWeights)
  class(cellWalk) = "cellWalk2"
  return(cellWalk)
}

