#' Generate combined graph
#'
#' \code{combineGraph()} combines weighted cell-to-label edges (dimensions
#' c-by-l) with cell-to-cell edges (dimensions c-by-c) into a single large graph
#'
#' @param labelEdges matrix with dimensions c-by-l
#' @param cellEdges matrix with dimensions c-by-c
#' @param weight numeric indicating label edge weight
#' @return combined matrix with dimensions (c+l)-by-(c+l)
#' @export
combineGraph = function(labelEdges, cellEdges, weight=1){
  l = dim(labelEdges)[2]

  combinedGraph = rbind(cbind(matrix(0,l,l),t(weight*labelEdges)),
        cbind(weight*labelEdges,cellEdges))

  combinedGraph
}

#' Generate combined graph with list of label edges
#'
#' \code{combineMultiLabelGraph()} combines list of weighted label edges with
#' cell edges into a single graph
#'
#' @param labelEdgesList list of matrices each with dimensions c-by-l
#' @param cellEdges matrix with dimensions c-by-c
#' @param weights numeric vector indicating label edge weight for each label edge in list
#' @return combined matrix with dimensions (c+L)-by-(c+L) where L is the sum of l's
#' @export
combineMultiLabelGraph = function(labelEdgesList, cellEdges, weights){
  numLabelLists = length(labelEdgesList)
  if(missing(weights)){
    weights = rep(1, numLabelLists)
  }

  labelEdges = weights[1]*labelEdgesList[[1]]
  if(numLabelLists>1){
    for(i in 2:numLabelLists){
      labelEdges = cbind(labelEdges, weights[i]*labelEdgesList[[i]])
    }
  }

  combinedGraph = combineGraph(labelEdges, cellEdges)

  combinedGraph
}

#' Random walk with restarts
#'
#' \code{randomWalk()} solves a random walk with restarts on an adjacency matrix A using
#' the closed form solution for the influence matrix
#' \ifelse{html}{\out{F = r(I - (1 - r)W)<sup>-1</sup>}}{\eqn{F = r(I - (1 - r)W)^{-1}}}
#' where \ifelse{html}{\out{W = D<sup>−1</sup>A}}{\eqn{W = D^{−1}A}} and D is a diagonal
#' matrix of the sums of edge weights for each node and r is the restart probability
#'
#' @param adj adjacency matrix
#' @param r restart probability
#' @return influence matrix, each column is the vector of influences on each row
#' @export
randomWalk = function(adj, r=0.5){
  len = dim(adj)[1]

  D = Matrix::Diagonal(x=Matrix::rowSums(adj))
  W = solve(D)%*%adj

  infMat = r*solve(Matrix::Diagonal(len)-r*W)

  infMat
}

#' Normalize influence matrix
#'
#' \code{normalizeInfluence()} normalizes the label-to-cell portions of the
#' influence matrix for cross-label comparisons
#'
#' @param labelCellInf label-to-cell portion of influence matrix
#' @return normalized influence matrix
#' @export
normalizeInfluence = function(labelCellInf){
  normMat = apply(labelCellInf, 2, function(x) (x-mean(x))/sd(x))
  normMat = apply(normMat, 2, function(x) x/max(x))

  normMat
}

#' Compute cell homogeneity
#'
#' \code{computeCellHomogeneity()} computes the median influence between cells of the same
#' type compared to cells of other types. This serves as a proxy for for how well cells
#' are labeled and can be used to tune parameters such as label edge weight
#'
#' @param cellWalk a cellWalk object
#' @param cellTypes list of names of cell types, if not provided unique labels are used
#' @return cell homogeneity
#' @export
computeCellHomogeneity = function(cellWalk, cellTypes){
  if(missing(cellWalk) || !is(cellWalk, "cellWalk")){
    stop("Must provide a cellWalk object")
  }
  infMat = cellWalk[["infMat"]]
  cellLabels = cellWalk[["cellLabels"]]

  if(missing(cellTypes)){
    cellTypes = unique(cellLabels)
  }

  l = length(cellTypes)
  log(median(unlist(sapply(cellTypes, function(cellType)
    median(sapply(which(cellLabels==cellType), function(x) {
      median(infMat[x+l,which(cellLabels==cellType)+l])/median(infMat[x+l,which(!cellLabels==cellType)+l])
    }))
  ))))
}

#' Compute Jaccard on sprase matrix
#'
#' \code{sparseJaccard()} computes the Jaccard similarity between the rows
#' of a sparse matrix.
#'
#' @param m sparse matrix
#' @return Jaccard similarity matrix between rows of m
#' @export
sparseJaccard = function(m) {
  A = Matrix::tcrossprod(m)
  im = Matrix::which(A>0, arr.ind=TRUE)
  b = Matrix::rowSums(m!=0)

  Aim = A[im]

  J = Matrix::sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  )

  J
}
