# check label Edges and make labelEdges having the same rows as cellGraph
checkLabelEdges <- function(labelEdges, groups, cellGraph, tips = NULL, suffix = NULL) {
  if (missing(labelEdges) ||
      !inherits(labelEdges, c("data.frame", "matrix", "Matrix"))) {
    stop("Must provide a dataframe or matrix of cell-to-label edges")
  }
  if (is.null(colnames(labelEdges)) || is.null(rownames(labelEdges))) {
    stop("labelEdges must have cell barcodes as rownames and labels as colnames")
  }
  if (!is.numeric(groups) || length(groups) != nrow(labelEdges)) {
    stop("groups must be numeric with length = nrow(labelEdges)")
  }

  names(groups) <- rownames(labelEdges)
  stopifnot("groups contain NA" = all(!is.na(groups)))

  # Early return if already aligned
  if (nrow(labelEdges) == nrow(cellGraph) &&
      all(rownames(labelEdges) == rownames(cellGraph))) {
    return(list(labelEdges, groups))
  }

  # Align to cellGraph
  rowIdx <- match(rownames(cellGraph), rownames(labelEdges))
  if (all(is.na(rowIdx))) {
    warning('No common cell barcodes between labelEdges and cellGraph')
  }

  groups <- groups[rowIdx]
  labelEdges <- labelEdges[rowIdx, , drop = FALSE]
  groups[is.na(groups)] <- 0
  labelEdges[is.na(labelEdges)] <- 0
  rownames(labelEdges) <- names(groups) <- rownames(cellGraph)

  # Filter to tree tips if provided
  if (!is.null(tips)) {
    if (!all(tips %in% colnames(labelEdges))) {
      stop('Missing cell-to-label edges for tree tips: ',
           paste(setdiff(tips, colnames(labelEdges)), collapse = ", "))
    }
    labelEdges <- labelEdges[, tips, drop = FALSE]
  }

  # Add suffix if provided
  if (!is.null(suffix)) {
    colnames(labelEdges) <- paste(colnames(labelEdges), suffix, sep = '_')
  }

  return(list(labelEdges, groups))
}


# convert ape tree to adjacency matrix and add internal node names
# adjacacency matrix can be asymmetric as assigning different weights going up/down the tree

#tree2Mat <- function(tr1, w_up = 1, w_down = 1, suffix = NULL)
#{
#  cellTypesM = matrix(0, tr1$Nnode * 2 + 1, tr1$Nnode * 2 + 1)
#  ntips = (tr1$Nnode+1)
#  nn = tr1$Nnode * 2 + 1
#  for(i in 1:nrow(tr1$edge)){
#    xx = tr1$edge[i, ]
#    k = xx[1]; l = xx[2]
#    if(k > ntips) k = nn  - k + ntips + 1
#    if(l > ntips) l = nn  - l + ntips + 1
#    cellTypesM[k, l] = w_down
#    cellTypesM[l, k] = w_up
#  }
#  allCellTypes = c(tr1$tip.label, rep(0, tr1$Nnode))
#  for(j in (ntips+1):nn)  # rename internal nodes
#  {
#    idx = which(cellTypesM[j, 1:j]!=0)
#    stopifnot(length(idx) == 2)
#    names = allCellTypes[idx]
#    if(any(is.na(names))) stop("name internal node error")
#    ly1 = ly2 = 0
#    if(grepl(':', names[1])) {
#      ly1 = as.numeric(strsplit(names[1], ':')[[1]][3])
#      names[1] = strsplit(names[1], ':')[[1]][1]
#    }
#    if(grepl(':', names[2])) {
#      ly2 = as.numeric(strsplit(names[2], ':')[[1]][3])
#      names[2] = strsplit(names[2], ':')[[1]][1]
#    }
#    allCellTypes[j] = paste(names[1], names[2], max(ly1,ly2)+1, sep=':')
#  }
#  if(!is.null(suffix))
#  {
#    allCellTypes = paste(allCellTypes, suffix, sep='_')
#  }
#  colnames(cellTypesM) = rownames(cellTypesM) = allCellTypes
#  return(list(cellTypesM, allCellTypes))
#}

tree2Mat = function(tr,  w_up = 1, w_down = 1, suffix = NULL) {
  CellTypes = tr$tip.label
  ll_all = 2*length(CellTypes) - 1
  A = matrix(0, ll_all, ll_all)
  allCellTypes = c(CellTypes, rep(NA, (length(CellTypes)- 1)))
  edge_sort = tr$edge[order(-tr$edge[,1]), ]
  for(i in 1:nrow(edge_sort)){
    children =  edge_sort[i,]
    A[children[1], children[2]] = w_down
    A[children[2], children[1]] = w_up
    if(is.na(allCellTypes[children[1]]) & sum(A[, children[1]]!=0)==2)
    {
      names = allCellTypes[which(A[, children[1]]!=0)]
      if(any(is.na(names))) stop("name internal node error")
      ly1 = ly2 = 0
      if(grepl(':', names[1])) {
        ly1 = as.numeric(strsplit(names[1], ':')[[1]][3])
        names[1] = strsplit(names[1], ':')[[1]][1]
      }
      if(grepl(':', names[2])) {
        ly2 = as.numeric(strsplit(names[2], ':')[[1]][3])
        names[2] = strsplit(names[2], ':')[[1]][2]
      }
      allCellTypes[children[1]] = paste(names[1], names[2], max(ly1,ly2)+1, sep=':')
    }
  }
  if(!is.null(suffix))
  {
    allCellTypes = paste(allCellTypes, suffix, sep='_')
  }
  colnames(A) = rownames(A) = allCellTypes
  list(A, allCellTypes)
}


# Optimized random permutation with proper edge case handling
simulate_rand0 <- function(labelEdges2) {
  # Quick returns for edge cases
  if (nrow(labelEdges2) <= 1 || ncol(labelEdges2) == 0) {
    return(labelEdges2)
  }

  # Vectorized margin computations
  labelEdges2_capped <- pmin(labelEdges2, 1)
  cuts <- seq(0, 1, by = 0.1)
  n_bins <- length(cuts)

  # Cell margins: probability distribution across bins for each cell
  cell_margin <- apply(labelEdges2_capped, 1, function(x) {
    cumsum_bins <- sapply(cuts, function(y) mean(x <= y))
    bin_probs <- diff(c(0, cumsum_bins))
    rep(bin_probs, each = 2) / 2
  })

  # Label margins: count of values in each bin for each label
  label_margin <- apply(labelEdges2_capped, 2, function(x) {
    cumsum_bins <- sapply(cuts, function(y) sum(x <= y))
    diff(c(0, cumsum_bins))
  })

  # Initialize output
  labelEdges2_rand <- matrix(0, nrow(labelEdges2), ncol(labelEdges2))

  # Assign values by label
  for (j in seq_len(ncol(labelEdges2))) {
    probs <- label_margin[, j]
    idx <- seq_len(nrow(labelEdges2))

    # Iterate through bins from highest to lowest
    for (k in length(probs):2) {
      n_needed <- probs[k]
      if (n_needed == 0 || length(idx) == 0) next

      # Sample cells for this bin
      if (length(idx) <= n_needed) {
        sampled_idx <- idx
      } else {
        probs_vec <- pmax(cell_margin[k, idx], 1e-10)
        sampled_idx <- sample(idx, n_needed, prob = probs_vec)
      }

      # Assign random values within bin range
      labelEdges2_rand[sampled_idx, j] <- runif(
        length(sampled_idx),
        cuts[k - 1],
        cuts[k]
      )

      idx <- setdiff(idx, sampled_idx)
    }
  }

  # Preserve names
  dimnames(labelEdges2_rand) <- dimnames(labelEdges2)
  return(labelEdges2_rand)
}


simulate_rand_binary0 <- function(labelEdges2)
{
  cell_margin = rowMeans(labelEdges2)
  label_margin = colSums(labelEdges2)

  labelEdges2_rand = matrix(0, nrow(labelEdges2), ncol(labelEdges2))

  for(j in 1:length(label_margin))
  {
    prob = label_margin[j]
    idx  = sample(length(cell_margin), prob, prob = cell_margin)
    labelEdges2_rand[idx, j] = 1
  }
  colnames(labelEdges2_rand) = colnames(labelEdges2)
  rownames(labelEdges2_rand) = rownames(labelEdges2)
  return(labelEdges2_rand)
}

simulate_rand <- function(labelEdges2) #quantile
{
  cuts = quantile(c(labelEdges2[labelEdges2>0]), seq(0,1,by = 0.1))
  cell_margin = apply(labelEdges2, 1, function(x){
    l = sapply(cuts[seq(1,11,by = 2)], function(y) mean(x <= y))
    c(l[1], rep(diff(l), each = 2)/2)
  })
  label_margin = apply(labelEdges2, 2, function(x){
    l = sapply(cuts, function(y) sum(x <= y))
    c(l[1], diff(l))
  })
  labelEdges2_rand = matrix(0, nrow(labelEdges2), ncol(labelEdges2))

  cuts = quantile(c(labelEdges2[labelEdges2>0]), c(seq(0,0.9,by = 0.1), 0.98))
  for(j in 1:ncol(label_margin))
  {
    probs = label_margin[, j]
    idx = 1:ncol(cell_margin)
    for(k in length(probs):2)
    {
      if(length(idx) == probs[k])
      {
        aa = idx
      }else{
        aa  = sample(idx, probs[k], prob = cell_margin[k, idx] + 0.01)
      }
      if(k > 1)
      {
        labelEdges2_rand[aa, j] = runif(length(aa), cuts[k-1], cuts[k])
      }
      idx = setdiff(idx, aa)
    }
  }
  colnames(labelEdges2_rand) = colnames(labelEdges2)
  rownames(labelEdges2_rand) = rownames(labelEdges2)
  return(labelEdges2_rand)
}

# Compute Z-score with vectorized operations
compute_zscore <- function(info, info_rand, nround) {
  # Input validation
  if (nround < 1) {
    stop("nround must be at least 1")
  }

  dims <- dim(info)
  if (any(sapply(info_rand, function(x) !identical(dim(x), dims)))) {
    stop("All matrices in info_rand must have same dimensions as info")
  }

  # Vectorized computation using array operations
  info_array <- array(unlist(info_rand), dim = c(dims, nround))
  info_mean <- rowMeans(info_array, dims = 2)
  info_sq_mean <- rowMeans(info_array^2, dims = 2)

  # Variance: E[X²] - E[X]²
  info_var <- info_sq_mean - info_mean^2

  # Handle numerical issues
  if (any(info_var < 0, na.rm = TRUE)) {
    message('Some variances are negative (numerical precision issue); setting to zero')
    info_var[info_var < 0] <- 0
  }

  # Z-score calculation
  info_std <- sqrt(info_var)
  info_std[info_std == 0] <- 1  # Avoid division by zero
  zscore <- (info - info_mean) / info_std
  zscore[zscore < 0 | is.na(zscore)] <- 0

  return(list('zscore' = zscore, 'mean' = info_mean, 'var' = info_var))
}


# Compute Cosine distance on sparse matrix between the rows
sparseCosine = function(m) {
  stopifnot(is(m, 'matrix') | is(m, 'Matrix'))
  if(!is(m, 'Matrix')) m = as(m, 'Matrix')
  if(!is(m, 'dgCMatrix'))
  {
    m = as(as(m, "generalMatrix"), "CsparseMatrix")
  }
  A = Matrix::tcrossprod(m)
  im = Matrix::summary(A) #which(A>0, arr.ind=TRUE)
  b = Matrix::rowSums(m!=0)

  Aim = im[,3] #A[im]

  J = Matrix::sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / sqrt(b[im[,1]] * b[im[,2]]),
    dims = dim(A),
    symmetric = T
  )

  J
}

# Compute Jaccard distance on sparse matrix between the rows
sparseJaccard = function(m) {
  stopifnot(is(m, 'matrix') | is(m, 'Matrix'))
  if(!is(m, 'Matrix')) m = as(m, 'Matrix')
  if(!is(m, 'dgCMatrix'))
  {
    m = as(as(m, "generalMatrix"), "CsparseMatrix")
  }
  A = Matrix::tcrossprod(m)
  im = Matrix::summary(A)
  b = Matrix::rowSums(m!=0)

  Aim = im[,3]

  J = Matrix::sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A),
    symmetric = T
  )

  J
}

# Compute Euclidean distance on latent semantics index
LSI = function(m, ndim = 30) {
  atac <- Seurat::CreateSeuratObject(counts = t(m))
  atac <- Signac::FindTopFeatures(atac, min.cutoff = 10)
  atac <- Signac::RunTFIDF(atac)
  atac <- Signac::RunSVD(atac)
  distance <- dist(atac@reductions$lsi@cell.embeddings[,1:ndim]) #include the weight of each PC
  distance = distance/max(distance)
  cellEdges = as.matrix(1 - distance)
  diag(cellEdges) = 1
  return(cellEdges)
}

# recompute ATAC counts to a new peakset
recomputeATACMat = function(peaks_all, ATAC_Mat, peaks)
{
  markerOverlaps = GenomicRanges::findOverlaps(peaks, peaks_all)
  ATAC_Mat = as.data.table(ATAC_Mat)
  ATAC_Mat = ATAC_Mat[markerOverlaps@from,]
  ATAC_Mat[, subH := markerOverlaps@to]
  ATAC_Mat = ATAC_Mat[, lapply(.SD, sum), by = subH]
  # fill in zero
  add_matrix = data.table('subH' = 1:length(peaks_all))
  add_matrix = merge(add_matrix, ATAC_Mat, all.x = T)
  add_matrix[,subH := NULL]
  add_matrix[is.na(add_matrix)] = 0
  as(add_matrix, 'Matrix')
}

computeEdgeLabels = function(
    labelGenes,
    ATACMat,
    ATACGenePeak,
    cellPeakCounts,
    method="None",
    drop = T)
{
  labels = unique(labelGenes[,2])
  cellsInMarkers = c()
  for(label in labels){
    if(!label %in% dimnames(ATACGenePeak)[[2]]){
      stop("A label is missing from the peak to region mapping. Is it the right mapping?")
    }
    peaksInMarkers = ATACGenePeak[["peak",label]]
    idx = !duplicated(peaksInMarkers) # if a peak overlap with multiple regions, the weight used will be one of them
    peaksInMarkers = peaksInMarkers[idx]
    genesInMarkers = ATACGenePeak[["gene",label]][idx]
    if(max(peaksInMarkers)>dim(ATACMat)[2]){
      stop("Indecies dont match in peak to region mapping. Is it the right mapping?")
    }
    if(length(which(genesInMarkers %in% labelGenes[,1]))!=length(genesInMarkers)){
      stop("Genes dont match in peak to region mapping. Is it the right mapping?")
    }
    if(method=="None"){
      geneExp = rep(1, length(genesInMarkers))
    }else if(method=="Weight"){
      if(dim(labelGenes)[2]<3){
        stop("Must provide a table with genes of interest in first column and corresponding labels in
                 second column and log fold gene expression in the third")
      }
      geneExp = labelGenes[labelGenes[,2]==label,3][match(genesInMarkers,labelGenes[labelGenes[,2]==label,1])]
    }
    else{
      stop("Not a recognized scaling method")
    }
    cellsInLabelMarkers = Matrix::colSums(Matrix::t(ATACMat[,peaksInMarkers])*geneExp)/cellPeakCounts
    cellsInLabelMarkers[cellsInLabelMarkers<0] = 0
    cellsInMarkers = cbind(cellsInMarkers,cellsInLabelMarkers)
  }
  colnames(cellsInMarkers) = labels
  if(drop) cellsInMarkers = cellsInMarkers[,colSums(cellsInMarkers)!=0]
  cellsInMarkers
}

computeLabelHomogeneity = function(cellWalk, cellTypes, bulkLabels){
  if(missing(cellWalk) || !is(cellWalk, "cellWalk")){
    stop("Must provide a cellWalk object")
  }
  infMat = cellWalk[["infMat"]]
  cellLabels = cellWalk[["cellLabels"]]

  if(missing(cellTypes) | missing(bulkLabels)){
    stop('MUst provide cellTypes and bulkLabels')
  }

  l = length(cellTypes)
  infMat = infMat[(l+1):(l+length(bulkLabels)), 1:l]
  y = sapply(1:length(bulkLabels), function(i){
    x = infMat[i, ]
    x = x/sum(x)
    x = x[x > 0]
    sum(x * log(x)) / log(l) + 1
  })
  median(y)
}

convertSeqnames = function(regions, sep = c(':', '-'))
{
  sep = paste(sep, collapse = '|')
  cols = setdiff(colnames(regions), 'sequence_name')
  loci = t(sapply(regions$sequence_name, function(x) strsplit(x, sep)[[1]]))
  colnames(loci) = c('seqnames', 'start', 'end')
  regions <- cbind(loci, regions[,cols])
  return(regions)
}


# Helper function for convergence checking
check_convergence <- function(info_rand_list, window = 20, threshold = 0.01) {
  n <- length(info_rand_list)
  if (n < 2 * window) return(FALSE)

  # Compare recent vs. previous window
  recent_idx <- (n - window + 1):n
  prev_idx <- (n - 2 * window + 1):(n - window)

  recent_mean <- Reduce('+', info_rand_list[recent_idx]) / window
  prev_mean <- Reduce('+', info_rand_list[prev_idx]) / window

  # Mean absolute difference
  diff <- mean(abs(recent_mean - prev_mean))

  message(sprintf(
    "  Convergence check (round %d): MAD = %.6f (threshold = %.6f)",
    n, diff, threshold
  ))

  return(diff < threshold)
}
