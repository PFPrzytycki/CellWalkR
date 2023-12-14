#' Construct cell-to-cell graph from multiple sources
#'
#' \code{constructCellGraph} construct cell-to-cell graph integrating RNASeq, ATACSeq and multiomic data.
#'
#' @param RNA_Mat a gene-by-cell count matrix, RNASeq part of the multiomic data.Rownames are genes and colnames are cells
#' @param ATAC_Mat a peak-by-cell count matrix, ATACSeq part of the multiomic data. Must have barcodes as colnames
#' @param peaks a dataframe of genomic coordinates of peaks of ATACSeq part of the multiomic data. Must contain at least three columns
#'  with names 'seqnames', 'start' and 'end'.
#' @param RNA_Mat2 a gene-by-cell count matrix, unpaired RNASeq data.Rownames are genes and colnames are cells
#' @param ATAC_Mat2 a peak-by-cell count matrix,unpaired ATACSeq data. Must have cell barcodes as colnames
#' @param peaks2 a dataframe of genomic coordinates of peaks of unpaired ATACSeq data. Must contain at least three columns
#'  with names 'seqnames', 'start' and 'end'. Must input peaks2 if \code{ATAC_mat2} is not NULL
#' @param distan distance metrics for compute similarity of ATACSeq profile. Either sparseJaccard, sparseCosine or LSI.
#' @param filterGene filter peaks close to gene body and promoters to compute cell-cell similarity
#' @param ext_region expand gene region by \code{ext_region} bp for filtering peaks. Default: 1e4
#' @param filterFreq filter peaks by proportion of cells having the peaks. a vector of low and high proportion cutoffs
#' @param ndim the PC dimensions to compute cell-cell distance using RNASeq data
#' @param logarithm Whether to take logarithm of the cell-cell distance using ATACSeq data
#' @param ATAC_weight the weight of ATACSeq similarity for combining ATACSeq and RNASeq profile of multiomic data. \code{ATAC_weight} ATAC similarity +
#' (1-\code{ATAC_weight}) RNA similarity. Must to between 0 and 1.
#' @param normalized Whether to normalize the cell-cell distances between 0 and 1
#' @param knn Defines k for the k-nearest neighbor algorithm
#' @param ... additional parameters used by method
#' @return matrix of cell-to-cell similarity
#' @export
constructCellGraph = function(RNA_Mat, ATAC_Mat, peaks, RNA_Mat2=NULL, ATAC_Mat2=NULL, peaks2=NULL, distan = 'sparseCosine', filterGene = T, ext_region = 1e4,
                              ATAC_weight = 0.3, ndim = 1:30, filterFreq = c(0.002, 0.2),  logarithm = T,  normalized = T, knn = 30 )
{

  if(missing(RNA_Mat) | !(is(RNA_Mat, 'Matrix') | is(RNA_Mat, 'matrix') | is(RNA_Mat, 'data.frame')))
  {
     stop('Must input RNASeq part of multiomic data as a matrix or dataframe')
  }
  if(is.null(colnames(RNA_Mat)) | is.null(rownames(RNA_Mat))){
    stop('RNA_Mat must have geneIDs as rownames and cell barcodes  as colnames')
  }
  if(missing(ATAC_Mat) | !(is(ATAC_Mat, 'Matrix') | is(ATAC_Mat, 'matrix') | is(ATAC_Mat, 'data.frame')))
  {
    stop('Must input ATACSeq part of multiomic data as a matrix or dataframe')
  }
  if(is.null(colnames(ATAC_Mat)) | any(colnames(ATAC_Mat) != colnames(RNA_Mat))){
    stop('ATAC_Mat must have cell barcodes as colnames, same order as colnames of RNA_Mat')
  }
  if(missing(peaks) | !is(peaks, 'data.frame'))
  {
    stop('Must input genomic coordinates of peaks for the ATACSeq part of multiomic data as a dataframe')
  }
  if(ncol(peaks) < 3 | is.null(colnames(peaks)) | any(colnames(peaks)[1:3]!=c('seqnames', 'start', 'end'))){
    stop('peaks must have at least 3 columns with names seqnames, start and end')
  }

  if(!is.null(RNA_Mat2) & (!(is(RNA_Mat2, 'Matrix') | is(RNA_Mat2, 'matrix') | is(RNA_Mat2, 'data.frame'))))
  {
    stop('Must input unpaired RNASeq data as a matrix or dataframe')
  }
  if(!is.null(RNA_Mat2) & (is.null(colnames(RNA_Mat2)) | is.null(rownames(RNA_Mat2)))){
    stop('RNA_Mat2 must have geneIDs as rownames and cell barcodes  as colnames')
  }
  if(!is.null(RNA_Mat2) & any(colnames(RNA_Mat2) %in% colnames(RNA_Mat))){
    stop('RNA_Mat2 must have different cell barcodes than RNA_Mat')
  }
  if(!is.null(ATAC_Mat2) & (!(is(ATAC_Mat2, 'Matrix') | is(ATAC_Mat2, 'matrix') | is(ATAC_Mat2, 'data.frame'))))
  {
    stop('Must input unpaired ATACSeq as a matrix or dataframe')
  }
  if(!is.null(ATAC_Mat2) & (is.null(colnames(ATAC_Mat2)) | any(colnames(ATAC_Mat2) %in% colnames(ATAC_Mat)))){
    stop('ATAC_Mat2 must have cell barcodes as colnames, different from colnames of ATAC_Mat')
  }
  if(!is.null(ATAC_Mat2) & (is.null(peaks2) | !is(peaks2, 'data.frame')))
  {
    stop('Must input genomic coordinates of peaks for the unpaired ATACSeq data as a dataframe')
  }

  message('compute cell-cell similarity for paired data')
  cellEdges = constructCellGraphFromATAC(ATAC_Mat, peaks, distan = distan, filterGene = filterGene, ext_region = ext_region,
                                    filterFreq = filterFreq,  logarithm = logarithm,  normalized = normalized)
  cellnames = colnames(cellEdges)
  RNA_results = processRNASeq(RNA_Mat, do.findMarkers = F, computeKNN = F, computeSimilarity = T, buildTree = F, dims = ndim)
  cellEdges2 = RNA_results$cellGraph
  if(logarithm){
    cellEdges2 = -log(1-cellEdges2 + 0.01)
    cellEdges2 = (cellEdges2 - mean(cellEdges2))/sd(cellEdges2)
  }
  stopifnot(all(colnames(cellEdges2) == colnames(cellEdges)))
  cellEdges = cellEdges * ATAC_weight + cellEdges2 * (1 - ATAC_weight)
  if(logarithm)
  {
    cellEdges = (cellEdges - mean(cellEdges))/sd(cellEdges)
  }

  # compute cell-cell similarity for unpaired data
  if(!is.null(RNA_Mat2))
  {
    message('compute cell-cell similarity for unpaired RNASeq data')
    mergeRNA = mergeRNASeq(list(RNA_Mat, RNA_Mat2), integrate = T, computeKNN = F) # generate cell-cell graph combining RNASeq datasets
    cellgraph_R = mergeRNA$cellGraph
  }
  if(!is.null(ATAC_Mat2))
  {
    message('compute cell-cell similarity for unpaired ATACSeq data')
    mergeATAC = mergeATACSeq(ATAC_Mat, peaks, ATAC_Mat2, peaks2)
    cellgraph_A = constructCellGraphFromATAC(mergeATAC$mat, mergeATAC$peaks,  distan = distan, filterGene = filterGene, ext_region = ext_region,
    filterFreq = filterFreq,  logarithm = logarithm,  normalized = normalized)
  }

  message('combine cell graphs')
  if(!is.null(ATAC_Mat2)) {
    cell_Auniq = setdiff(colnames(cellgraph_A), cellnames)
    cellgraph_comb = cellgraph_A[c(cellnames, cell_Auniq), c(cellnames, cell_Auniq)]
    cellgraph_comb[1:length(cellnames), 1:length(cellnames)] = cellEdges
  }else {
    cell_Auniq = NULL
    cellgraph_comb = cellEdges
  }
  if(!is.null(RNA_Mat2)) {
    cell_Runiq = setdiff(colnames(cellgraph_R), cellnames)
    temp = rbind(cellgraph_R[cellnames, cell_Runiq], matrix(0, length(cell_Auniq), length(cell_Runiq)))
    cellgraph_comb = rbind(cbind(cellgraph_comb, temp), cbind(t(temp), cellgraph_R[cell_Runiq, cell_Runiq]))
  }else cell_Runiq = NULL
  cellgraph_comb = as.matrix(cellgraph_comb)

  if(logarithm)
  {
    cellgraph_comb = max(cellgraph_comb) - cellgraph_comb # distance
    diag(cellgraph_comb) = 0
    knn_graph = Seurat::FindNeighbors(cellgraph_comb, distance.matrix = T, k.param = knn)
  }else{
    knn_graph = Seurat::FindNeighbors(1-cellgraph_comb, distance.matrix = T, k.param = knn)
  }
  knn_graph$snn

}

#' Construct cell similarity graph based on ATACSeq data
#'
#' \code{constructCellGraphFromATAC} construct cell-to-cell graph from ATACSeq count data.
#'
#' @param ATAC_Mat a peak-by-cell count matrix. Must have barcodes as colnames
#' @param peaks a dataframe or GRanges object of genomic coordinates of peaks. If it's a dataframe, uust contain at least three columns
#'  with names 'seqnames', 'start' and 'end'
#' @param distan distance metrics for compute similarity of ATACSeq profile. Either sparseJaccard, sparseCosine or LSI.
#' @param filterGene filter peaks close to gene body and promoters to compute cell-cell similarity
#' @param ext_region expand gene region by \code{ext_region} bp for filtering peaks. Default: 1e4
#' @param filterFreq filter peaks by proportion of cells having the peaks. a vector of low and high proportion cutoffs
#' @param logarithm Whether to take logarithm of the cell-cell distance using ATACSeq data
#' @param normalized Whether to normalize the cell-cell distances between 0 and 1
#' @return matrix of cell-to-cell similarity
#' @export
constructCellGraphFromATAC = function(ATAC_Mat, peaks, distan = 'sparseCosine', filterGene = T, ext_region = 1e4,
                                  filterFreq = c(0.002, 0.2),  logarithm = T,  normalized = T)
{
  if(missing(ATAC_Mat) | !(is(ATAC_Mat, 'Matrix') | is(ATAC_Mat, 'matrix') | is(ATAC_Mat, 'data.frame')))
  {
    stop('Must input ATACSeq data as a matrix or dataframe')
  }
  if(is.null(colnames(ATAC_Mat))){
    stop('ATAC_Mat must have cell barcodes as colnames')
  }
  if(missing(peaks))
  {
    stop('Must input genomic coordinates of peaks for the ATACSeq data')
  }else if(is(peaks, 'data.frame')){
    if(ncol(peaks) < 3 | is.null(colnames(peaks)) | any(colnames(peaks)[1:3]!=c('seqnames', 'start', 'end'))){
      stop('peaks must have at least 3 columns with names seqnames, start and end')
    }
    peaks = as(peaks, 'GRanges')
  }else{
    stopifnot('peaks must to a dataframe or granges: '= is(peaks, 'GRanges'))
  }
  seqlevelsStyle(peaks) <- 'UCSC' # add 'chr'
  if(filterGene)
  {
    # filter for peaks near genes
    gene_regions <- getRegions(geneBody = TRUE, genome = "hg38", names = "Entrez") # include gene body and promoter
    gene_regions = GenomicRanges::trim(gene_regions + ext_region) # extend regions both direction
    markerOverlaps = GenomicRanges::findOverlaps(peaks, gene_regions)
    hits = unique(markerOverlaps@from)
    ATAC_Mat = ATAC_Mat[hits, ]
    peaks = peaks[hits]
  }
  if(!is.null(filterFreq))
  {
    # filter by frequency
    idf = Matrix::rowMeans(ATAC_Mat > 0)
    ind = which(idf > filterFreq[1] & idf < filterFreq[2])
    ATAC_Mat = ATAC_Mat[ind,]
    peaks = peaks[ind]
  }
  message('number of peaks to compute cell-cell similarity by ATACSeq:', nrow(ATAC_Mat))

  message('computing cell-cell similarity for ATACSeq data...')
  cellEdges = computeCellSim(Matrix::t(ATAC_Mat),distan) # celll-to-cell similarity
  rownames(cellEdges) = colnames(cellEdges) = colnames(ATAC_Mat)
  cellEdges = as.matrix(cellEdges)
  cellEdges
}

#' Merge two ATACSeq data
#'
#' \code{mergeATACSeq} merge peaks from two ATACSeq data and combine cells to generate a single peak-by-cell matrix
#'
#' @param ATAC_Mat a peak-by-cell count matrix. Must have barcodes as colnames
#' @param peaks a dataframe of genomic coordinates. Must contain at least three columns
#'  with names 'seqnames', 'start' and 'end'.
#' @param ATAC_Mat2 another peak-by-cell count matrix. Must have cell barcodes as colnames different from colnames of \code{ATAC_Mat}
#' @param peaks2 a dataframe of genomic coordinates of peaks the second count matrix. Must contain at least three columns
#'  with names 'seqnames', 'start' and 'end'.
#' @return a peak-by-cell matrix and a GRanges object of peaks
#' @export
mergeATACSeq <-function(ATAC_Mat, peaks, ATAC_Mat2, peaks2)
{
  if(missing(ATAC_Mat) | !(is(ATAC_Mat, 'Matrix') | is(ATAC_Mat, 'matrix') | is(ATAC_Mat, 'data.frame')))
  {
    stop('Must input ATACSeq data as a matrix or dataframe')
  }
  if(is.null(colnames(ATAC_Mat))){
    stop('ATAC_Mat must have cell barcodes as colnames')
  }
  if(missing(peaks) | !is(peaks, 'data.frame'))
  {
    stop('Must input genomic coordinates of peaks for the ATACSeq part of multiomic data as a dataframe')
  }
  if(ncol(peaks) < 3 | is.null(colnames(peaks)) | any(colnames(peaks)[1:3]!=c('seqnames', 'start', 'end'))){
    stop('peaks must have at least 3 columns with names seqnames, start and end')
  }
  if(missing(ATAC_Mat2) & !(is(ATAC_Mat2, 'Matrix') | is(ATAC_Mat2, 'matrix') | is(ATAC_Mat2, 'data.frame')))
  {
    stop('Must input ATACSeq data 2 must be a matrix or dataframe')
  }
  if(is.null(colnames(ATAC_Mat2)) | any(colnames(ATAC_Mat2) %in% colnames(ATAC_Mat))){
    stop('ATAC_Mat2 must have cell barcodes as colnames and distinct from colnames of ATAC_Mat')
  }
  if(missing(peaks2) | !is(peaks2, 'data.frame'))
  {
    stop('Must input genomic coordinates of peaks for the unpaired ATACSeq data as a dataframe')
  }
  if(ncol(peaks2) < 3 | is.null(colnames(peaks2)) | any(colnames(peaks2)[1:3]!=c('seqnames', 'start', 'end'))){
    stop('peaks2 must have at least 3 columns with names seqnames, start and end')
  }
  peaks = as(peaks, 'GRanges')
  peaks2 = as(peaks2, 'GRanges')
  seqlevelsStyle(peaks) <- 'NCBI'
  seqlevelsStyle(peaks2) <- 'NCBI'

  peak.ranges <- Reduce(f = c, x = list(peaks, peaks2))
  peaks_all = GenomicRanges::reduce(x = peak.ranges, ignore.strand=T)
  ATAC_Mat = recomputeATACMat(peaks_all, ATAC_Mat, peaks)
  ATAC_Mat2 = recomputeATACMat(peaks_all, ATAC_Mat2, peaks2)
  ATAC_Mat = cbind(ATAC_Mat, ATAC_Mat2)

  return(list('mat' = ATAC_Mat, 'peaks' = peaks_all))

}

#' Identify motifs in given genomic regoins
#'
#' \code{findMotifs} merge peaks from two ATACSeq data and combine cells to generate a single peak-by-cell matrix
#'
#' @param pRE a dataframe of genomic coordinates. Must contain at least three columns
#'  with names 'seqnames', 'start' and 'end'.
#' @return a dataframe with the first column pREs and second column motifs
#' @export
findMotifs = function(pRE)
{
  # Get a list of motif position frequency matrices from the JASPAR database
  pfm <- TFBSTools::getMatrixSet(
    x = JASPAR2020::JASPAR2020,
    opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
  )

  pRE = GenomicRanges::makeGRangesFromDataFrame(pRE)
  pRE$name = paste0('pRE_', 1:length(pRE))

  chrom_assay <- Signac::CreateChromatinAssay(
    counts = matrix(0, length(pRE), 10, dimnames = list(NULL, 1:10)),
    ranges = pRE,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = NULL,
    min.cells = 0,
    annotation = NULL
  )

  pRE_obj <- Seurat::CreateSeuratObject(counts = chrom_assay, assay = "peaks")
  # add motif information
  pRE_obj <- Signac::AddMotifs(
    object = pRE_obj,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    pfm = pfm
  )
  motifs1 = pRE_obj@assays$peaks@motifs@data
  # output motif matrix
  motifs1_idx = Matrix::summary(motifs1)
  motifs1_idx$i = rownames(motifs1)[motifs1_idx$i]
  motifs1_idx$j = colnames(motifs1)[motifs1_idx$j]
  motifs1_idx$j = pRE_obj@assays$peaks@motifs@motif.names[motifs1_idx$j]
  motifs1_idx$j = gsub("\\(var.[0-9]+\\)",'', motifs1_idx$j)
  colnames(motifs1_idx)[1:2] = c('sequence_name', 'cluster')
  motifs1_idx
}
