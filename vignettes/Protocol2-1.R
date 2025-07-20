library(CellWalkR)
library(ggplot2)
library(data.table)
library(Matrix)
library(ape) # for computation on tree
library(ggtree)
library(dplyr)

library(foreach)
library(doParallel)
cl<-makeCluster(8)  # change to the number of cores you would like to use
registerDoParallel(cl)

ATAC_Mat = SampleCortexSingleCellData$ATAC_Mat
peaks = SampleCortexSingleCellData$peaks
colnames(peaks)[1:3] = c('seqnames', 'start', 'end')
# ATAC-Seq part of multi-omic data: count and peak coordinates
ATAC_Mat0 = SampleCortexSingleCellData$ATAC_Mat0
peaks0 = SampleCortexSingleCellData$peaks0
# RNA-Seq part of multi-omic data
counts = SampleCortexSingleCellData$counts
# scRNA-Seq data
counts2 = SampleCortexSingleCellData$counts2
# cell type tree and markers
RNA_markers = SampleCortexSingleCellData$RNA_markers
tr = SampleCellTypeTree
tr <- root(tr, outgroup = c('Peric.','EC'), resolve.root = TRUE) # reroot the original cell type tree and make non-neuronal cells as outgroup for better visualization

# RNASeq part of multiomic data, normalize data
dataset1 = processRNASeq(counts, do.findMarkers = F, computeKNN = F, computeSimilarity = F, buildTree = F)
# unpaired RNASeq, normalize data
dataset2 = processRNASeq(counts2, do.findMarkers = F, computeKNN = F, computeSimilarity = F, buildTree = F)

#### mapping pREs ####
labelEdges1 = computeTypeEdges(dataset1$expr_norm, RNA_markers)
labelEdges2 = computeTypeEdges(dataset2$expr_norm, RNA_markers)
labelEdges = rbind(labelEdges1, labelEdges2) # combine cells
cellgraph = constructCellGraph(counts, ATAC_Mat0, peaks0, counts2, ATAC_Mat, peaks)

# input genomic coordinates of region-specific pREs
pRE = read.csv(system.file("extdata", "pRE_region_bg_cortex.csv", package = "CellWalkR"))
#rename some column names to be readable by CellWalker2
colnames(pRE)[c(1,11)] = c('seqnames', 'cluster')

# compute cell-to-annotation edge weights from ATAC-Seq
labelEdges1 = computeBulkEdges(pRE, peaks0, ATAC_Mat0)
labelEdges2 = computeBulkEdges(pRE, peaks, ATAC_Mat)
# connect all ATAC cells to bulk annotations
labelEdges2 = rbind(labelEdges1, labelEdges2)

cellWalk2 = annotateBulkRegion(cellgraph, labelEdges, labelEdges2, tr1 = tr, wtree = c(1, 0.1), labelEdgeWeights = NULL,
                               sampleDepth = 2000,  parallel = T, numCores = 8) # with tuning edgeWeights

# add internal node names to the cell type tree
tr$node.label = colnames(cellWalk2$zscore)[-1:-(tr$Nnode+1)]
# plot Z-scores on the cell type tree, only Z-score > 3 is shown
p1 = plotZscoreTree(tr, cellWalk2$zscore, cutoff = 3)
jpeg('/Users/zhu/Gladstone Dropbox/Zhirui Hu/cell_hierarchy/CellGenomics_Submission/STARProtocols/Figures/Figure4.jpg',
     width = 16, height = 13, res = 300, quality = 100, units = 'cm')
p1
dev.off()

#### mapping TFs ####
pRE = read.table(system.file("extdata", "pRE-hg38.bed", package = "CellWalkR"))
colnames(pRE) = c('seqnames', 'start', 'end')
motifs = findMotifs(pRE)
# select motifs appear in more than 1000 pREs
motifs = data.table::as.data.table(motifs)
motifs[, count:= .N, by = cluster]
motifs = as.data.frame(motifs[count > 1000])
regionMat = convertToMatrix(motifs) # a data.table with sequence name of pRE as the first column and TFs as the following columns
labelEdges1 = computeBulkEdges(motifs, peaks0, ATAC_Mat0)
labelEdges2 = computeBulkEdges(motifs, peaks, ATAC_Mat)
labelEdges2 = rbind(labelEdges1, labelEdges2)
# no permutation between cell-to cell type labels edges
groups1 =  rep(0, nrow(labelEdges))
# permutation between all motifs and regions
groups2 = rep(1, nrow(regionMat))
cellWalk2 = annotateBulkRegion(cellgraph, labelEdges, labelEdges2, groups1, groups2, regionMat,
                               list(ATAC_Mat0, ATAC_Mat), list(peaks0, peaks), tr1 = tr, wtree = c(1, 0.1),
                               labelEdgeWeights = NULL, sampleDepth = 2000, parallel = T, numCores = 8)
# reorder cell types so that cell types closer on the hierarchy will appear in adjacent rows on the heatmap
Zscore = cellWalk2$zscore[, tr$edge[,2]]
# Z-score > 5 is shown
p1 = plotZscoreDotplot(Zscore, th = 5)
jpeg('/Users/zhu/Gladstone Dropbox/Zhirui Hu/cell_hierarchy/CellGenomics_Submission/STARProtocols/Figures/Figure5.jpg',
     width = 16, height = 17, res = 300, quality = 100, units = 'cm')
p1
dev.off()

# Load log normalized count for each TF in each cell
tf_exp = SampleCortexSingleCellData$tf_exp
label = SampleCortexSingleCellData$label # cell labels
# remove cell types less than 5 cells
ct_num = xtabs(~label)
clusters = names(ct_num)[which(ct_num > 4)]
cell_ind = which(label %in% clusters)
label = label[cell_ind]
tf_exp = tf_exp[, cell_ind]
res = computeTFexp(tf_exp, Zscore, label)
# compute correlation between Z-score and expression
corr = sapply(1:nrow(res[[1]]), function(x) cor(unlist(res[[1]][x,]), unlist(res[[2]][x,]), method= 'spearman'))
names(corr) = rownames(res[[1]])
posTF = corr[corr > 0.2] # less strict cutoff for demo
# For each cell type, select TFs with larger Z-score and expression among all the positively correlated TFs
nn = sapply(1:ncol(res[[1]]), function(i) { # iterate over all the cell types
  y = sum(res[[1]][,i] >qnorm(0.995)) # select TFs with Z-score > 3
  x = order(-res[[1]][,i])[1:min(50,y)] # further select top 50 TFs if too many to plot
  x = rownames(res[[1]])[x]
  intersect(intersect(x, rownames(res[[2]])[res[[2]][,i] > 0.5]), names(posTF)) # select positive correlated TF with standardized expression (>0.5)
})
names(nn) = colnames(res[[1]])
tfs = unique(unlist(nn))
ord2 = na.omit(tr$tip.label[tr$edge[,2]]) # reorder cell types shown on the heatmap so that cell types closer on the hierarchy will appear in adjacent rows on the heatmap
# convert the format of Z-scores and gene expression for ggplot
res2 = convert2plot(res[[1]][tfs,], res[[2]][tfs,], ord2 = ord2, th = qnorm(0.995))
jpeg('/Users/zhu/Gladstone Dropbox/Zhirui Hu/cell_hierarchy/CellGenomics_Submission/STARProtocols/Figures/Figure6.jpg',
     width = 12, height = 15, res = 300, quality = 100, units = 'cm')
ggplot(res2[[1]], aes(x=enhancer , y=celltype, group=enhancer)) +
  geom_tile(data = res2[[2]], mapping =  aes(fill= expression)) + geom_point(aes(size = zscore),  alpha = 0.8) + xlab('') +
  theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  scale_fill_gradient2(low = "mediumblue",  high = "red2", space = "Lab") + scale_size(range = c(0.5, 6))
dev.off()

# compute standardized gene expression per cell type including internal nodes
res = computeTFexp(tf_exp, Zscore, label, tr = tr, levels=c(1:8))

# For each cell type, select TFs with larger Z-score and expression among all the positively correlated TFs
nn = sapply(1:ncol(res[[1]]), function(i) { # iterate over all the cell types
  y = sum(res[[1]][,i] >qnorm(0.995)) # select TFs with Z-score > 3
  x = order(-res[[1]][,i])[1:min(50,y)] # further top 50 TFs if too many to plot
  x = rownames(res[[1]])[x]
  intersect(intersect(x, rownames(res[[2]])[res[[2]][,i] > 0.5]), names(posTF)) # select TF with standardized expression (>0.5) and correlation
})
names(nn) = colnames(res[[1]])

# rename the nodes of the cell type tree to be the select TFs for plotting
tr1 = tr
tr1$node.label = sapply(nn[tr$node.label], function(x) paste(x, collapse = ','))
tr1$tip.label = sapply(nn[tr$tip.label], function(x) paste(x, collapse = ','))
# reduce TF labels on the cell type tree: each TF will only appear in the most upstream node in which it is selected
nodes = which(!is.na(names(tr1$node.label)))
for(e in nodes) # traverse the nodes from bottom to top of the tree
{
  st = subtrees(tr1)
  a = tr1$node.label[e]
  a = strsplit(a, ',')[[1]]
  for(i in 1:length(st[[e]]$tip.label))
  {
    b = st[[e]]$tip.label[i]
    if(b == '') next
    tr1$tip.label[names(b)] = paste(setdiff(strsplit(b, ',')[[1]], a), collapse = ',')
  }
  if(length(st[[e]]$node.label) < 2) next
  for(i in 2:length(st[[e]]$node.label))
  {
    b = st[[e]]$node.label[i]
    if(b == '') next
    tr1$node.label[names(b)] = paste(setdiff(strsplit(b, ',')[[1]], a), collapse = ',')
  }
}
td <- data.frame(node = 1:length(tr$tip.label),
                 'TFs' = tr1$tip.label, check.names = F)
nd <- data.frame(node = (1:length(tr$node.label)) + length(tr$tip.label),
                 'TFs' = tr1$node.label, check.names = F)
d <- rbind(td, nd)
d[d$TFs=='', 'TFs'] = NA
tree <- full_join(tr, d, by = 'node')

jpeg('/Users/zhu/Gladstone Dropbox/Zhirui Hu/cell_hierarchy/CellGenomics_Submission/STARProtocols/Figures/Figure7.jpg',
     width = 17.2, height = 15, res = 300, quality = 100, units = 'cm')
ggtree(tree, branch.length = 'none') + geom_tiplab(hjust = -0.3) + xlim(0, 16) + geom_label(aes(x = branch, label=TFs)) +
  theme(text = element_text(size = 9))
dev.off()

#### input ChIP-Seq data ####
chipseq = readRDS(system.file("extdata", 'remap2022_nr_macs2_hg38_brain_tfs.rds', package = "CellWalkR"))
# overlap ChIPSeq peaks with pRE
chip_peaks = GRanges(chipseq[,1:3])
pRE = GRanges(pRE)
overlaps = findOverlaps(chip_peaks, pRE)
chipseq2 = cbind(chipseq[queryHits(overlaps), ], as.data.frame(pRE[subjectHits(overlaps),]))
chipseq2 = unique(chipseq2[, c(5:7,4)]) # select columns: seqnames, start, end, cluster.
regionMat = convertToMatrix(chipseq2) # a data.table with sequence name of pRE as the first column and clusters (TFs) as the following columns

labelEdges1 = computeBulkEdges(chipseq2, peaks0, ATAC_Mat0)
labelEdges2 = computeBulkEdges(chipseq2, peaks, ATAC_Mat)
labelEdges2 = rbind(labelEdges1, labelEdges2)

groups1 =  rep(0, nrow(labelEdges)) # no permutation between cell-to cell type labels edges
groups2 = rep(1, nrow(regionMat)) # permutation between all TFs and regions
cellWalk2_chip = annotateBulkRegion(cellgraph, labelEdges, labelEdges2, groups1, groups2, regionMat, list(ATAC_Mat0, ATAC_Mat), list(peaks0, peaks), tr1 = tr, wtree = c(1, 0.1), labelEdgeWeights = NULL, sampleDepth = 2000, parallel = T, numCores = 8)

Zscore_chip = cellWalk2_chip$zscore[, tr$edge[,2]] # reorder cell types shown on the heatmap
p1 = plotZscoreDotplot(Zscore_chip, th = 2)
jpeg('/Users/zhu/Gladstone Dropbox/Zhirui Hu/cell_hierarchy/CellGenomics_Submission/STARProtocols/Figures/Figure8.jpg',
     width = 10, height = 17, res = 300, quality = 100, units = 'cm')
p1
dev.off()







