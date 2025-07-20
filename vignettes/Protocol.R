# Protocol
library(CellWalkR)
library(ggplot2)
library(data.table)
library(Matrix)
library(ape) # for computation on tree

library(foreach)
library(doParallel)
cl<-makeCluster(8)  # change to the number of cores you would like to use
registerDoParallel(cl)

# count matrices
counts1 = SampleSingleCellRNASeq$counts1 # reference data
counts2 = SampleSingleCellRNASeq$counts2 # query data
# metadata matrices
meta.data1 = SampleSingleCellRNASeq$meta.data1 # reference
meta.data2 = SampleSingleCellRNASeq$meta.data2 # query

dataset1= processRNASeq(counts1, meta.data1, group.col = 'Original_annotation', do.findMarkers = F, computeKNN = F, buildTree = F)

markers = SampleSingleCellRNASeq$RNA_markers1
tree = SampleSingleCellRNASeq$tree1
# only keep the cell types (tips of the tree) in the dataset
tree = ape::keep.tip(tree, unique(meta.data1$Original_annotation))
dataset1$tr = tree
dataset1$markers = markers

#### labeling cells ####
dataset2 = processRNASeq(counts2, do.findMarkers = F, computeKNN = T)
labelEdges = computeTypeEdges(dataset2$expr_norm, dataset1$markers)
results = annotateCells(dataset2$cellGraph, labelEdges, weight1 = NULL, sampleDepth = 3000, labelEdgeOpts = 10^seq(-5,3,1), tr1 = dataset1$tr, wtree = 10)
cellLabel = results[[1]]$cellLabels
cellLabel = cellLabel[rownames(meta.data2)] # make sure the cell barcodes are matched between the results and the meta data

aa = data.table('Yoshida' = meta.data2$Original_annotation, 'Ren' = cellLabel) # combine with original annotation
aa = aa[, list(count= .N), by = c('Yoshida', 'Ren')]
aa[, prob:= count/sum(count), by =  'Yoshida']
# reorder the cell types for better visualization
aa$Ren = factor(aa$Ren, levels = SampleSingleCellRNASeq$label_ord1)
aa$Yoshida = factor(aa$Yoshida, levels = SampleSingleCellRNASeq$label_ord2)

jpeg('/Users/zhu/Gladstone Dropbox/Zhirui Hu/cell_hierarchy/CellGenomics_Submission/STARProtocols/Figures/Figure1.jpg',
     width = 17.2, height = 13, res = 300, quality = 100, units = 'cm')
ggplot(aa, aes(x= Ren, y=Yoshida, size=count, color=prob, group=Ren)) +
  geom_point(alpha = 0.8) +
  theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1), text = element_text(size = 10), legend.position = "top") +
  scale_color_gradient(low = "mediumblue",  high = "red2", space = "Lab", limit = c(0, 1))+scale_size(range = c(0.5, 6))
dev.off()

#### comparing cell types ###
dataset2$markers = SampleSingleCellRNASeq$RNA_markers2
tree = SampleSingleCellRNASeq$tree2
# only keep the cell types (tips of the tree) in the dataset
tree = ape::keep.tip(tree, unique(meta.data2$Original_annotation))
dataset2$tr = tree

labelEdges1 = computeTypeEdges(dataset1$expr_norm, dataset1$markers)
labelEdges2 = computeTypeEdges(dataset2$expr_norm, dataset2$markers, log2FC.cutoff = 0.25) # adjust log2FC to be smaller to increase number of markers
labelEdgesList = list(labelEdges1, labelEdges2)
mergeResult =  mergeRNASeq(list(counts1, counts2), nfeatures = 5000)
treeList = list(dataset1$tr, dataset2$tr)
# different weights going up/down the tree
wtrees = matrix(c(1,0.1,0.1,1), ncol=2) # first row for dataset1 and second for dataset2
groupsList = list(rep(0, nrow(labelEdges1)), rep(1, nrow(labelEdges2)))
cellWalk2_map = mapCellTypes(mergeResult$cellGraph, labelEdgesList, labelEdgeWeights = NULL, treeList = treeList, wtrees = wtrees,
                             groupsList = groupsList, compute.Zscore = TRUE, nround = 55, parallel = T, sampleDepth = 1000, numCores = 8)
Zscore = cellWalk2_map$zscore[[2]]

bb = reshape2::melt(Zscore[1:((nrow(Zscore) + 1)/2),1:((ncol(Zscore) + 1)/2)]) # only plot tip cell types for illustration
# removing suffix to cell type labels adding in mapCellTypes
bb$Var1 = sub( '_[0-9]$', '',bb$Var1)
bb$Var2 = sub( '_[0-9]$', '',bb$Var2)
colnames(bb) = c('Yoshida','Ren', 'Zscores')
# reorder cell types for better visualization
bb$Ren = factor(bb$Ren, levels = SampleSingleCellRNASeq$label_ord1)
bb$Yoshida = factor(bb$Yoshida, levels = SampleSingleCellRNASeq$label_ord2)

jpeg('/Users/zhu/Gladstone Dropbox/Zhirui Hu/cell_hierarchy/CellGenomics_Submission/STARProtocols/Figures/Figure2.jpg',
     width = 17.2, height = 13, res = 300, quality = 100, units = 'cm')
# generate heatmap for Z-scores
ggplot(bb, aes(x= Ren, y=Yoshida, group=Ren)) +
  geom_tile(aes(fill=Zscores), alpha = 0.6) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 10),legend.position = "top") +
  scale_fill_gradient(low = "white",  high = "red2", space = "Lab")
dev.off()

tr = treeList[[1]] # reference cell type tree
# add internal node labels to the tree
tr$node.label = colnames(Zscore)[-1:-(tr$Nnode+1)]
tr$node.label = sub('_[0-9]$', '',tr$node.label)
# select two example cell types to plot
dat = Zscore[c('NK_2', 'T CD8 EMRA_2'),, drop=F]
colnames(dat) =  sub('_[0-9]$', '',colnames(dat))
rownames(dat) =  sub('_[0-9]$', '',rownames(dat))
# extract subtree with positive Z-scores of example cell types for better visualization
cl1 = ape::extract.clade(tr, 'T_CD8_c07-TYROBP:T_gdT_c14-TRDV2:3')
cl2 = ape::extract.clade(tr, 'T_CD8_c04-COTL1:T_CD4_c09-GZMK-FOS_l:3')
subtr = ape::keep.tip(tr, c(cl1$tip.label, cl2$tip.label))
subtr = ape::keep.tip(tr, c(cl1$tip.label, cl2$tip.label))
# plot Z-score on the subtree, only Z-score > 15 is shown
jpeg('/Users/zhu/Gladstone Dropbox/Zhirui Hu/cell_hierarchy/CellGenomics_Submission/STARProtocols/Figures/Figure3.jpg',
     width = 16, height = 11, res = 300, quality = 100, units = 'cm')
pp = CellWalkR::plotZscoreTree(subtr, dat[,c(subtr$tip.label, subtr$node.label)], cutoff = 15)
plot(pp)
dev.off()





