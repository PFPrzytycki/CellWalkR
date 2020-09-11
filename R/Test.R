# $ R
# > library(devtools)
# > install_github("PFPrzytycki/CellWalkR")

# #Testing
#
# #load sample data
# ATACMat = Matrix::readMM("../CellWalkRSampleData/SamplePeakMat.mtx")
# peaks = as(data.table::fread("../CellWalkRSampleData/SamplePeaks.txt", header = FALSE)$V1, "GRanges")
# #load data in snap form
# # snapData = readRDS("../../Data/CellWalkRSampleData/SampleSnap.rds")
# # ATACMat = snapData@pmat
#
# #cell sim
# # cellEdges = computeCellSim(ATACMat, method="Jaccard")
#
# #regions
# regions = getRegions()
#
# #map peaks to genes
# labelGenes = data.table::fread("../CellWalkRSampleData/SampleMarkers1.txt")
# ATACGenePeak = mapPeaksToGenes(labelGenes, ATACMat, peaks, regions)
#
# #label weights
# labelEdges = computeLabelEdges(labelGenes, ATACMat, ATACGenePeak)
# labelEdgesList = list(labelEdges)
# edgeWeights = tuneEdgeWeights(cellEdges, labelEdgesList, labelEdgeOpts=10^seq(-2,7,1),
#                                       sampleDepth=1000)
# optWeight = 10000 #estimated log cell hom = 0.06618983
#
# filter = data.table::fread("../CellWalkRSampleData/SampleFilter.bed")
# filter = GRanges(filter$V1, IRanges(filter$V2, filter$V3))
# filters = list(filter)
# labelGenesList = list(labelGenes)
# filterWeights = tuneFilterWeights(cellEdges, labelGenesList, labelEdgesList, optWeight,
                                  # ATACMat, ATACGenePeak,
                                  # filters=filters, regions=regions, sampleDepth=1000)


# labelEdgesFilter = computeLabelEdges(labelGenes, ATACMat, ATACGenePeak,
#                                      filters = list(filter), filterWeights = 1, regions = regions)
# labelEdgesFilterList = list(labelEdgesFilter)
# edgeWeightsFilter = tuneEdgeWeights(cellEdges, labelEdgesFilterList, labelEdgeOpts=10^seq(-2,7,1),
#                                                           sampleDepth=1000)
# optWeight = 100000 #estimated log cell hom = 0.1249883
# cellWalk = walkCells(cellEdges, labelEdgesList, labelEdgeWeights = 1e+05)

# head(cellWalk$cellLabels)
# cellWalk = clusterLabels(cellWalk, plot=TRUE)
# cellWalk = findUncertainLabels(cellWalk, plot=TRUE)

# sampleEnhancers = data.table::fread("../CellWalkRSampleData/sampleEnhancers.bed")
# sampleEnhancers = GRanges(sampleEnhancers$V1, IRanges(sampleEnhancers$V2, sampleEnhancers$V3))
# mappedLabel = labelBulk(cellWalk, sampleEnhancers[1:60], ATACMat, peaks) #warning, slow

#BONUS
# labelGenesB = data.table::fread("../CellWalkRSampleData/SampleMarkers2.txt")
# ATACGenePeakB = mapPeaksToGenes(labelGenesB, ATACMat, peaks, regions)
# labelEdgesB = computeLabelEdges(labelGenesB, ATACMat, ATACGenePeakB)
# labelEdgesListB = list(labelEdges, labelEdgesB)
# edgeWeightsB = tuneEdgeWeights(cellEdges, labelEdgesListB, labelEdgeOpts=10^seq(-2,7,1),
#                                sampleDepth=1000)
# cellWalkB = walkCells(cellEdges, labelEdgesListB, labelEdgeWeights = c(1000, 1e+06))
#
# also splitting labels by cells?
# turn plotting into function?
# make code work with a gmat (should be its own section later)
# section headings: Using cell-by-gene matrices
