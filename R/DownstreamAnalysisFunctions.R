#' Cluster labels
#'
#' \code{clusterLabels()} Computes hierarchical clustering of labels
#'
#' @param cellWalk a cellWalk object
#' @param cellTypes character, vector of labels to use, all labels used by default
#' @param distMethod character, method used to compute distance
#' @param plot boolean, plot output matrix
#' @return cellWalk object with label clustering stored in "cluster"
#' @export
clusterLabels = function(cellWalk, cellTypes, distMethod="euclidean", plot=FALSE){
  if(missing(cellWalk) || !is(cellWalk, "cellWalk")){
    stop("Must provide a cellWalk object")
  }
  numLabels = dim(cellWalk[["normMat"]])[2]
  typeTypeInf = cellWalk[["infMat"]][1:numLabels,1:numLabels]
  colnames(typeTypeInf) = colnames(cellWalk[["normMat"]])
  rownames(typeTypeInf) = colnames(cellWalk[["normMat"]])

  if(missing(cellTypes)){
    cellTypes = colnames(cellWalk[["normMat"]])
  }
  keepLabels = which( colnames(cellWalk[["normMat"]]) %in% cellTypes)
  if(length(keepLabels)<2){
    stop("Fewer than two labels match cell types")
  }
  typeTypeInf = typeTypeInf[keepLabels,keepLabels]

  typeTypeClust = hclust(dist(typeTypeInf, method = distMethod))

  if(plot){
    print(plot(typeTypeClust, main="Hierarchical Clustering of Labels", sub=NA, xlab=NA))
  }

  cellWalk[["cluster"]] = typeTypeClust
  cellWalk
}


#' Find uncertain labels
#'
#' \code{findUncertainLabels()} generates an l-by-l matrix for how often each label is
#' confused for each other label at a given threshold, optinally normalized for total
#' counts of cells for each label (default FALSE)
#'
#' @param cellWalk a cellWalk object
#' @param cellTypes character, vector of labels to use, all labels used by default
#' @param threshold numeric, quantile threshold for uncertain label
#' @param labelThreshold numeric, set a threshold below which cells aren't labeled (e.g. 0)
#' @param plot boolean, plot output matrix
#' @param normalize boolean, normalize plot scale
#' @return cellWalk object with label uncertainty matrix (l-by-l) stored in "uncertaintyMatrix"
#' @export
findUncertainLabels = function(cellWalk, cellTypes, threshold=.1, labelThreshold, plot=FALSE, normalize=FALSE){
  if(missing(cellWalk) || !is(cellWalk, "cellWalk")){
    stop("Must provide a cellWalk object")
  }
  normMat = cellWalk[["normMat"]]

  if(missing(cellTypes)){
    cellTypes = colnames(normMat)
  }
  keepLabels = which( colnames(normMat) %in% cellTypes)
  if(length(keepLabels)<2){
    stop("Fewer than two labels match cell types")
  }
  cellTypes = cellTypes[cellTypes %in% colnames(normMat)]
  normMat = normMat[,match(cellTypes, colnames(normMat))]
  cellLabels = cellWalk[["cellLabels"]]

  if(!missing(labelThreshold)){
    cellLabels = cellLabels[apply(normMat, 1, max)>labelThreshold]
    normMat = normMat[apply(normMat, 1, max)>labelThreshold,]
  }

  #difference between top two label scores
  uncertaintyScore = apply(normMat, 1, function(x) sort(x, decreasing = TRUE)[1]-sort(x, decreasing = TRUE)[2])
  uncertaintyTypes = apply(normMat, 1, function(x) c(cellTypes[order(x, decreasing = TRUE)][1], cellTypes[order(x, decreasing = TRUE)][2]))
  #all pairs of labels where difference in scores below threshold
  uncertainPairs = uncertaintyTypes[,uncertaintyScore < quantile(uncertaintyScore, threshold)]
  uncertainMat = sapply(cellTypes, function(c1) sapply(cellTypes, function(c2) length(which((uncertainPairs[1,]==c1 & uncertainPairs[2,]==c2) | (uncertainPairs[1,]==c2 & uncertainPairs[2,]==c1)))))
  uncertainMatNorm = sapply(cellTypes, function(c1)
    sapply(cellTypes, function(c2)
      length(which((uncertainPairs[1,]==c1 & uncertainPairs[2,]==c2) |
                     (uncertainPairs[1,]==c2 & uncertainPairs[2,]==c1)))
      /length(which(cellLabels==c1 | cellLabels==c2))))

  if(plot){
    if(!requireNamespace("ggplot2", quietly = TRUE)){
      stop("Must install ggplot2")
    }
    if(!requireNamespace("reshape2", quietly = TRUE)){
      stop("Must install reshape2")
    }
    if(normalize){
      uncertainMatMelt = reshape2::melt(uncertainMatNorm)
      labelText = "Fraction of \nnearly equal \nscoring cells"
    }
    else{
      uncertainMatMelt = reshape2::melt(uncertainMat)
      labelText = "Number of \nnearly equal \nscoring cells"
    }
    uncertainMatMelt$Var1 = factor(uncertainMatMelt$Var1, levels = cellTypes)
    uncertainMatMelt$Var2 = factor(uncertainMatMelt$Var2, levels = cellTypes)
    print(ggplot2::ggplot(uncertainMatMelt, ggplot2::aes(Var1, Var2)) +
      ggplot2::geom_tile(ggplot2::aes(fill = value), color = "black") +
      ggplot2::scale_fill_gradient(low = "white",high = "steelblue") +
      ggplot2::labs(fill=labelText) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust=0),
                     axis.title.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(angle = 45),
                     axis.title.y = ggplot2::element_blank()) +
      ggplot2::scale_x_discrete(position = "top"))
  }

  cellWalk[["uncertaintyMatrix"]] = uncertainMat
  cellWalk

}

#' Plot Cells
#'
#' \code{plotCells()} generates a tSNE plot based on cell-to-cell influence
#'
#' @param cellWalk a cellWalk object
#' @param cellTypes character, optional vector of labels to use, all labels used by default. Most likely label among those listed
#'  will be used. If only a single label is provided, all cells will be colored by their score for that label. If two labels are
#'  given, the difference in score for each cell is shown.
#' @param labelThreshold numeric, set a threshold below which cells aren't labeled (e.g. 0)
#' @param embedding method with which to embed data, either "tSNE" or "UMAP"
#' @param initial_dims numeric, number of PCA dims to use for tSNE
#' @param perplexity numeric, perplexity parameter for tSNE
#' @param recompute boolean, recompute tSNE
#' @param plot boolean, optionally don't plot output and only compute the embedding
#' @param seed numeric, random seed
#' @return cellWalk object with embedding stored in "tSNE" or "UMAP"
#' @export
plotCells = function(cellWalk, cellTypes, labelThreshold, embedding = "tSNE", initial_dims = 10, perplexity = 50, recompute = FALSE, plot = TRUE, seed){
  if(missing(cellWalk) || !is(cellWalk, "cellWalk")){
    stop("Must provide a cellWalk object")
  }

  if(!missing(seed)){
    set.seed(seed)
  }

  numLabels = dim(cellWalk[["normMat"]])[2]
  cellCellInf = cellWalk[["infMat"]][-(1:numLabels),-(1:numLabels)]

  if(embedding=="tSNE"){
    if(!requireNamespace("Rtsne", quietly = TRUE)){
      stop("Must install 'Rtsne' to compute tSNE")
    }

    if(recompute | is.null(cellWalk[["tSNE"]])){
      celltSNE = Rtsne::Rtsne(cellCellInf, initial_dims = initial_dims, perplexity = perplexity)
      cellWalk[["tSNE"]] = celltSNE$Y
    }
    celltSNE = cellWalk[["tSNE"]]
  }else{
    if(!requireNamespace("uwot", quietly = TRUE)){
      stop("Must install 'uwot' to compute UMAP")
    }

    if(recompute | is.null(cellWalk[["UMAP"]])){
      celltSNE = uwot::umap(cellCellInf)
      cellWalk[["UMAP"]] = celltSNE
    }
    celltSNE = cellWalk[["UMAP"]]
  }

  plotColor = cellWalk$cellLabels
  if(!missing(labelThreshold)){
    plotColor[apply(cellWalk[["normMat"]], 1, max)<=labelThreshold] = "Other"
  }
  labelText = "Label"
  if(!missing(cellTypes)){
    cellTypes = cellTypes[cellTypes %in% colnames(cellWalk[["normMat"]])]
    if(length(cellTypes)==0){
      stop("No labels match cell types")
    }
    else if(length(cellTypes)==1){
      plotColor = cellWalk[["normMat"]][,cellTypes]
      labelText = paste(cellTypes,"Score")
    }
    else if(length(cellTypes)==2){
      plotColor = cellWalk[["normMat"]][,cellTypes[1]]-cellWalk[["normMat"]][,cellTypes[2]]
      labelText = paste0(cellTypes[1]," vs\n",cellTypes[2]," Score")
    }
    else{
      normMatTrim = cellWalk[["normMat"]][,cellTypes]
      plotColor = apply(normMatTrim, 1, function(x) cellTypes[order(x, decreasing = TRUE)][1])
      if(!missing(labelThreshold)){
        plotColor[apply(normMatTrim, 1, max)<=labelThreshold] = "Other"
      }
    }
  }

  if(plot){
    print(ggplot2::ggplot() +
            ggplot2::geom_point(ggplot2::aes(celltSNE[,1],celltSNE[,2], color=plotColor), size = 1) +
            ggplot2::xlab(paste0(embedding,"_1"))+
            ggplot2::ylab(paste0(embedding,"_2"))+
            ggplot2::labs(color=labelText)+
            ggplot2::theme_classic())
  }

  cellWalk
}

#' Comptute Minimum Spanning Tree
#'
#' \code{computeMST()} computes a minimum spanning tree based on cell-to-cell influence
#'
#' @param cellWalk a cellWalk object
#' @param cellTypes character, optional vector of labels to use, all labels used by default. Most likely label among those listed
#'  will be used. If only a single label is provided, all cells will be colored by their score for that label. If two labels are
#'  given, the difference in score for each cell is shown.
#' @param labelThreshold numeric, set a threshold below which cells aren't labeled (e.g. 0)
#' @param recompute boolean, recompute layout for MST plot
#' @param plot boolean, optionally don't plot output and only compute the embedding
#' @param seed numeric, random seed
#' @return cellWalk object with MST stored in "MST" and layout stored in "MST_layout"
#' @export
computeMST = function(cellWalk, cellTypes, labelThreshold, recompute = FALSE, plot = TRUE, seed){
  if(missing(cellWalk) || !is(cellWalk, "cellWalk")){
    stop("Must provide a cellWalk object")
  }

  if(!missing(seed)){
    set.seed(seed)
  }

  numLabels = dim(cellWalk[["normMat"]])[2]
  cellCellInf = cellWalk[["infMat"]][-(1:numLabels),-(1:numLabels)]

  if(!requireNamespace("igraph", quietly = TRUE)){
    stop("Must install igraph")
  }

  if(recompute | is.null(cellWalk[["MST"]]) | is.null(cellWalk[["MST_layout"]])){
    graph = igraph::graph_from_adjacency_matrix(cellCellInf, weighted = TRUE)
    MST = igraph::mst(graph, weights = max(igraph::E(graph)$weight)-igraph::E(graph)$weight)
    cellWalk[["MST"]] = MST
    cellWalk[["MST_layout"]] = igraph::layout_with_lgl(MST)
  }

  if(plot){
    plotColor = cellWalk$cellLabels
    if(!missing(labelThreshold)){
      plotColor[apply(cellWalk[["normMat"]], 1, max)<=labelThreshold] = "Other"
    }
    labelText = "Label"
    if(!missing(cellTypes)){
      cellTypes = cellTypes[cellTypes %in% colnames(cellWalk[["normMat"]])]
      if(length(cellTypes)==0){
        stop("No labels match cell types")
      }
      else if(length(cellTypes)==1){
        plotColor = cellWalk[["normMat"]][,cellTypes]
        labelText = paste(cellTypes,"Score")
      }
      else if(length(cellTypes)==2){
        plotColor = cellWalk[["normMat"]][,cellTypes[1]]-cellWalk[["normMat"]][,cellTypes[2]]
        labelText = paste0(cellTypes[1]," vs\n",cellTypes[2]," Score")
      }
      else{
        normMatTrim = cellWalk[["normMat"]][,cellTypes]
        plotColor = apply(normMatTrim, 1, function(x) cellTypes[order(x, decreasing = TRUE)][1])
        if(!missing(labelThreshold)){
          plotColor[apply(normMatTrim, 1, max)<=labelThreshold] = "Other"
        }
      }
    }
    else{cellTypes = colnames(cellWalk[["normMat"]])}
    if(length(cellTypes)==1){
      legendScores = round(quantile(plotColor, seq(0,1,1/5)), 2)
      plotColor = cut(plotColor, 100)
      print(igraph::plot.igraph(cellWalk$MST, layout=cellWalk$MST_layout, vertex.size=2, vertex.label=NA, edge.width=2, edge.arrow.size=.1, vertex.color=colorRampPalette(c("white","blue"))(100)[as.factor(plotColor)]))
      legend(title=labelText, x=par()$usr[2],y=par()$usr[4], legend = legendScores, bty = "n", fill=(colorRampPalette(c("white","blue"))(100))[seq(1,100,19)], xpd=TRUE)
    } else if(length(cellTypes)==2){
      legendScores = round(quantile(plotColor, seq(0,1,1/5)), 2)
      plotColor = cut(plotColor, 100)
      print(igraph::plot.igraph(cellWalk$MST, layout=cellWalk$MST_layout, vertex.size=2, vertex.label=NA, edge.width=2, edge.arrow.size=.1, vertex.color=colorRampPalette(c("red","blue"))(100)[as.factor(plotColor)]))
      legend(title=labelText, x=par()$usr[2],y=par()$usr[4], legend = legendScores, bty = "n", fill=(colorRampPalette(c("red","blue"))(100))[seq(1,100,19)], xpd=TRUE)
    }
    else{
      print(igraph::plot.igraph(cellWalk$MST, layout=cellWalk$MST_layout, vertex.size=2, vertex.label=NA, edge.width=2,
                                edge.arrow.size=.1, vertex.color=scales::hue_pal()(length(unique(plotColor)))[as.factor(plotColor)]))
      legend(title=labelText, x=par()$usr[2],y=par()$usr[4],legend=unique(plotColor), bty = "n", fill=scales::hue_pal()(length(unique(plotColor)))[as.factor(unique(plotColor))], xpd=TRUE)
    }
  }

  cellWalk
}


#' Label bulk data
#'
#' \code{labelBulk()} Determines labels for bulk data by mapping them via
#' the calculated information matrix
#'
#' @param cellWalk a cellWalk object
#' @param bulkPeaks GRanges of peaks in bulk data or GRangesList of sets of peaks
#' @param ATACMat cell-by-peak matrix
#' @param peaks GRanges of peaks in ATACMat
#' @param extendRegion GRanges defining where to extend mapping to consider peaks in a larger region (e.g. in LD) with bulk data
#' @param extendDistance numeric maximum distance to extend region by (if region is missing, just distance is used, if distance is missing whole region is used)
#' @param cellTypes character, vector of labels to use, all labels used by default
#' @param allScores return full table of scores
#' @param parallel execute in parallel
#' @param numCores number of cores to use for parallel execution
#' @return labels for each region in bulk data
#' @export
labelBulk = function(cellWalk, bulkPeaks, ATACMat, peaks, extendRegion, extendDistance, cellTypes, allScores=FALSE, parallel=FALSE, numCores=1){
  if(missing(cellWalk) || !is(cellWalk, "cellWalk")){
    stop("Must provide a cellWalk object")
  }
  infMat = cellWalk[["infMat"]]
  normMat = cellWalk[["normMat"]]

  if(missing(cellTypes)){
    cellTypes = colnames(normMat)
  }
  keepLabels = which( colnames(normMat) %in% cellTypes)
  if(length(keepLabels)<2){
    stop("Fewer than two labels match cell types")
  }
  cellTypes = cellTypes[cellTypes %in% colnames(normMat)]

  if(length(which(!colnames(normMat) %in% cellTypes))>0){
    infMat = infMat[-which(!colnames(normMat) %in% cellTypes),
                  -which(!colnames(normMat) %in% cellTypes)]
  }

  if(missing(peaks)){
    stop("Must provide a GRanges object of peaks")
  }
  if(missing(bulkPeaks) || !(is(bulkPeaks, "GRanges") | is(bulkPeaks, "GRangesList"))){
    stop("Must provide a GRanges or GRangesList object of peaks to map")
  }

  if(!missing(extendRegion)){
    if(!is(extendRegion, "GRanges")){
      stop("extendRegion must be a GRanges object")
    }
    whichExtends = GenomicRanges::findOverlaps(bulkPeaks, extendRegion)
    extensions = extendRegion[whichExtends@to]
    #restore original range if it was trimmed
    GenomicRanges::start(extensions) = sapply(1:length(extensions), function(x) min(GenomicRanges::start(extensions)[x], GenomicRanges::start(bulkPeaks[whichExtends@from[x]])))
    GenomicRanges::end(extensions) = sapply(1:length(extensions), function(x) max(GenomicRanges::end(extensions)[x], GenomicRanges::end(bulkPeaks[whichExtends@from[x]])))
    #cut extension to distance
    if(!missing(extendDistance)){
      GenomicRanges::start(extensions) = sapply(1:length(extensions), function(x) max(GenomicRanges::start(extensions)[x], GenomicRanges::start(bulkPeaks[whichExtends@from[x]])-extendDistance))
      GenomicRanges::end(extensions) = sapply(1:length(extensions), function(x) min(GenomicRanges::end(extensions)[x], GenomicRanges::end(bulkPeaks[whichExtends@from[x]])+extendDistance))
    }
    bulkPeaks[whichExtends@from] = extensions
  }

  if(missing(extendDistance) | !missing(extendRegion)){
    extendDistance = -1
  }

  peakOverlaps = GenomicRanges::findOverlaps(peaks, bulkPeaks, maxgap = extendDistance)

  if(parallel){
    if(!requireNamespace("parallel", quietly = TRUE)){
      stop("Must install parallel")
    }
    if(numCores >= parallel::detectCores()) {
      numCores = parallel::detectCores() - 1
      warning("numCores greater than avaialble cores, using available cores")
    }
    infCellOnType = parallel::mcmapply(function(e) {
      whichPeaks = peakOverlaps@from[peakOverlaps@to==e]
      if(length(whichPeaks)==0){
        rep(NA, length(cellTypes))
      } else{
        testCells = ATACMat[,whichPeaks]
        if(length(whichPeaks)==1){
          whichCells = which(testCells>0)
        } else{
          whichCells = which(Matrix::rowSums(testCells)>0)
        }
        if(length(whichCells)==1){
          infMat[1:length(cellTypes),length(cellTypes)+whichCells]
        }
        else{
          Matrix::rowSums(infMat[1:length(cellTypes),length(cellTypes)+whichCells])/length(whichCells)
        }
      }
    }, e=1:length(bulkPeaks), mc.cores = numCores)
  }
  else{
    infCellOnType = sapply(1:length(bulkPeaks), function(e) {
      whichPeaks = peakOverlaps@from[peakOverlaps@to==e]
      if(length(whichPeaks)==0){
        rep(NA, length(cellTypes))
      } else{
        testCells = ATACMat[,whichPeaks]
        if(length(whichPeaks)==1){
          whichCells = which(testCells>0)
        } else{
          whichCells = which(Matrix::rowSums(testCells)>0)
        }
        if(length(whichCells)==1){
          infMat[1:length(cellTypes),length(cellTypes)+whichCells]
        }
        else{
          Matrix::rowSums(infMat[1:length(cellTypes),length(cellTypes)+whichCells])/length(whichCells)
        }
      }
    })
  }

  if(allScores){
    rownames(infCellOnType) = cellTypes
    t(apply(infCellOnType, 2, function(x) (x-mean(x))/sd(x)))
  }
  else{
    mappedLabel = apply(infCellOnType, 2, function(x) ifelse(length(which(is.na(x)))==0,cellTypes[order(x, decreasing = TRUE)[1]],NA))
    mappedLabel
  }
}

#' Store bulk label data
#'
#' \code{storeBulk()} stores labels for bulk data in the cellWalk object
#'
#' @param cellWalk a cellWalk object
#' @param bulkPeaks GRanges of peaks in bulk data or GRangesList of sets of peaks
#' @param labelScores a matrix of label scores, output of labelBulk with allScores set to TRUE
#' @param bulkID character string assigning name to bulk data
#' @return cellWalk object with labels for each region in bulk data
#' @export
storeBulk = function(cellWalk, bulkPeaks, labelScores, bulkID){
  if(missing(cellWalk) || !is(cellWalk, "cellWalk")){
    stop("Must provide a cellWalk object")
  }
  if(missing(bulkPeaks) || !(is(bulkPeaks, "GRanges") | is(bulkPeaks, "GRangesList"))){
    stop("Must provide a GRanges or GRangesList object of bulk peaks")
  }
  if(missing(labelScores) || !is(labelScores, "matrix")){
    stop("Must provide a matrix of label scores (the output of labelBulk with allScores set to TRUE)")
  }
  if(dim(labelScores)[2] != dim(cellWalk$normMat)[2]){
    stop("label scores has wrong number of labels")
  }
  if(dim(labelScores)[1] != length(bulkPeaks)){
    stop("label scores has wrong number of mapped ranges")
  }

  if(missing(bulkID) || !is(bulkID, "character")){
    bulkID = deparse(substitute(bulkPeaks))
  }

  cellWalk[["labelBulk"]][[bulkID]] = list(bulkPeaks=bulkPeaks, labelScores=labelScores)

  cellWalk
}

#' Store ATACMat data
#'
#' \code{storeMat()} stores the original ATAC matrix
#'
#' @param cellWalk a cellWalk object
#' @param ATACMat cell-by-peak matrix
#' @param peaks GRanges of peaks in ATACMat
#' @return cellWalk object with original ATAC matrix stored
#' @export
storeMat = function(cellWalk, ATACMat, peaks){
  if(missing(cellWalk) || !is(cellWalk, "cellWalk")){
    stop("Must provide a cellWalk object")
  }
  if(missing(ATACMat)){
    stop("Must provide a cell-by-peak matrix")
  }
  if(missing(peaks)){
    stop("Must provide a GRanges object of peaks")
  }

  cellWalk[["ATACMat"]] = ATACMat
  cellWalk[["peaks"]] = peaks

  cellWalk
}

#' Select significant labels
#'
#' \code{selectLabels()} Determines labels for bulk data by mapping them via
#' the calculated information matrix
#'
#' @param labelScores a matrix of label scores, output of labelBulk with allScores set to TRUE
#' @param z numeric z-score threshold for significance
#' @return list of signficant labels for each region in bulk data
#' @export
selectLabels = function(labelScores, z=2){
  if(missing(labelScores) || !is(labelScores, "matrix")){
    stop("Must provide a matrix of label scores (the output of labelBulk with allScores set to TRUE)")
  }
  sigTypes = apply(labelScores, 1, function(x) colnames(labelScores)[x>z & !is.na(x)])
  sigTypes
}

#' Select significant labels across heirarchy
#'
#' \code{selectMultiLevelLabels()} Determines labels across heirarchy for bulk data by mapping them via
#' the calculated information matrix
#'
#' @param cellWalk a cellWalk object
#' @param labelScores a matrix of label scores, output of labelBulk with allScores set to TRUE
#' @param z numeric z-score threshold for significance
#' @return list of signficant labels for each region in bulk data
#' @export
selectMultiLevelLabels = function(cellWalk, labelScores, z=2){
  if(missing(cellWalk) || !is(cellWalk, "cellWalk")){
    stop("Must provide a cellWalk object")
  }
  if(missing(labelScores) || !is(labelScores, "matrix")){
    stop("Must provide a matrix of label scores (the output of labelBulk with allScores set to TRUE)")
  }

  sigTypes = apply(labelScores, 1, function(x) colnames(labelScores)[x>z & !is.na(x)])

  #check for cluster
  if(!any(names(cellWalk) == "cluster")){
    stop("Must first run clusterLabels on cellWalk")
  }

  #compute merge scores
  mergeScore = c()
  for(i in 1:dim(cellWalk$cluster$merge)[1]){
    a = if(cellWalk$cluster$merge[i,1]<0){labelScores[,-cellWalk$cluster$merge[i,1]]}else{mergeScore[,cellWalk$cluster$merge[i,1]]}
    b = if(cellWalk$cluster$merge[i,2]<0){labelScores[,-cellWalk$cluster$merge[i,2]]}else{mergeScore[,cellWalk$cluster$merge[i,2]]}
    c = (a+b)/2
    d = (2-cellWalk$cluster$height[i])*c
    mergeScore = cbind(mergeScore, d)
  }

  sigTypesLoc = sapply(1:length(sigTypes), function(x){
    as.numeric(which(mergeScore[x,]>z))
  })
  sigTypesList = sapply(1:length(sigTypes), function(x){
    sigGroups = rep(list(c()), length(sigTypesLoc[[x]]))
    if(length(sigTypesLoc[[x]])>0){
      for(i in 1:length(sigTypesLoc[[x]])){
        theseTypes = c()
        mergeList = sigTypesLoc[[x]][i]
        while(length(mergeList)>0) {
          thisMerge = mergeList[1]
          mergeList = mergeList[-1]
          if(thisMerge<0){
            theseTypes = c(theseTypes, colnames(labelScores)[-thisMerge])
          } else{
            mergeList = c(mergeList, cellWalk$cluster$merge[thisMerge,])
          }
        }
        sigGroups[[i]] = theseTypes
      }
    }
    sigGroups
  })
  sigTypesList

}

#' Plot significant labels across heirarchy
#'
#' \code{plotMultiLevelLabels()} Plots labels across heirarchy for bulk data by mapping them via
#' the calculated information matrix
#'
#' @param cellWalk a cellWalk object
#' @param labelScores a matrix of label scores, output of labelBulk with allScores set to TRUE
#' @param z numeric z-score threshold for significance
#' @param whichBulk, numeric optionally identify which bulk peak to plot
#' @return plot object p of either counts of signifcant labels or enrichment/depletion for single peak
#' @export
plotMultiLevelLabels = function(cellWalk, labelScores, z=2, whichBulk){
  if(missing(cellWalk) || !is(cellWalk, "cellWalk")){
    stop("Must provide a cellWalk object")
  }
  if(missing(labelScores) || !is(labelScores, "matrix")){
    stop("Must provide a matrix of label scores (the output of labelBulk with allScores set to TRUE)")
  }

  #check for cluster
  if(!any(names(cellWalk) == "cluster")){
    stop("Must first run clusterLabels on cellWalk")
  }

  #compute merge scores
  mergeScore = c()
  for(i in 1:dim(cellWalk$cluster$merge)[1]){
    a = if(cellWalk$cluster$merge[i,1]<0){labelScores[,-cellWalk$cluster$merge[i,1]]}else{mergeScore[,cellWalk$cluster$merge[i,1]]}
    b = if(cellWalk$cluster$merge[i,2]<0){labelScores[,-cellWalk$cluster$merge[i,2]]}else{mergeScore[,cellWalk$cluster$merge[i,2]]}
    c = (a+b)/2
    d = (2-cellWalk$cluster$height[i])*c
    mergeScore = cbind(mergeScore, d)
  }

  plotScores = mergeScore[,dim(cellWalk$cluster$merge)[1]]
  testList = cellWalk$cluster$merge[dim(cellWalk$cluster$merge)[1],]
  while(length(testList)>0){
    thisTest = testList[1]
    testList = testList[-1]
    if(thisTest<0){
      plotScores = cbind(plotScores, labelScores[,-thisTest])
    } else{
      plotScores = cbind(plotScores, mergeScore[,thisTest])
      testList = c(cellWalk$cluster$merge[thisTest,], testList)
    }
  }

  if(!requireNamespace("dendextend", quietly = TRUE)){
    stop("Must install dendextend")
  }
  if(!requireNamespace("circlize", quietly = TRUE)){
    stop("Must install circlize")
  }

  adjClust = cellWalk$cluster
  adjClust$height = adjClust$height-min(adjClust$height)
  if(missing(whichBulk)){
    weights = colSums(plotScores>z, na.rm=TRUE)/max(colSums(plotScores>z, na.rm=TRUE))
    weights[is.na(weights)] = 0
    # p = adjClust %>% as.dendrogram %>%
    #   dendextend::set("branches_lwd", weights*5+1) %>%
    #   dendextend::set("branches_col", colorRampPalette(c("gray", "#1F77B4"))(8)[cut(weights*5, breaks=c(-Inf,0,.5,1,1.5,2,2.5,3,Inf))])
    p = dendextend::set(dendextend::set(as.dendrogram(adjClust),
            what="branches_lwd", value=weights*5+1),
            what="branches_col", value=colorRampPalette(c("gray", "#1F77B4"))(8)[cut(weights*5, breaks=c(-Inf,0,.5,1,1.5,2,2.5,3,Inf))])

  }
  else{
    # p = adjClust %>% as.dendrogram %>%
    #   dendextend::set("branches_lwd", abs(plotScores[whichBulk,])) %>%
    #   dendextend::set("branches_col", colorRampPalette(c("red", "white", "blue"))(8)[cut(plotScores[whichBulk,], breaks=c(-Inf,-3,-2,-1,0,1,2,3,Inf))])
    p = dendextend::set(dendextend::set(as.dendrogram(adjClust),
      what="branches_lwd", abs(plotScores[whichBulk,])+2),
      what="branches_col", value=colorRampPalette(c("red", "white", "blue"))(8)[cut(plotScores[whichBulk,], breaks=c(-Inf,-3,-2,-1,0,1,2,3,Inf))])
  }

  print(dendextend::circlize_dendrogram(p, dend_track_height = .85))

  p
}
