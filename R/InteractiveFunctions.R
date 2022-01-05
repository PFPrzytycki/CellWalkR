#' Launch Interactve Visualization
#'
#' \code{launchViz()} launches an interactve visualization
#'
#' @param cellWalk a cellWalk object
#' @export
launchViz = function(cellWalk){

  if(!requireNamespace("shiny", quietly = TRUE)){
    stop("Must install shiny")
  }
  if(!requireNamespace("plotly", quietly = TRUE)){
    stop("Must install plotly")
  }

  options(shiny.maxRequestSize=500*1024^2)

  if(missing(cellWalk) || !is(cellWalk, "cellWalk")){
    stop("Must provide a cellWalk object")
  }

  cellTypes = colnames(cellWalk$normMat)

  if (interactive()) {

    ui <- fluidPage(
      titlePanel('CellWalkR'),
      sidebarLayout(
        sidebarPanel(
          div(h3("Visualize Results")),
          uiOutput("tSNE_box"),
          uiOutput("tSNE_opt"),
          tags$hr(),
          div(h3("Label Bulk Data")),
          uiOutput("bulk_select"),
          downloadButton("downloadLabel", "Save bulk label scores"),
          downloadButton("downloadSigLabel", "Save bulk cell types")
          ),
        mainPanel(
          tabsetPanel(id = "tabs",
            tabPanel("hClust", plotOutput("hClust", 600,600), tags$hr(),
                     "Using label-to-label influence, CelWalkR builds a hierachical clustering of cell types. This clustering can be used to visualize how bulk features map to the full heirarchy of cell types. Either select a GRanges object or upload a bed file and hit run to see how often those ranges map to each cell type. Select a single one of those ranges to see how that single featue scores for each cell type, with blue indicating  enrichment and red indicating depletion. All label scores for all bulk regions can be downloaded using 'Save label scores' and cell types for each bulk region can be downloaded using 'Save bulk cell types'",
                     tags$br(),downloadButton('savehClust', "Save image")),
            tabPanel("Confusion Matrix", plotOutput("cMat", 600,600), tags$hr(),
                     "Labels assigned by CellWalkR are 'fuzzy' meaning each cell actually has a distribution of scores from each label. Thanks to this, we can examine how often labels are confused for each other. This confusion matrix shows the number of cells that have nearly identical sores for each pair of labels.",
                     tags$br(),downloadButton('savecMat', "Save image"))
          )
        )
      )
    )

    server <- function(input, output, session) {

      if(!is.null(cellWalk[["tSNE"]])){
        output$tSNE_opt <- renderUI({
          list(checkboxGroupInput("cellType", label = "Cell types:", cellTypes, cellTypes, inline = TRUE),
               sliderInput("labelThreshold", "Label Threshold:",-1,1,0,.1, width='200px'))
        })
        if(!is.null(cellWalk[["UMAP"]])){
          prependTab("tabs", tabPanel("tSNE",plotlyOutput("tSNE"),plotlyOutput("UMAP"), tags$hr(),
                                      "To generate a two-dimensional embeddings of cells, CellWalkR directly embeds cell-to-cell influence scores. Importantly, this portion of the influence matrix is not used in establishing the labeling of cells. Thus this can serve as distinct way to explore how labels were distributed accross cells. Because of this, unlike with most scATAC-seq analysis pipelines, clusters observed in this embedding may not directly correspond to labels. This can help understand cell diversity, as well as assist in identifying rare cell types. By defeault, this plot shows all cell types and uses a threshold of 0 to label cells as not associated with any included cell type (other). Using the options on the left, a subset of cell types can be plotted, which aids with identifying rare cell populations. If just one or two cell types are selected, the plot will show the scores for that cell types, or the difference between scores for the two cell types."))
        }
        else{
          prependTab("tabs", tabPanel("tSNE",plotlyOutput("tSNE"),tags$hr(),
                                      "To generate a two-dimensional embeddings of cells, CellWalkR directly embeds cell-to-cell influence scores. Importantly, this portion of the influence matrix is not used in establishing the labeling of cells. Thus this can serve as distinct way to explore how labels were distributed accross cells. Because of this, unlike with most scATAC-seq analysis pipelines, clusters observed in this embedding may not directly correspond to labels. This can help understand cell diversity, as well as assist in identifying rare cell types. By defeault, this plot shows all cell types and uses a threshold of 0 to label cells as not associated with any included cell type (other). Using the options on the left, a subset of cell types can be plotted, which aids with identifying rare cell populations. If just one or two cell types are selected, the plot will show the scores for that cell types, or the difference between scores for the two cell types. For comparison: ",
                                      actionButton("renderUMAP", "Generate UMAP")))
        }
        updateTabsetPanel(session, "tabs", selected = "tSNE")
      } else{
        output$tSNE_box <- renderUI({
          actionButton("rendertSNE", "Generate tSNE")
        })
      }

      if(!is.null(cellWalk[["MST"]])){
        output$tSNE_opt <- renderUI({
          list(checkboxGroupInput("cellType", label = "Cell types:", cellTypes, cellTypes, inline = TRUE),
               sliderInput("labelThreshold", "Label Threshold:",-1,1,0,.1, width='200px'))
        })
        appendTab("tabs", tabPanel("MST", plotOutput("MST", 600, 600), tags$hr(),
                                   "An alternative to two-dimensional embeddings of cells is to generate a Minimum Spanning Tree of the cell-to-cell influence scores. This way of plotting the underlying graph can help emphasize distinct cell types. By defeault, this plot shows all cell types and uses a threshold of 0 to label cells as not associated with any included cell type (other). Using the options on the left, a subset of cell types can be plotted, which aids with identifying rare cell populations. If just one or two cell types are selected, the plot will show the scores for that cell types, or the difference between scores for the two cell types.",
                                   tags$br(),downloadButton("saveMST", "Save image")))
      } else{
        output$tSNE_box <- renderUI({
          actionButton("renderMST", "Generate MST")
        })
      }

      if(is.null(cellWalk[["tSNE"]]) & is.null(cellWalk[["MST"]])){
        output$tSNE_box <- renderUI({
          list(actionButton("rendertSNE", "Generate tSNE"),
          actionButton("renderMST", "Generate MST"))
        })
      }

      observe({
        if(is.null(input$bulkMap)){}
        if(length(cellWalk$labelBulk)>0){
        if(!is.null(input$bulkMap)){
        validRegions = unname(as.character(cellWalk$labelBulk[[input$bulkMap]][["bulkPeaks"]]))
        validRegions = validRegions[!is.na(rowSums(cellWalk$labelBulk[[input$bulkMap]][["labelScores"]]))]
        output$bulk_select <- renderUI({
          list(selectInput("bulkMap", "Which data: ", names(cellWalk$labelBulk), selected = input$bulkMap, width='200px'),
               selectInput("whichBulk", "Which mapped region: ", c("All", validRegions) , width='200px'),
               tags$hr(),
               fixedRow(column(width = 9, selectInput("bedRanges","Label new GRanges:", unlist(sapply(objects(1), function(x) if("GRanges" %in% is(get(x))){x}))),
                               fileInput("bedFile","Label a bed file:", accept = ".bed")),
                        column(width = 3,tags$br(), actionButton("runRanges","Run"),tags$br(),tags$br(),tags$br(),
                               actionButton("runBed","Run"))))
        })
        }else{
        validRegions = unname(as.character(cellWalk$labelBulk[[names(cellWalk$labelBulk)[1]]][["bulkPeaks"]]))
        validRegions = validRegions[!is.na(rowSums(cellWalk$labelBulk[[names(cellWalk$labelBulk)[1]]][["labelScores"]]))]
        output$bulk_select <- renderUI({
          list(selectInput("bulkMap", "Which data: ", names(cellWalk$labelBulk), width='200px'),
               selectInput("whichBulk", "Which mapped region: ", c("All", validRegions) , width='200px'),
               tags$hr(),
               fixedRow(column(width = 9, selectInput("bedRanges","Label new GRanges:", unlist(sapply(objects(1), function(x) if("GRanges" %in% is(get(x))){x}))),
                               fileInput("bedFile","Label a bed file:", accept = ".bed")),
                        column(width = 3,tags$br(), actionButton("runRanges","Run"),tags$br(),tags$br(),tags$br(),
                               actionButton("runBed","Run"))))
        })
        }
      } else{
        output$bulk_select <- renderUI({
         fixedRow(column(width = 9, selectInput("bedRanges","Label new GRanges:", unlist(sapply(objects(1), function(x) if("GRanges" %in% is(get(x))){x}))),
                         fileInput("bedFile","Label a bed file:", accept = ".bed")),
                  column(width = 3,tags$br(), actionButton("runRanges","Run"),tags$br(),tags$br(),tags$br(),
                         actionButton("runBed","Run")))
        })
      }
      })

      observeEvent(input$runRanges, {
        if(is.null(cellWalk$ATACMat) | is.null(cellWalk$peaks)){
          showNotification("cellWalk needs ATACMat and peaks, can be added using storeMat()")
        } else{
        withProgress(message = "Labeling GRanges", {
          labelScores = labelBulk(cellWalk, get(input$bedRanges), cellWalk$ATACMat, cellWalk$peaks, allScores = TRUE) #need to store these
        })
        cellWalk <<- storeBulk(cellWalk, get(input$bedRanges), labelScores, input$bedRanges)
        if(!is.null(input$bulkMap)){
          validRegions = unname(as.character(cellWalk$labelBulk[[input$bedRanges]][["bulkPeaks"]]))
          validRegions = validRegions[!is.na(rowSums(cellWalk$labelBulk[[input$bedRanges]][["labelScores"]]))]
          output$bulk_select <- renderUI({
            list(selectInput("bulkMap", "Which data: ", names(cellWalk$labelBulk), selected = input$bedRanges, width='200px'),
                 selectInput("whichBulk", "Which mapped region: ", c("All", validRegions) , width='200px'),
                 tags$hr(),
                 fixedRow(column(width = 9, selectInput("bedRanges","Label new GRanges:", unlist(sapply(objects(1), function(x) if("GRanges" %in% is(get(x))){x}))),
                                 fileInput("bedFile","Label a bed file:", accept = ".bed")),
                          column(width = 3,tags$br(), actionButton("runBed","Run"),tags$br(),tags$br(),tags$br(),
                                 actionButton("runRanges","Run"))))
          })
        }else{
          validRegions = unname(as.character(cellWalk$labelBulk[[names(cellWalk$labelBulk)[1]]][["bulkPeaks"]]))
          validRegions = validRegions[!is.na(rowSums(cellWalk$labelBulk[[names(cellWalk$labelBulk)[1]]][["labelScores"]]))]
          output$bulk_select <- renderUI({
            list(selectInput("bulkMap", "Which data: ", names(cellWalk$labelBulk), width='200px'),
                 selectInput("whichBulk", "Which mapped region: ", c("All", validRegions) , width='200px'),
                 tags$hr(),
                 fixedRow(column(width = 9, selectInput("bedRanges","Label new GRanges:", unlist(sapply(objects(1), function(x) if("GRanges" %in% is(get(x))){x}))),
                                 fileInput("bedFile","Label a bed file:", accept = ".bed")),
                          column(width = 3,tags$br(), actionButton("runRanges","Run"),tags$br(),tags$br(),tags$br(),
                                 actionButton("runBed","Run"))))
          })
        }
        updateTabsetPanel(session, "tabs", selected = "hClust")
        }
      })

      observeEvent(input$rendertSNE, {
        withProgress(message = "Generating tSNE", {
          cellWalk <<- plotCells(cellWalk, labelThreshold = 0, plot=FALSE)
        })
        if(!is.null(cellWalk[["MST"]])){
          output$tSNE_box <- renderUI({})
        }else{
          output$tSNE_box <- renderUI({
            actionButton("renderMST", "Generate MST")
          })
        }
        output$tSNE_opt <- renderUI({
          list(checkboxGroupInput("cellType", label = "Cell types:", cellTypes, cellTypes, inline = TRUE),
               sliderInput("labelThreshold", "Label Threshold:",-1,1,0,.1, width='200px'))
        })
        if(!is.null(cellWalk[["UMAP"]])){
          prependTab("tabs", tabPanel("tSNE",plotlyOutput("tSNE"),plotlyOutput("UMAP"), tags$hr(),
                                      "To generate a two-dimensional embeddings of cells, CellWalkR directly embeds cell-to-cell influence scores. Importantly, this portion of the influence matrix is not used in establishing the labeling of cells. Thus this can serve as distinct way to explore how labels were distributed accross cells. Because of this, unlike with most scATAC-seq analysis pipelines, clusters observed in this embedding may not directly correspond to labels. This can help understand cell diversity, as well as assist in identifying rare cell types. By defeault, this plot shows all cell types and uses a threshold of 0 to label cells as not associated with any included cell type (other). Using the options on the left, a subset of cell types can be plotted, which aids with identifying rare cell populations. If just one or two cell types are selected, the plot will show the scores for that cell types, or the difference between scores for the two cell types."))
        }
        else{
          prependTab("tabs", tabPanel("tSNE",plotlyOutput("tSNE"),tags$hr(),
                                      "To generate a two-dimensional embeddings of cells, CellWalkR directly embeds cell-to-cell influence scores. Importantly, this portion of the influence matrix is not used in establishing the labeling of cells. Thus this can serve as distinct way to explore how labels were distributed accross cells. Because of this, unlike with most scATAC-seq analysis pipelines, clusters observed in this embedding may not directly correspond to labels. This can help understand cell diversity, as well as assist in identifying rare cell types. By defeault, this plot shows all cell types and uses a threshold of 0 to label cells as not associated with any included cell type (other). Using the options on the left, a subset of cell types can be plotted, which aids with identifying rare cell populations. If just one or two cell types are selected, the plot will show the scores for that cell types, or the difference between scores for the two cell types. For comparison: ",
                                      actionButton("renderUMAP", "Generate UMAP")))
        }
        updateTabsetPanel(session, "tabs", selected = "tSNE")
      })

      observeEvent(input$renderUMAP, {
        withProgress(message = "Generating UMAP", {
          cellWalk <<- plotCells(cellWalk, embedding = "UMAP", labelThreshold = 0, plot = FALSE)
        })
        if(!is.null(cellWalk[["MST"]])){
          output$tSNE_box <- renderUI({})
        }else{
          output$tSNE_box <- renderUI({
            actionButton("renderMST", "Generate MST")
          })
        }
        output$tSNE_opt <- renderUI({
          list(checkboxGroupInput("cellType", label = "Cell types:", cellTypes, cellTypes, inline = TRUE),
               sliderInput("labelThreshold", "Label Threshold:",-1,1,0,.1, width='200px'))
        })
        prependTab("tabs", tabPanel("UMAP",plotlyOutput("UMAP"), tags$hr(),
                                    "To generate a two-dimensional embeddings of cells, CellWalkR directly embeds cell-to-cell influence scores. Importantly, this portion of the influence matrix is not used in establishing the labeling of cells. Thus this can serve as distinct way to explore how labels were distributed accross cells. Because of this, unlike with most scATAC-seq analysis pipelines, clusters observed in this embedding may not directly correspond to labels. This can help understand cell diversity, as well as assist in identifying rare cell types. By defeault, this plot shows all cell types and uses a threshold of 0 to label cells as not associated with any included cell type (other). Using the options on the left, a subset of cell types can be plotted, which aids with identifying rare cell populations. If just one or two cell types are selected, the plot will show the scores for that cell types, or the difference between scores for the two cell types."))
        removeTab("tabs","tSNE")
        prependTab("tabs", tabPanel("tSNE",plotlyOutput("tSNE"), tags$hr(),
                                    "To generate a two-dimensional embeddings of cells, CellWalkR directly embeds cell-to-cell influence scores. Importantly, this portion of the influence matrix is not used in establishing the labeling of cells. Thus this can serve as distinct way to explore how labels were distributed accross cells. Because of this, unlike with most scATAC-seq analysis pipelines, clusters observed in this embedding may not directly correspond to labels. This can help understand cell diversity, as well as assist in identifying rare cell types. By defeault, this plot shows all cell types and uses a threshold of 0 to label cells as not associated with any included cell type (other). Using the options on the left, a subset of cell types can be plotted, which aids with identifying rare cell populations. If just one or two cell types are selected, the plot will show the scores for that cell types, or the difference between scores for the two cell types."))
        updateTabsetPanel(session, "tabs", selected = "UMAP")
      })

      observeEvent(input$renderMST, {
        withProgress(message = "Generating MST", {
          cellWalk <<- computeMST(cellWalk, plot=FALSE)
        })
        if(!is.null(cellWalk[["tSNE"]])){
          output$tSNE_box <- renderUI({})
        }else{
          output$tSNE_box <- renderUI({
            actionButton("rendertSNE", "Generate tSNE")
          })
        }
        output$tSNE_opt <- renderUI({
          list(checkboxGroupInput("cellType", label = "Cell types:", cellTypes, cellTypes, inline = TRUE),
               sliderInput("labelThreshold", "Label Threshold:",-1,1,0,.1, width='200px'))
        })
        appendTab("tabs", tabPanel("MST", plotOutput("MST", 600,600), tags$hr(),
                                   "An alternative to two-dimensional embeddings of cells is to generate a Minimum Spanning Tree of the cell-to-cell influence scores. This way of plotting the underlying graph can help emphasize distinct cell types. By defeault, this plot shows all cell types and uses a threshold of 0 to label cells as not associated with any included cell type (other). Using the options on the left, a subset of cell types can be plotted, which aids with identifying rare cell populations. If just one or two cell types are selected, the plot will show the scores for that cell types, or the difference between scores for the two cell types.",
                                   tags$br(),downloadButton("saveMST", "Save image")))
        updateTabsetPanel(session, "tabs", selected = "MST")
      })

      output$tSNE = renderPlotly({
        theseTypes = input$cellType
        labelText = "Label"
        if(length(theseTypes)==1){
          plotColor = cellWalk[["normMat"]][,theseTypes]
          labelText = paste(theseTypes,"Score")
        } else if(length(theseTypes)==2){
          plotColor = cellWalk[["normMat"]][,theseTypes[1]]-cellWalk[["normMat"]][,theseTypes[2]]
          labelText = paste0(theseTypes[1]," vs\n",theseTypes[2]," Score")
        } else{
          normMatTrim = cellWalk[["normMat"]][,theseTypes]
          plotColor = apply(normMatTrim, 1, function(x) theseTypes[order(x, decreasing = TRUE)][1])
          labelThreshold = as.numeric(input$labelThreshold)
          plotColor[apply(normMatTrim, 1, max)<=labelThreshold] = "Other"
        }
        plotTable = data.frame(tSNE_1=cellWalk$tSNE[,1], tSNE_2=cellWalk$tSNE[,2], label=plotColor, score=apply(cellWalk$normMat, 1, max))
        if(length(theseTypes)==1){
          p = ggplot(plotTable) + geom_point(aes(tSNE_1, tSNE_2, color=label, alpha=score), size=1) + labs(color=labelText) + theme_classic() + guides(alpha="none") +
            scale_color_gradient(low="gray")
        } else if(length(theseTypes)==2){
          p = ggplot(plotTable) + geom_point(aes(tSNE_1, tSNE_2, color=label, alpha=score), size=1) + labs(color=labelText) + theme_classic() + guides(alpha="none") +
            scale_color_gradient(low="red", high = "blue")
        }
        else{
          p = ggplot(plotTable) + geom_point(aes(tSNE_1, tSNE_2, color=label, alpha=score), size=1) + labs(color=labelText) + theme_classic() + guides(alpha="none")
        }
        ggplotly(p)
      })

      output$UMAP = renderPlotly({
        theseTypes = input$cellType
        labelText = "Label"
        if(length(theseTypes)==1){
          plotColor = cellWalk[["normMat"]][,theseTypes]
          labelText = paste(theseTypes,"Score")
        } else if(length(theseTypes)==2){
          plotColor = cellWalk[["normMat"]][,theseTypes[1]]-cellWalk[["normMat"]][,theseTypes[2]]
          labelText = paste0(theseTypes[1]," vs\n",theseTypes[2]," Score")
        } else{
          normMatTrim = cellWalk[["normMat"]][,theseTypes]
          plotColor = apply(normMatTrim, 1, function(x) theseTypes[order(x, decreasing = TRUE)][1])
          labelThreshold = as.numeric(input$labelThreshold)
          plotColor[apply(normMatTrim, 1, max)<=labelThreshold] = "Other"
        }
        plotTable = data.frame(UMAP_1=cellWalk$UMAP[,1], UMAP_2=cellWalk$UMAP[,2], label=plotColor, score=apply(cellWalk$normMat, 1, max))
        if(length(theseTypes)==1){
          p = ggplot(plotTable) + geom_point(aes(UMAP_1, UMAP_2, color=label, alpha=score), size=1) + labs(color=labelText) + theme_classic() + guides(alpha="none") +
            scale_color_gradient(low="gray")
        } else if(length(theseTypes)==2){
          p = ggplot(plotTable) + geom_point(aes(UMAP_1, UMAP_2, color=label, alpha=score), size=1) + labs(color=labelText) + theme_classic() + guides(alpha="none") +
            scale_color_gradient(low="red", high = "blue")
        }
        else{
          p = ggplot(plotTable) + geom_point(aes(UMAP_1, UMAP_2, color=label, alpha=score), size=1) + labs(color=labelText) + theme_classic() + guides(alpha="none")
        }
        ggplotly(p)
      })

      output$hClust = renderPlot({
        if(is.null(cellWalk[["cluster"]])){
          cellWalk <<- clusterLabels(cellWalk)
        }
        if(is.null(input$whichBulk)){}
        if(length(cellWalk$labelBulk)==0){
          plotMultiLevelLabels(cellWalk, matrix(0, nrow = 1, ncol=length(cellWalk$cellLabels)))
        }else{
          if(input$whichBulk!="All"){
            whichBulk = which(cellWalk$labelBulk[[input$bulkMap]]$bulkPeaks ==  input$whichBulk)
            plotMultiLevelLabels(cellWalk, cellWalk$labelBulk[[input$bulkMap]][['labelScores']], whichBulk = whichBulk)
          } else{
            plotMultiLevelLabels(cellWalk, cellWalk$labelBulk[[input$bulkMap]][['labelScores']])
          }
        }
      }, width = 600, height = 600)

      output$cMat = renderPlot({
        labelThreshold = as.numeric(input$labelThreshold)
        if(length(labelThreshold)==0){
          findUncertainLabels(cellWalk, plot = TRUE)
        }
        else{
          findUncertainLabels(cellWalk, labelThreshold = labelThreshold, plot = TRUE)
        }
      }, width = 600, height = 600)

      output$MST = renderPlot({
        theseTypes = input$cellType
        # if(is.null(theseTypes)){theseTypes=cellTypes}
        labelText = "Label"
        if(length(theseTypes)==1){
          plotColor = cellWalk[["normMat"]][,theseTypes]
          labelText = paste(theseTypes,"Score")
        } else if(length(theseTypes)==2){
          plotColor = cellWalk[["normMat"]][,theseTypes[1]]-cellWalk[["normMat"]][,theseTypes[2]]
          labelText = paste0(theseTypes[1]," vs\n",theseTypes[2]," Score")
        } else{
          normMatTrim = cellWalk[["normMat"]][,theseTypes]
          plotColor = apply(normMatTrim, 1, function(x) theseTypes[order(x, decreasing = TRUE)][1])
          labelThreshold = as.numeric(input$labelThreshold)
          plotColor[apply(normMatTrim, 1, max)<=labelThreshold] = "Other"
        }
        if(length(theseTypes)==1){
          legendScores = round(quantile(plotColor, seq(0,1,1/5)), 2)
          plotColor = cut(plotColor, 100)
          par(mar=c(5,4,4,6))
          igraph::plot.igraph(cellWalk$MST, layout=cellWalk$MST_layout, vertex.size=2, vertex.label=NA, edge.width=2, edge.arrow.size=.1, vertex.color=colorRampPalette(c("white","blue"))(100)[as.factor(plotColor)])
          legend(title=labelText, x=par()$usr[2],y=par()$usr[4], legend = legendScores, bty = "n", fill=(colorRampPalette(c("white","blue"))(100))[seq(1,100,19)], xpd=TRUE)
        } else if(length(theseTypes)==2){
          legendScores = round(quantile(plotColor, seq(0,1,1/5)), 2)
          plotColor = cut(plotColor, 100)
          par(mar=c(5,4,4,6))
          igraph::plot.igraph(cellWalk$MST, layout=cellWalk$MST_layout, vertex.size=2, vertex.label=NA, edge.width=2, edge.arrow.size=.1, vertex.color=colorRampPalette(c("red","blue"))(100)[as.factor(plotColor)])
          legend(title=labelText, x=par()$usr[2],y=par()$usr[4], legend = legendScores, bty = "n", fill=(colorRampPalette(c("red","blue"))(100))[seq(1,100,19)], xpd=TRUE)
        }
        else{
          par(mar=c(5,4,4,6))
          igraph::plot.igraph(cellWalk$MST, layout=cellWalk$MST_layout, vertex.size=2, vertex.label=NA, edge.width=2, edge.arrow.size=.1, vertex.color=scales::hue_pal()(length(unique(plotColor)))[as.factor(plotColor)])
          legend(title=labelText, x=par()$usr[2],y=par()$usr[4],legend=unique(plotColor), bty = "n", fill=scales::hue_pal()(length(unique(plotColor)))[as.factor(unique(plotColor))], xpd=TRUE)
        }
      }, width = 600, height = 600)

      observeEvent(input$runBed, {
        bedFile = input$bedFile
        validate(need(!is.null(bedFile), "Please select a bed file"))

        bedData = data.table::fread(bedFile$datapath, header=FALSE)
        bedRanges = GRanges(bedData$V1, IRanges(bedData$V2, bedData$V3))

        if(is.null(cellWalk$ATACMat) | is.null(cellWalk$peaks)){
          showNotification("cellWalk needs ATACMat and peaks, can be added using storeMat()")
        } else{
        withProgress(message = "Labeling bed file", {
          labelScores = labelBulk(cellWalk, bedRanges, cellWalk$ATACMat, cellWalk$peaks, allScores = TRUE) #need to store these
        })
        cellWalk <<- storeBulk(cellWalk, bedRanges, labelScores, basename(bedFile$datapath))

        if(!is.null(input$bulkMap)){
          validRegions = unname(as.character(cellWalk$labelBulk[[basename(bedFile$datapath)]][["bulkPeaks"]]))
          validRegions = validRegions[!is.na(rowSums(cellWalk$labelBulk[[basename(bedFile$datapath)]][["labelScores"]]))]
          output$bulk_select <- renderUI({
            list(selectInput("bulkMap", "Which data: ", names(cellWalk$labelBulk), selected = basename(bedFile$datapath), width='200px'),
                 selectInput("whichBulk", "Which mapped region: ", c("All", validRegions) , width='200px'),
                 tags$hr(),
                 fixedRow(column(width = 9, selectInput("bedRanges","Label new GRanges:", unlist(sapply(objects(1), function(x) if("GRanges" %in% is(get(x))){x}))),
                                 fileInput("bedFile","Label a bed file:", accept = ".bed")),
                          column(width = 3,tags$br(), actionButton("runBed","Run"),tags$br(),tags$br(),tags$br(),
                                 actionButton("runRanges","Run"))))
          })
        }else{
          validRegions = unname(as.character(cellWalk$labelBulk[[names(cellWalk$labelBulk)[1]]][["bulkPeaks"]]))
          validRegions = validRegions[!is.na(rowSums(cellWalk$labelBulk[[names(cellWalk$labelBulk)[1]]][["labelScores"]]))]
          output$bulk_select <- renderUI({
            list(selectInput("bulkMap", "Which data: ", names(cellWalk$labelBulk), width='200px'),
                 selectInput("whichBulk", "Which mapped region: ", c("All", validRegions) , width='200px'),
                 tags$hr(),
                 fixedRow(column(width = 9, selectInput("bedRanges","Label new GRanges:", unlist(sapply(objects(1), function(x) if("GRanges" %in% is(get(x))){x}))),
                                 fileInput("bedFile","Label a bed file:", accept = ".bed")),
                          column(width = 3,tags$br(), actionButton("runRanges","Run"),tags$br(),tags$br(),tags$br(),
                                 actionButton("runBed","Run"))))
          })
        }
        updateTabsetPanel(session, "tabs", selected = "hClust")
        }
      })

      output$downloadLabel = downloadHandler(
        filename = function() {paste0(input$bulkMap, '-labelScores.csv')},
        content = function(file) {
          if(!is.null(input$bulkMap)){
            regions = unname(as.character(cellWalk$labelBulk[[input$bulkMap]][["bulkPeaks"]]))
            labelScores = cellWalk$labelBulk[[input$bulkMap]][["labelScores"]]
            writeFrame = data.frame(regions, labelScores)
          } else{writeFrame = data.frame()}
          write.csv(writeFrame, file, row.names=FALSE, quote=FALSE)
        }
      )

      output$downloadSigLabel = downloadHandler(
        filename = function() {paste0(input$bulkMap, '-sigLabels.csv')},
        content = function(file) {
          if(!is.null(input$bulkMap)){
            regions = unname(as.character(cellWalk$labelBulk[[input$bulkMap]][["bulkPeaks"]]))
            labelScores = cellWalk$labelBulk[[input$bulkMap]][["labelScores"]]
            mappedLabel = selectLabels(labelScores, z = 1.5)
            mappedLabelKeep = sapply(mappedLabel, function(x) length(x)>0)
            sigTypes = mappedLabel[mappedLabelKeep]
            regionsTable = data.frame(regions=as.character(regions[mappedLabelKeep]))
            regionsTable = data.frame(regions=regionsTable[rep(1:length(sigTypes), times=sapply(sigTypes, function(x) length(x))),])
            regionsTable$Type = unlist(sigTypes)
          } else{regionsTable = data.frame()}
          write.csv(regionsTable, file, row.names=FALSE, quote=FALSE)
        }
      )

      output$savecMat = downloadHandler(
        filename = function() {'Confusion_Matrix.png'},
        content = function(file) {
          png(file, 1200, 800)
          cellWalk = findUncertainLabels(cellWalk, cellTypes = input$cellType, plot = TRUE)
          dev.off()
        }
      )

      output$savehClust = downloadHandler(
        filename = function() {'hClust.png'},
        content = function(file) {
          png(file, 1200, 800)
          if(is.null(cellWalk[["cluster"]])){
            cellWalk <<- clusterLabels(cellWalk)
          }
          if(is.null(input$whichBulk)){}
          if(length(cellWalk$labelBulk)==0){
            plotMultiLevelLabels(cellWalk, matrix(0, nrow = 1, ncol=length(cellWalk$cellLabels)))
          }else{
            if(input$whichBulk!="All"){
              whichBulk = which(cellWalk$labelBulk[[input$bulkMap]]$bulkPeaks ==  input$whichBulk)
              plotMultiLevelLabels(cellWalk, cellWalk$labelBulk[[input$bulkMap]][['labelScores']], whichBulk = whichBulk)
            } else{
              plotMultiLevelLabels(cellWalk, cellWalk$labelBulk[[input$bulkMap]][['labelScores']])
            }
          }
          dev.off()
        }
      )

      output$saveMST = downloadHandler(
        filename = function() {'MST.png'},
        content = function(file) {
          png(file, 1200, 800)
          cellWalk = computeMST(cellWalk, cellTypes = input$cellType, labelThreshold = input$labelThreshold)
          dev.off()
        }
      )

    }

    shinyApp(ui, server)
  }
}
