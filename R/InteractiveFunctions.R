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
  require("shiny")
  if(!requireNamespace("plotly", quietly = TRUE)){
    stop("Must install plotly")
  }
  require("plotly")

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
          div(h3("Visualze Results")),
          uiOutput("tSNE_box"),
          uiOutput("tSNE_opt"),
          tags$hr(),
          div(h3("Label Bulk Data")),
          uiOutput("bulk_select"),
          downloadButton("downloadLabel", "Save label scores")
          ),
        mainPanel(
          tabsetPanel(id = "tabs",
            tabPanel("hClust", plotOutput("hClust")),
            tabPanel("Confusion Matrix", plotOutput("cMat"))
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
        prependTab("tabs", tabPanel("tSNE",plotlyOutput("tSNE")))
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
        appendTab("tabs", tabPanel("MST", plotOutput("MST")))
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
          cellWalk <<- plotCells(cellWalk, labelThreshold = 0, seed = 1)
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
        prependTab("tabs", tabPanel("tSNE",plotlyOutput("tSNE")))
        updateTabsetPanel(session, "tabs", selected = "tSNE")
      })

      observeEvent(input$renderMST, {
        withProgress(message = "Generating MST", {
          cellWalk <<- computeMST(cellWalk)
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
        appendTab("tabs", tabPanel("MST", plotOutput("MST")))
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
          p = ggplot(plotTable) + geom_point(aes(tSNE_1, tSNE_2, color=label, alpha=score), size=1) + labs(color=labelText) + theme_classic() + guides(alpha=FALSE) +
            scale_color_gradient(low="gray")
        } else if(length(theseTypes)==2){
          p = ggplot(plotTable) + geom_point(aes(tSNE_1, tSNE_2, color=label, alpha=score), size=1) + labs(color=labelText) + theme_classic() + guides(alpha=FALSE) +
            scale_color_gradient(low="red", high = "blue")
        }
        else{
          p = ggplot(plotTable) + geom_point(aes(tSNE_1, tSNE_2, color=label, alpha=score), size=1) + labs(color=labelText) + theme_classic() + guides(alpha=FALSE)
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
          plotColor = cut(plotColor, 100)
          igraph::plot.igraph(cellWalk$MST, layout=cellWalk$MST_layout, vertex.size=2, vertex.label=NA, edge.arrow.size=.1, vertex.color=colorRampPalette(c("white","blue"))(100)[as.factor(plotColor)])
        } else if(length(theseTypes)==2){
          plotColor = cut(plotColor, 100)
          igraph::plot.igraph(cellWalk$MST, layout=cellWalk$MST_layout, vertex.size=2, vertex.label=NA, edge.arrow.size=.1, vertex.color=colorRampPalette(c("red","blue"))(100)[as.factor(plotColor)])
        }
        else{
          igraph::plot.igraph(cellWalk$MST, layout=cellWalk$MST_layout, vertex.size=2, vertex.label=NA, edge.arrow.size=.1, vertex.color=scales::hue_pal()(length(unique(plotColor)))[as.factor(plotColor)])
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

    }

    shinyApp(ui, server)
  }
}
