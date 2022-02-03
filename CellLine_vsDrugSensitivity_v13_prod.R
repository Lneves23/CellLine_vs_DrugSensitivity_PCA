library("shiny")
library("gplots")
library("reshape2")
library("coRdon")
library("readxl")
library("data.table")
library("FactoMineR")
library("factoextra")
library("plotly")
library("shinythemes")
library("BiocManager")
library("stringr")
options(repos = BiocManager::repositories())


load("./gdsc_data_v3.Robj")
load("./drug_names.Robj")
load("./example_file.Robj")


target_pathway <- gdsc[,c("DRUG_ID", "PUTATIVE_TARGET", "PATHWAY_NAME")]
target_pathway <- unique(target_pathway)
target_pathway$DRUG_NAME <- drug_names$DRUG_NAME[match(target_pathway$DRUG_ID, drug_names$DRUG_ID)]
target_pathway <- target_pathway[,-c(which(colnames(target_pathway)=="DRUG_ID"))]
target_pathway <- unique(target_pathway)


get.drug.ID <- function(data){
  if(length(unique(data$DRUG_ID)) > 1){
    return(FALSE)
  } else {
    drugID <- as.character(unique(data$DRUG_ID))
    return(drugID)
  }
}

check.headers <- function(data){
  headers <- colnames(data)
  if("DRUG_ID" %in% headers  == FALSE){
    return(FALSE)
  } else if("DepMap_ID" %in% headers  == FALSE){
    return(FALSE)
  } else if("IC50" %in% headers  == FALSE){
    return(FALSE)
  } else{
    return(TRUE)
  }
}

delete.na <- function(data, n = 0, col_row) {
  if(col_row == "row"){
    #data <- data[rowSums(is.na(data)) <= n,]
    return(data[rowSums(is.na(data)) <= n,])
  } else if(col_row == "col"){
    #data <- data[,colSums(is.na(data)) <= n]
    return(data[,colSums(is.na(data)) <= n])
  }
}

scale.data <- function(data){
  data <- scale(data)
  return(data)
}


#FUNCTION 1
#User input data formatting
formatData <- function(data, drugName){
  data <- data[c("DRUG_ID", "DepMap_ID", "IC50")]
  data <- as.data.table(data)
  
  gdsc_limited <- gdsc[c("DRUG_ID", "DepMap_ID", "IC50")]
  gdsc_limited <- as.data.table(gdsc_limited)
  
  
  #combine input data to drug list
  drugcompare <- rbind(gdsc_limited, data)
  
  drugcompare$DRUG_ID <- as.factor(drugcompare$DRUG_ID)
  drugcompare$DepMap_ID <- as.factor(drugcompare$DepMap_ID)
  drugcompare$IC50 <- as.numeric(drugcompare$IC50)
  drugcompare <- drugcompare[!duplicated(drugcompare[,c("DRUG_ID","DepMap_ID")]),]
  
  
  #casts file in wide format - include with data formating above
  drugcompare <- dcast.data.table(drugcompare, DepMap_ID ~ DRUG_ID)
  
  #using drug_names pre-loaded file to replace drug_ID with drugnames
  colnames(drugcompare) <- drug_names$DRUG_NAME[pmatch(colnames(drugcompare), drug_names$DRUG_ID, duplicates.ok = TRUE)]
  
  #removes cell line datapoints that don't include drug of interest, include with data formatting above
  drugcompare <- drugcompare[!is.na(drugcompare[[drugName]]),]
  
  drugcompare <- as.data.frame(drugcompare)
  
  return(drugcompare)
}



#FUNCTION 2
#remove rows with all NAs
removeNAs <- function(drugcompare, NArows, NAcols){
  
  drugcompare_NAremoved <- delete.na(drugcompare, NArows, "row")
  drugcompare_NAremoved <- delete.na(drugcompare_NAremoved, NAcols, "col")
  
  #Check if length of data frame - n > 100 datapoints for all rows
  if(nrow(drugcompare_NAremoved) < 100){
    return("insufficientData")
  } 
  
  else{
    rownames(drugcompare_NAremoved) <- drugcompare_NAremoved[,1]
    drugcompare_NAremoved <- drugcompare_NAremoved[,-c(1)]
    drugcompare_NAremoved[is.na(drugcompare_NAremoved)] <- 10000
    drugcompare_NAremoved <- log10(drugcompare_NAremoved)
    return(drugcompare_NAremoved)
  }
}



#FUNCTION 3
#get eigenvalue table
getPCAEigValue <- function(drugcompare){
  
  #get PCA
  drugcompare_pca <- PCA(t(drugcompare), scale.unit = TRUE, ncp = 6, graph = FALSE)
  
  #get eigan values
  drugcompare_pca_eigval <- as.data.frame(get_eigenvalue(drugcompare_pca))
  
  return(drugcompare_pca_eigval)
}


#FUNCTION 4
#Generate pca
generatePCA <- function(drugcompare, input_drugname){
  
  #get PCA
  drugcompare_pca <- PCA(t(drugcompare), scale.unit = TRUE, ncp = 6, graph = FALSE)
  
  #Adjust coordinates for PCS
  drugcompare_pca_coord <- as.data.frame(drugcompare_pca$ind$coord)
  drugcompare_pca_cos2 <- as.data.frame(drugcompare_pca$ind$cos2)
  drugcompare_pca_forPlot <- merge.data.frame(drugcompare_pca_coord, drugcompare_pca_cos2, suffixes = c("_coord","_cos2"), by=0)
  
  #filter on drug of interest
  drugcompare_pca_forPlot$symb_inputDrug <- "Other"
  drugcompare_pca_forPlot$symb_inputDrug[which(drugcompare_pca_forPlot$Row.names == input_drugname)] <- input_drugname
  drugcompare_pca_forPlot$symb_inputDrug <- factor(drugcompare_pca_forPlot$symb_inputDrug, levels = c("Other", input_drugname))
  
  #add plot brush key 
  drugcompare_pca_forPlot$key <- row.names(drugcompare_pca_forPlot)
  
  #some drugs have multiple PCA points and are denoted with a '.1', '.2', etc suffix. 
  #create new column that only includes base drug name without unique suffix
  drugcompare_pca_forPlot$SimpleName <- str_remove(drugcompare_pca_forPlot$Row.names, coll(".1"))
  drugcompare_pca_forPlot$SimpleName <- str_remove(drugcompare_pca_forPlot$SimpleName, coll(".2"))
  drugcompare_pca_forPlot$SimpleName <- str_remove(drugcompare_pca_forPlot$SimpleName, coll(".3"))
  
  #add PCA data (cosine value) for brush data
  drugcompare_pca_forPlot$point_size <- sqrt(drugcompare_pca_forPlot$Dim.1_cos2^2 + drugcompare_pca_forPlot$Dim.2_cos2^2)
  
  #add data from gdsc table (putative target and pathway)
  drugcompare_pca_forPlot$PUTATIVE_TARGET <- target_pathway$PUTATIVE_TARGET[match(drugcompare_pca_forPlot$SimpleName, target_pathway$DRUG_NAME)]
  drugcompare_pca_forPlot$PATHWAY_NAME <- target_pathway$PATHWAY_NAME[match(drugcompare_pca_forPlot$SimpleName, target_pathway$DRUG_NAME)]
  

  
  return(drugcompare_pca_forPlot)
}


#FUNCTION 5
#generate pca plot
generatePlot <- function(drugcompare_pca_forPlot, drugcompare_pca_eigval, xDim, yDim){
  xaxis <- paste0("Dim.", xDim, "_coord")
  yaxis <- paste0("Dim.", yDim, "_coord")
  pcaPlot <- ggplot(drugcompare_pca_forPlot, 
                    aes(x = get(xaxis),
                        y = get(yaxis),
                        size = point_size,
                        label = symb_inputDrug, 
                        fill = symb_inputDrug, 
                        text = paste0(Row.names, "\nCosine2: ", signif(point_size,2)),
                        key = key
                    )) +
    geom_point(shape = 21, color = "black") + 
    theme_bw() +
    ylab(paste0("PCA dimension ", yDim, ": ",round(drugcompare_pca_eigval$variance.percent[yDim], digits=1), "% of variance", sep = "")) +
    xlab(paste0("PCA dimension ", xDim, ": ",round(drugcompare_pca_eigval$variance.percent[xDim], digits=1), "% of variance", sep = "")) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    geom_hline(yintercept = 0, alpha = 0.5) +
    theme(axis.title.x = element_text(size=rel(1.5)),axis.text.x  = element_text(angle=45, vjust=0.5, size=rel(1.25))) +
    theme(axis.title.y = element_text(size=rel(1.5)),axis.text.y  = element_text(angle=45, vjust=0.5, size=rel(1.25)))  +
    scale_size_continuous(name = NULL, range = c(1,3)) +
    scale_fill_manual(values = c("#E7B800", "red"), name = "Drug\n") + 
    theme(legend.text = element_text(size = rel(1)), legend.title = element_text(size = rel(1.2)))
  
  pcaPlot <- ggplotly(pcaPlot,
                      tooltip = "text",
                      #autosize = TRUE,
                      height = 550,
                      width = 750
  )
  
  return(pcaPlot)
}




ui <- 
  fluidPage(theme = shinytheme("sandstone"),
            
            navbarPage(title = "Cell line sensitivity drug comparison",
                       tabPanel("Home",
                              
                                fixedRow(column(4,
                                                fileInput("upload", strong("Upload a cell line sensitivity file for one drug:"),
                                                          accept = c(".xlsx", ".xls", ".csv"),
                                                          placeholder = ".xlsx, .xls, .csv")),
                                         column(5),
                                         column(3,
                                                br(),
                                                downloadButton("downloadExample", "Download an example"))
                                ),
                                br(),
                                fixedRow(column(1),
                                         column(7,
                                                plotlyOutput("pcaPlot")),
                                         column(1),
                                         column(3,
                                                conditionalPanel(condition = "output.upload",
                                                                 wellPanel(
                                                                   tags$b("PLOT OPTIONS"),
                                                                   br(),
                                                                   br(),
                                                                   tags$b("Max number of missing values:"),
                                                                   numericInput("row", "Cell line sensitivity:", 100, min = 10, max = 500, step = 10),
                                                                   numericInput("col", "Drug:", 60, min = 10, max = 500, step = 10),
                                                                   br(),
                                                                   tags$b("PCA dimensions displayed:"),
                                                                   selectInput("xaxis", "X-axis PCA dimension:", choices = c(1:6), selected = 1),
                                                                   selectInput("yaxis", "Y-axis PCA dimension:", choices = c(1:6), selected = 2),
                                                                   br(),
                                                                   radioButtons(inputId = "scaled",
                                                                                label = strong("How would you like to view data?"),
                                                                                choices = c("Unscaled", "Scaled")))
                                                ))),
                              #  br(),
                                br(),
                                fixedRow(column(1),
                                         column(8,
                                                verbatimTextOutput("brush")),
                                         column(2,
                                                uiOutput("downloadButton"))
                                )
                       ),
                       tabPanel("About")
            )
  )




# Define server logic required to draw a histogram
server <- 
  function(input, output) {
    
    output$upload <- reactive({
      #return(!is.null(input$upload))
      #return(!is.null(plot()))
      return(!is.null(drugCompareFile()))
    })
    
    outputOptions(output, "upload", suspendWhenHidden = FALSE)
    
    upload <- reactive({
      inFile <- input$upload
      ext <- tools::file_ext(inFile$datapath)
      if (is.null(inFile))
        return(NULL)
      
      if(ext == "xlsx" | ext == "xls"){
        df <- as.data.frame(read_excel(inFile$datapath, col_names = TRUE))
      } else if(ext == "csv"){
        df <- as.data.frame(read.csv(inFile$datapath, header = TRUE, comment.char = "#"))
      }
      
      if(check.headers(df) == FALSE){
        return("ImproperHeaders")
      } else if(get.drug.ID(df)==FALSE){
        return("MultipleDrugs")
      } else{
        return(df)
      }
    })
    
    drugID <- function(){
      req(upload())
      validate(
        need(upload() != "ImproperHeaders", 
             "Required column names not found. \nPlease upload a file with 3 columns: 'DRUG_ID', 'DepMap_ID', 'IC50'."))
      validate(
        need(upload() != "MultipleDrugs", 
             "Multiple drugs detected in input file. At this time, only 1 drug is supported. \nPlease upload a file with a single DRUG_ID value."))
      get.drug.ID(upload())
    }
    
    drugCompareFile <- reactive(
      formatData(upload(), drugID())
    )
    
    drugCompareFile_NAremoved <- reactive({
      validate(
        need(removeNAs(drugCompareFile(), input$row, input$col) != "insufficientData", "Too few data points. Please increase the maximum number of NAs allowed for sensitivity values."))
      if(input$scaled == "Scaled"){
        return(scale.data(removeNAs(drugCompareFile(), input$row, input$col)))
      } else{
        return(removeNAs(drugCompareFile(), input$row, input$col))
      }
    })
    
    
    
    pcaEigValues <- reactive(
      getPCAEigValue(drugCompareFile_NAremoved())
    )
    
    pca <- reactive({
      req(drugID())
      generatePCA(drugCompareFile_NAremoved(), drugID())
    })
    
    
    selectedPoints <- reactive({
      req(event_data("plotly_selected"))
      brushData <- event_data("plotly_selected")
      brushDataKey <- as.character(brushData$key)
      table <- as.data.frame(pca()[brushDataKey,"Row.names"])
      colnames(table)[1] <- "Selected Points"
      table$Cosine2 <- signif(pca()[brushDataKey,"point_size"],4)
      table$Pathway <- pca()[brushDataKey,"PATHWAY_NAME"]
      table$"Putative target" <- pca()[brushDataKey,"PUTATIVE_TARGET"]
      return(table)
    })
    
    
    plot <- reactive({
      validate(
        need(removeNAs(drugCompareFile(), input$row, input$col) != "insufficientData", ""))
      generatePlot(pca(), pcaEigValues(), as.numeric(input$xaxis), as.numeric(input$yaxis))
    })
    
    
    output$pcaPlot <- renderPlotly({
      req(pca())
      req(pcaEigValues())
      plot()
    })
    
    
    output$brush <- renderPrint({
      validate(
        need(upload() != "ImproperHeaders", ""))
      validate(
        need(upload() != "MultipleDrugs", ""))
      req(plot())
      if(is.null(event_data("plotly_selected")) == FALSE){
        print(selectedPoints(), right = FALSE, width = 200)
      } else{
        cat("Use the Box Select or Lasso Select tool at the upper \nright of the graph to select some data points")
      }
    })
    
    output$downloadExample <- downloadHandler(
      filename = function() {
        paste("CellLine_IC50_Example", ".csv" , sep = "")
      },
      content = function(file) {
        write.table(Example, file, row.names = FALSE, col.names = FALSE, sep=",", quote = FALSE)
      })
    
    output$downloadButton <- renderUI({
      req(selectedPoints())
      downloadButton("downloadSelected", "Download selected data")
    })
    
    output$downloadSelected <- downloadHandler(
      filename = function() {
        req(selectedPoints())
        paste("CellLineSensitivity_SelectedDataPoints", ".csv" , sep = "")
      },
      content = function(file) {
        req(selectedPoints())
        write.csv(selectedPoints(), file, row.names = FALSE)
      })
    
  }


# Run the application 
shinyApp(ui = ui, server = server)
