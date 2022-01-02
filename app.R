#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
## For Shiny and plots
library(shiny)
library(shinythemes)
library(DT)
library(plotly)


## Read in functions
library(readxl)
source("./direct_functions.R", local = TRUE)


##### Read in the original data file
mydata <- prepare.data("data/piano.csv")


ui <- navbarPage(
        title = "Mass Spectrometry", 
        theme = shinytheme("cerulean"),
        
        tabPanel("Setup",
                 

                 
                 fluidRow(
                   column( 
                   textInput('element.1','Element 1','C'),
                   width = 3),        
                   column(
                     numericInput("ele.1.mass", "Mass", 12, 0, 1e10,1e-4) ,
                     width = 3,
                          offset = 0),
                   column(
                     numericInput("ele.1.max", "Maximum number", 80, 0, 80,1) ,
                     width = 4
                   )
                   ),
                 fluidRow(
                   column( 
                     textInput('element.2','Element 2','H'),
                     width = 3),        
                   column(
                     numericInput("ele.2.mass", "Mass", 1.00783, 0, 1e10,1e-4) ,
                     width = 3,
                     offset = 0),
                   column(
                     numericInput("ele.2.max", "Maximum number", 80, 0, 80,1) ,
                     
                     width = 4
                   )
                 ),

                 fluidRow(
                   column( 
                     textInput('element.3','Element 3','O'),
                     width = 3),        
                   column(
                     numericInput("ele.3.mass", "Mass", 15.99491, 0, 1e10,1e-4) ,
                     width = 3,
                     offset = 0),
                   column(
                     numericInput("ele.3.max", "Maximum number", 10, 0, 80,1) ,
                     
                     width = 4
                   )
                 ),
                 fluidRow(
                   column( 
                     textInput('element.4','Element 4',''),
                     width = 3),        
                   column(
                     numericInput("ele.4.mass", "Mass", 0, 0, 1e10,1e-4) ,
                     width = 3,
                     offset = 0),
                   column(
                     numericInput("ele.4.max", "Maximum number", 1, 0, 80,1) ,
                     
                     width = 4
                   )
                 ),
                
                 
                 fluidRow(
                   column( 
                     textInput('element.5','Element 5',''),
                     width = 3),        
                   column(
                     numericInput("ele.5.mass", "Mass", 0, 0, 1e10,1e-4) ,
                     width = 3,
                     offset = 0),
                   column(
                     numericInput("ele.5.max", "Maximum number", 1, 0, 80,1) ,
                     
                     width = 4
                   )
                 ),
                 fluidRow(
                   column(
                     fileInput("file1", "Upload your CSV file",
                               multiple = TRUE,
                               accept = c(
                                 ".csv"
                               )),
                     width=3
                   ),
                   column(
                     numericInput("initial.offest", "Offest", 1.0071074, 0, 1e10,1e-4) ,
                     width = 2,
                     offset = 0),
                   
                   column(
                     numericInput("iteration", "Iterations", 3, 1, 1e2,1) ,
                     width = 2,
                     offset = 0),
            
                   column(
                     numericInput("ppm", "ppm threshold", 10, 0, 1e2,1) ,
                     width = 2,
                     offset = 0)
                 ),
                 fluidRow(
                 column(
                 downloadButton("downloadData_ex", "Download an example file"),
                 width=4,align= "center"),
                 column(
                   textOutput("text_memory"), 
                   width=6
                 )
                 
                 )

        ), ## end of tab 1

        tabPanel("Result",
                 
                 fluidPage(
                   fluidRow(
                   column(
                     "Please hit the button after a table appears below.",
                     br(),
                     downloadButton("downloadData", "Download Results"),
                     align = "center", width = 8, offset = 2
                   )
                   ),
                   hr(),
                 fluidRow(
                   column(
                   dataTableOutput("user_table"),
                   width=12)
                 )
                 )
                 
  
      
        ), ## end of tab 2

        tabPanel("Plot",
                 fluidPage(

                   plotlyOutput("scatterPlot"),
                   
                 )
        ) ## end of tab 3

    ) # end of navbarPage






server <- function(input, output) {
    

    inputVar <- reactive({
        
        
        ele.names <- c(input$element.1, input$element.2, input$element.3, input$element.4, input$element.5)
        ele.masses <- c(input$ele.1.mass, input$ele.2.mass, input$ele.3.mass, input$ele.4.mass, input$ele.5.mass)
        ele.max <- c(input$ele.1.max, input$ele.2.max, input$ele.3.max, input$ele.4.max, input$ele.5.max)
        non.na.loc <- which(ele.names != '')
        
        value <- list(ele.names[non.na.loc], 
                      ele.masses[non.na.loc], 
                      ele.max[non.na.loc], input$initial.offest, input$iteration,
                      input$ppm)

    })
    
    
    output$text_memory <- renderText({
      value <- inputVar()
      
      file <- input$file1
      ext <- tools::file_ext(file$datapath)
      # validate(need(ext == "csv", "Please upload a csv file"))
      
      req(input$file1)
      
      user_data <- read.csv(input$file1$datapath,
                            header = T,
                            sep = ",",
      )
      feature.ID<-user_data[,1]
      m.z<-user_data[,2]
      if (ncol(user_data)>2){second.dim<-user_data[,3]}else{second.dim<-NULL}
      n_user_data <- dim(user_data)[1]
      
      
      
      block.name<-value[[1]]
      block.mass<-value[[2]]
      block.max<-value[[3]]
      
      
      offset<- value[[4]] # positive ion mode
      ppm.thresh<-value[[6]]
      
      ### Approximate memory allocation
      num_cell_1 <- prod(block.max)
      num_cell_2 <- num_cell_1 * n_user_data
      
      memo_1 <- 400 + 8 * num_cell_1
      memo_2 <- 400 + 8 * num_cell_2
      
      memo <- sum(memo_1,memo_2) / 1e+9
      
      # paste("If ", memo, " is greater than 1, it will take more time to get the results.", sep="")
      # if( memo > 1e+9){
      #   paste("The largest memory usage of your setting will be ", round(memo,2), "GB (>1GB) and it will take more time to get the output.", sep="")
      # } else {
      #   paste("The largest memory usage of your setting will be ", round(memo,2), "GB (<1GB) and it will not take much time to get the output.", sep="")
      #   
      # }
      paste("The largest memory usage of your setting will be ", round(memo,2), "GB. If it is greater than 1GB, it will take more time to get the output.", sep="")
      
    })
    
    output$user_table <- renderDataTable({
        result_data <- data_down_reactive()
        datatable(result_data)
    })
    
    
    
    
    # Downloadable csv of selected dataset ----
    
    data_down_reactive <- reactive({
      value <- inputVar()
      
      file <- input$file1
      ext <- tools::file_ext(file$datapath)
      validate(need(ext == "csv", "Please upload a csv file"))
      
      req(input$file1)
      
      user_data <- read.csv(input$file1$datapath,
                            header = T,
                            sep = ",",
                            )
      feature.ID<-user_data[,1]
      m.z<-user_data[,2]
      if (ncol(user_data)>2){second.dim<-user_data[,3]}else{second.dim<-NULL}
      n_user_data <- dim(user_data)[1]
      
      
      
      block.name<-value[[1]]
      block.mass<-value[[2]]
      block.max<-value[[3]]
      

      offset<- value[[4]] # positive ion mode
      ppm.thresh<-value[[6]]
      
      ### Approximate memory allocation
      num_cell_1 <- prod(block.max)
      num_cell_2 <- num_cell_1 * n_user_data
      
      memo_1 <- 400 + 8 * num_cell_1
      memo_2 <- 400 + 8 * num_cell_2
      
      if( sum(memo_1, memo_2) < 1e+9){
        temp<-findbest(feature.ID,m.z,offset,ppm.thresh,block.name,block.mass,block.max,second.dim,numiter=value[[5]])
        
      } else {
        temp<-findbest2(feature.ID,m.z,offset,ppm.thresh,block.name,block.mass,block.max,second.dim,numiter=value[[5]])
        
      }

      temp      
    })
    
    output$downloadData <- downloadHandler(
      
      filename = function() {
        "result_file.csv"
      },
      content = function(file) {
        write.csv(data_down_reactive(), file, row.names = FALSE)
      }
    )
    
    
    output$downloadData_ex <- downloadHandler(
      
      filename = "example_data.csv",
      content = function(file) {
        file.copy('data/piano.csv', file)
      }
    )
    
    output$scatterPlot <- renderPlotly({
      finaldf <- data_down_reactive()

      # plot_ly(data = finaldf, x = ~m.z.corrected, y = ~second.dim,text = ~paste("ID:",feature.ID,
      #                                                                                "<br>predicted formula:",our.formula,
      #                                                                                "<br>obs. corrected mass:",round(m.z.corrected,6),
      #                                                                                "<br>predicted mass:",round(masses.predicted,6),
      #                                                                                "<br>ppm diff:",ppm.diff),
      #              color=~mycolornum,colors=c("red","black","blue"))
      plot_ly(data = finaldf, x = ~m.z.corrected, y = ~second.dim,text = ~paste("ID:",feature.ID,
                                                                                     "<br>predicted formula:",our.formula,
                                                                                     "<br>obs. corrected mass:",round(m.z.corrected,6),
                                                                                     "<br>predicted mass:",round(masses.predicted,6),
                                                                                     "<br>ppm diff:",ppm.diff),
                   color=~mycolornum)
    })
    

    
    
} # end of server

# Run the application 
shinyApp(ui = ui, server = server)
