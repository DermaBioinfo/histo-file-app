library(shiny)
library(ggplot2)
library(reshape2)

# Define UI for application that draws a histogram
ui <- fluidPage(
   # Application title
   titlePanel("Histo Analysis"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
        mainPanel(
            textOutput("result.text", container = tags$h3),
            tabsetPanel(
                tabPanel("Staining Distribution (All)", plotOutput("hist.plot", height="600px")),
                tabPanel("Staining Distribution (Positive)", plotOutput("hist.plot.pos", height="600px")),
                tabPanel("PCA", 
                         fluidRow(
                           column(9, plotOutput("pca")),
                           column(3, plotOutput("location.frequ"))
                         )
                ),
                tabPanel("Phenotypes",
                         fluidRow(
                           column(12, tags$h4("Currently stored phenotypes:")),
                           column(12, tableOutput("phenotype.table")),
                           column(12, tags$hr()),
                           column(12, tags$h4("Add new phenotype:")),

                           column(3,
                                   textOutput("n1.ab1", container = tags$b, inline = T),
                                   selectInput("phen.ab1", label = NULL, choices = c("Ignore", "Up", "Down"))
                           ),
                           column(3,
                                  textOutput("n1.ab2", container = tags$b, inline = T),
                                  selectInput("phen.ab2", label = NULL, choices = c("Ignore", "Up", "Down"))
                           ),
                           column(3,
                                  textOutput("n1.ab3", container = tags$b, inline = T),
                                  selectInput("phen.ab3", label = NULL, choices = c("Ignore", "Up", "Down"))
                           ),
                           column(3,
                                  textOutput("n1.ab4", container = tags$b, inline = T),
                                  selectInput("phen.ab4", label = NULL, choices = c("Ignore", "Up", "Down"))
                           ),
                           column(3,
                                  textOutput("n1.ab5", container = tags$b, inline = T),
                                  selectInput("phen.ab5", label = NULL, choices = c("Ignore", "Up", "Down"))
                           ),
                           column(3,
                                  textOutput("n1.ab6", container = tags$b, inline = T),
                                  selectInput("phen.ab6", label = NULL, choices = c("Ignore", "Up", "Down"))
                           ),
                           column(3,
                                    textOutput("n1.ab7", container = tags$b, inline = T),
                                    selectInput("phen.ab7", label = NULL, choices = c("Ignore", "Up", "Down"))
                           ),
                           column(3, 
                                  textInput("phen.name", "Name", placeholder = "B cells"),
                                  actionButton("save.phenotype", "Save"))
                         )
                )
            )
        ),       
      sidebarPanel(
          fileInput(inputId = "core.file", label = "Input text file", buttonLabel = "Browse"),
          tags$p("When uploading a file, this will overwrite the currently loaded data!"),
          tags$hr(),
          textOutput("n.ab1", container = tags$b),
          sliderInput("t.ab1", NULL, min = 0, max = 20, value = 0, step = 0.1),
          selectInput("s.ab1", NULL, choices = c("Membrane", "Cytoplasm", "Nucleus", "Entire Cell")),
          textOutput("n.ab2", container = tags$b),
          sliderInput("t.ab2", NULL, min = 0, max = 20, value = 0, step = 0.1),
          selectInput("s.ab2", NULL, choices = c("Membrane", "Cytoplasm", "Nucleus", "Entire Cell")),
          textOutput("n.ab3", container = tags$b),
          sliderInput("t.ab3", NULL, min = 0, max = 20, value = 0, step = 0.1),
          selectInput("s.ab3", NULL, choices = c("Membrane", "Cytoplasm", "Nucleus", "Entire Cell")),
          textOutput("n.ab4", container = tags$b),
          sliderInput("t.ab4", NULL, min = 0, max = 20, value = 0, step = 0.1),
          selectInput("s.ab4", NULL, choices = c("Membrane", "Cytoplasm", "Nucleus", "Entire Cell")),
          textOutput("n.ab5", container = tags$b),
          sliderInput("t.ab5", NULL, min = 0, max = 20, value = 0, step = 0.1),
          selectInput("s.ab5", NULL, choices = c("Membrane", "Cytoplasm", "Nucleus", "Entire Cell")),
          textOutput("n.ab6", container = tags$b),
          sliderInput("t.ab6", NULL, min = 0, max = 20, value = 0, step = 0.1),
          selectInput("s.ab6", NULL, choices = c("Membrane", "Cytoplasm", "Nucleus", "Entire Cell")),
          textOutput("n.ab7", container = tags$b),
          sliderInput("t.ab7", NULL, min = 0, max = 20, value = 0, step = 0.1),
          selectInput("s.ab7", NULL, choices = c("Membrane", "Cytoplasm", "Nucleus", "Entire Cell"))
      )
   )
)

fixData <- function(data) {
    if ("Confidence" %in% colnames(data)) {
        data$Confidence <- gsub(" %", "", data$Confidence)
        data$Confidence <- as.numeric(data$Confidence) / 100
    }
    
    return(data)
}

extractRelevantColumns <- function(data) {
    cols.to.keep <- colnames(data) %in% c("Cell.ID", "Sample.Name", "Cell.X.Position", "Cell.Y.Position", "Confidence", "Phenotype")
    cols.to.keep <- cols.to.keep |
        grepl("^Cytoplasm\\..*.*Mean\\.\\.Normalized\\.Counts\\..*", colnames(data)) |
        grepl("^Membrane\\..*.*Mean\\.\\.Normalized\\.Counts\\..*", colnames(data)) |
        grepl("^Nucleus\\..*.*Mean\\.\\.Normalized\\.Counts\\..*", colnames(data)) |
        grepl("^Entire\\.Cell\\..*.*Mean\\.\\.Normalized\\.Counts\\..*", colnames(data)) |
        grepl("Tissue\\.Category", colnames(data))
    
    return(data[, cols.to.keep])
}

shortenColnames <- function(data) {
    colnames(data) <- gsub("([^.]*)\\.(.*)\\.\\.Opal.*", "\\2.\\1", colnames(data))
    colnames(data) <- gsub("([^.]*)\\.DAPI.*", "DAPI.\\1", colnames(data))
    #colnames(data) <- gsub("\\.", "-", colnames(data))
    colnames(data) <- gsub("^Cell\\.", "", colnames(data))
    colnames(data)[grepl("Entire.DAPI.Cell", colnames(data))] <- "DAPI.Entire"
    colnames(data)[colnames(data) == "Phenotype"] <- "Phenotype.org"
    
    return(data)
}

#' Stores the current list of phenotypes
#' 
#' @param phenotypes A list where the names are the phenotype names and vector of antibodies their description.
#' @example 
#' phenotypes <- list()
#' # + / - indicate whether the molecule is up- or down-regulated
#' phenotypes[["B.cell"]] <- c("MS4A1+", "CD19+")
#' savePhenotypes(phenotypes)
savePhenotypes <- function(phenotypes) {
  # currently, only save as an RDS file
  saveRDS(phenotypes, "phenotypes.rds")
}

#' Retrieves the stored phenotypes
#' 
#' @return A named list with the phenotype names and a vector of molecules
loadPhenotypes <- reactive({
  if (file.exists("phenotypes.rds")) {
    return(readRDS("phenotypes.rds"))
  } else {
    return(list())
  }
})

# increase the upload size
options(shiny.maxRequestSize=100*1024^2)

# Define server logic required todraw a histogram
server <- function(input, output) {
  v <- reactiveValues(phenotypes = list())
  
  # check whether a new file was uploaded by the user
  loadCoreData <- reactive({
    if ("core.file" %in% names(input) && !is.null(input$core.file)) {
      message("Loading file..")
      data <- read.csv(input$core.file$datapath, sep = "\t", na.strings = c("NA", "N/A", "#N/A"))
      data <- fixData(data)
      data <- extractRelevantColumns(data)
      data <- shortenColnames(data)
      
      # also update the antibody names
      ab.names <<- gsub("\\.Membrane", "", colnames(data)[grep("\\.Membrane", colnames(data))])
      
      cores <<- data
      
      message(paste0(colnames(cores), collapse = ", "))
    }
    
    if (exists("cores")) {
      return(cores)
    } else {
      return(NULL)
    }
  })
  
  # update the labels
  getAbLabel <- function(n.ab) {
    cores <- loadCoreData()
    
    if (exists("ab.names")) {
      return(ab.names[n.ab])
    } else {
      return("")
    }
  }
  
  output$n.ab1 <- renderText(getAbLabel(1))
  output$n.ab2 <- renderText(getAbLabel(2))
  output$n.ab3 <- renderText(getAbLabel(3))
  output$n.ab4 <- renderText(getAbLabel(4))
  output$n.ab5 <- renderText(getAbLabel(5))
  output$n.ab6 <- renderText(getAbLabel(6))
  output$n.ab7 <- renderText(getAbLabel(7))
  
  output$n1.ab1 <- renderText(getAbLabel(1))
  output$n1.ab2 <- renderText(getAbLabel(2))
  output$n1.ab3 <- renderText(getAbLabel(3))
  output$n1.ab4 <- renderText(getAbLabel(4))
  output$n1.ab5 <- renderText(getAbLabel(5))
  output$n1.ab6 <- renderText(getAbLabel(6))
  output$n1.ab7 <- renderText(getAbLabel(7))
  
  # build the colnames
  getColnames <- reactive({
    ab.names <- getAbLabel()
    req(ab.names)
    
    ab.colnames <- c()
    for (n in 1:7) {
      staining.name <- input[[paste0("s.ab", n)]]
      
      if (staining.name == "Entire Cell") {
        staining.name <- "Entire"
      }
      
      ab.colnames <- c(ab.colnames, paste0(ab.names[n], ".", staining.name))
    }
    
    return(ab.colnames)
  })
  
  # ---- apply the thresholds ----
  getPositiveCells <- reactive({
    cores <- loadCoreData()
    ab.colnames <- getColnames()
    req(cores, ab.colnames)
    
    is.positive <- rep(T, nrow(cores))
  
    for (n in 1:7) {
      ab.col <- ab.colnames[n]
      
      if (!ab.col %in% colnames(cores)) {
        stop("Error: No values for ", ab.col)
      }
      
      ab.threshold <- input[[paste0("t.ab", n)]]
      
      is.positive <- is.positive & !(is.na(cores[, ab.col])) & cores[, ab.col] >= ab.threshold
    }
    
    return(is.positive)
  })
  
  # ---- Create the outputs ----
  # add the frequency table
  output$result.text <- renderText({
    cores <- loadCoreData()
    req(cores)
    is.positive <- getPositiveCells()
    
    paste0(sum(is.positive, na.rm = T), " / ", nrow(cores), " (",
           round(sum(is.positive, na.rm = T) / nrow(cores) * 100, 1), "%) positive cells")
  })
  
  # render the histograms
  output$hist.plot <- renderPlot({
    cores <- loadCoreData()
    ab.colnames <- getColnames()
    req(cores, ab.colnames)
    
    par(mfrow=c(3,3))
    # apply the thresholds
    for (n in 1:7) {
      threshold <- input[[paste0("t.ab", n)]]
      hist(cores[, ab.colnames[n]], main = ab.colnames[n], breaks = 50)
      abline(v = threshold, col = "red", lty = 2)
    }
  })
  
  output$hist.plot.pos <- renderPlot({
    cores <- loadCoreData()
    ab.colnames <- getColnames()
    req(cores, ab.colnames)
    
    # create the histograms
    par(mfrow=c(3,3))
    for (n in 1:7) {
      # get the abs threshold
      ab.threshold <- input[[paste0("t.ab", n)]]
      
      # only plot positive cells
      is.positive <- cores[, ab.colnames[n]] >= ab.threshold
      
      # create the histogram
      hist(cores[is.positive, ab.colnames[n]], main = ab.colnames[n], breaks = 50)
    }
  })
  
  # location.frequ
  output$location.frequ <- renderPlot({
    cores <- loadCoreData()
    is.positive <- getPositiveCells()
    req(cores, is.positive)
    
    if (!"Tissue.Category" %in% colnames(cores)) {
      return(NULL)
    }
    
    ggplot(cores, aes(x = Tissue.Category)) +
      geom_bar() +
      theme_bw() +
      labs(title = "Positive cells")
  })
  
  # perform the PCA
  output$pca <- renderPlot({
    cores <- loadCoreData()
    is.positive <- getPositiveCells()
    req(cores, is.positive)
    
    staining.cols <- grepl("Membrane", colnames(cores)) |
      grepl("Nucleus", colnames(cores)) |
      grepl("Cytoplasm", colnames(cores))
    
    pca.exprs <- t(cores[, staining.cols])
    pca.exprs[is.na(pca.exprs)] <- 0
    
    fit <- prcomp(pca.exprs)
    plot.data <- data.frame(fit$rotation)
    plot.data$is.positive <- is.positive
    
    ggplot(plot.data, aes(x = PC1, y = PC2, color = is.positive)) +
      geom_point() +
      scale_color_manual(values = c("grey", "red"))
  })
  
  # ---- Phenotypes ----
  # create the table of currently stored phenotypes
  output$phenotype.table <- renderTable({
    v$phenotypes <- loadPhenotypes()
    
    if (length(v$phenotypes) < 0) {
      return(NULL)
    }
    
    # create the data.frame
    phen.data <- do.call("rbind", lapply(names(v$phenotypes), function(name) {
      data.frame(
        name = name,
        molecules = paste0(v$phenotypes[[name]], collapse = ", ")
      )
    }
    ))
    
    return(phen.data)
  })
  
  # observe the button to add a new phenotype
  observeEvent(input$save.phenotype, {
    req(input$phen.name)
    ab.names <- getColnames()
    
    # build the phenotype string
    phenotype.vector <- c()
    
    for (n in 1:7) {
      # get the selection
      sel <- input[[paste0("phen.ab", n)]]
      
      if (sel == "Up") {
        phenotype.vector <- c(phenotype.vector, paste0(ab.names[n], "+"))
      } else if (sel == "Down") {
        phenotype.vector <- c(phenotype.vector, paste0(ab.names[n], "-"))
      }
    }
    
    if (length(phenotype.vector) > 0) {
      message("Saving phenotypes")
      
      v$phenotypes <- loadPhenotypes()
      v$phenotypes[[input$phen.name]] <- phenotype.vector
      savePhenotypes(v$phenotypes)
    }
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)