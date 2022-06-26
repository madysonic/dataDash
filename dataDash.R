# # # # # # # # # # # # # # # # # # # # # # # # #
#              ALL Data Dashboard               #
# # # # # # # # # # # # # # # # # # # # # # # # #

# Set Up ------------------------------------------------------------------

### Load Libraries ###

library(lessR)
library(shinyWidgets)
library(shinythemes)
library(shiny)
library(shinydashboard)
library(plotly)
library(dplyr)
library(ggplot2)
library(RSQLite)
library(edgeR)
library(magrittr)
library(Rtsne)
library(reshape2)
library(tibble)

### Source Functions ###
source("dashFunctions.R")

# Establish a read-only connection
LBank <- dbConnect(RSQLite::SQLite(), "Leukaemia_Bank", flags = SQLITE_RO)

### Get DGE Object for Wrangling ###
dge <- getLBankDGE()
saveRDS(dge, "LBankDGE.rds")
#dge <- readRDS("LBankDGE.rds")


# UI ----------------------------------------------------------------------


ui <- navbarPage(title = "ALL Dashboard",
                 
                 ### Overview ### 
                 tabPanel(title = "Cohort Overview", 
                          fluidPage(theme = shinytheme("yeti"),
                                    fluidRow(#box(title = "Filter by subtype:",
                                      pickerInput("Subtype",
                                                  label = "Choose subtype/s:",
                                                  choices = unique(sort(dge$samples$col)),
                                                  multiple = T, 
                                                  options= list(`actions-box` = TRUE),
                                                  width = "300",
                                                  inline = T,
                                                  selected = unique(sort(dge$samples$col))
                                      ),
                                      #fluidRow(box(title = "Filter by age group:",
                                      pickerInput("pickAge",
                                                  label = "Choose age group/s:",
                                                  choices = unique(sort(dge$samples$ageGp)),
                                                  multiple = T, 
                                                  options= list(`actions-box` = TRUE),
                                                  width = "300",
                                                  inline = T,
                                                  selected = unique(sort(dge$samples$ageGp))
                                      ),
                                      
                                      #fluidRow(box(title = "Filter by timepoint:",
                                      pickerInput("pickTime",
                                                  label = "Choose timepoint:",
                                                  choices = unique(sort(dge$samples$sample)),
                                                  multiple = T, 
                                                  options= list(`actions-box` = TRUE),
                                                  width = "300",
                                                  inline = T,
                                                  selected = unique(sort(dge$samples$sample))
                                      ),
                                      pickerInput("pickType",
                                                  label = "Choose sample type:",
                                                  choices = unique(sort(dge$samples$sampleType)),
                                                  multiple = T, 
                                                  options= list(`actions-box` = TRUE),
                                                  width = "300",
                                                  inline = T,
                                                  selected = unique(sort(dge$samples$sampleType))
                                      ),
                                      
                                      fluidRow(title = "Cohort Subtypes", 
                                               plotOutput("colPie", height = 500, width = 600),
                                               title = "Cohort Age Groups",
                                               plotOutput("agePie", height = 500, width = 600),
                                               title = "Cohort Timepoints",
                                               plotOutput("timePie", height = 500, width = 600),
                                               title = "Cohort Sample Type",
                                               plotOutput("sampPie", height = 500, width = 600)),
                                      
                                    ))),
                 
                 
                 ### Download Data ###      
                 tabPanel(title = "Browse and Download Data",
                          fluidPage(theme = shinytheme("yeti"),
                                    fluidRow(h3("Select Data:"),
                                             pickerInput("subSelect",
                                                         label = "Choose primary subtype/s:",
                                                         choices = unique(sort(dge$samples$col)),
                                                         multiple = T, 
                                                         options= list(`actions-box` = TRUE),
                                                         width = "300",
                                                         inline = T,
                                                         selected = "Hyperdiploid"),
                                             
                                             #fluidRow(box(title = "Filter by age group:",
                                             pickerInput("subAgeSelect",
                                                         label = "Choose age group/s:",
                                                         choices = unique(sort(dge$samples$ageGp)),
                                                         multiple = T, 
                                                         options= list(`actions-box` = TRUE),
                                                         width = "300",
                                                         inline = T,
                                                         selected = unique(sort(dge$samples$ageGp))
                                             ),
                                             
                                             #fluidRow(box(title = "Filter by timepoint:",
                                             pickerInput("subTimeSelect",
                                                         label = "Choose timepoint:",
                                                         choices = unique(sort(dge$samples$sample)),
                                                         multiple = T, 
                                                         options= list(`actions-box` = TRUE),
                                                         width = "300",
                                                         inline = T,
                                                         selected = "Dx")),
                                    
                                    fluidRow(dataTableOutput("subTable")),
                                    
                                    fluidRow(h4("Download Table, Counts or SNVs for Selected Samples:"),
                                             h5(""),
                                             downloadButton("downloadMeta", "Download Table"),
                                             downloadButton("downloadSubCounts", "Download Gene Counts"),
                                             downloadButton("downloadSubSNVs", "Download SNVs")),
                                    fluidRow(h1("")),
                                    fluidRow(h3("Or Search For Patient:"),
                                             selectInput("patSelect", 
                                                         label = "Select Patients", 
                                                         choice = unique(sort(dge$samples$LB_ID)),
                                                         multiple = T,
                                                         selected = "ADI_002")),
                                    
                                    fluidRow(dataTableOutput("patTable")),
                                    
                                    fluidRow(h4("Download Table, Counts or SNVs for Selected Patients:"),
                                             h5(""),
                                             downloadButton("downloadPats", "Download Table"),
                                             downloadButton("downloadPatCounts", "Download Gene Counts"),
                                             downloadButton("downloadPatSNVs", "Download SNVs")),
                                    
                                    
                          )),
                 ### TSNE ###     
                 tabPanel(title = "TSNE",
                          fluidPage(theme = shinytheme("yeti"),
                                    fluidRow(title = "tSNE", 
                                             plotlyOutput("tsne", height = 800, width=1000)))),
                 
                 ### Explore Genes ###  
                 tabPanel(title = "Explore Genes",
                          fluidPage(theme = shinytheme("yeti"),
                                    fluidRow(textInput("geneNames", 
                                                       label = "Enter a gene name",
                                                       value = "TP53"),
                                             
                                             actionButton(inputId = "plotClickGene", label = "Update plot")),
                                    
                                    fluidRow(plotOutput("plotGene", height = 500))
                          )),
                 
                 ### Explore SNVs ###  
                 tabPanel(title = "Explore SNVs",
                          fluidPage(theme = shinytheme("yeti"),
                                    fluidRow(textInput("snvNames", 
                                                       label = "Enter a gene name",
                                                       value = "KRAS"),
                                             pickerInput("snvSubtype",
                                                         label = "Choose subtype/s:",
                                                         choices = unique(sort(dge$samples$col)),
                                                         multiple = T, 
                                                         options= list(`actions-box` = TRUE),
                                                         width = "fit",
                                                         inline = T,
                                                         selected = "Hyperdiploid"
                                             ),
                                             
                                             actionButton(inputId = "plotClickSNVs", label = "Update plot")),
                                    
                                    fluidRow(plotOutput("percPlotSNVs"))
                          ))
                 
                 
)



# Server ------------------------------------------------------------------


server <- function(input, output, session) {
  
  ### Overview ###  
  output$colPie <- renderPlot({
    PieChart(col,
             rows =(col %in% input$Subtype & ageGp %in% input$pickAge & sample %in% input$pickTime & sampleType %in% input$pickType),
             data = dge$samples,
             #data = dbGetQuery(LBank,"SELECT * FROM Sample_Anno WHERE [exclude] != 'Yes';"),
             hole = 0.5,
             values = "input",
             values_position = "out",
             values_color = c("black"),
             labels_cex = 0.9,
             main = NULL,
             quiet = T)
  })
  
  output$agePie <- renderPlot({
    PieChart(ageGp,
             rows = (col %in% input$Subtype & ageGp %in% input$pickAge & sample %in% input$pickTime & sampleType %in% input$pickType),
             data = dge$samples,
             #data = dbGetQuery(LBank,"SELECT * FROM Sample_Anno WHERE [exclude] != 'Yes';"),
             hole = 0.5,
             values = "input",
             values_position = "out",
             values_color = c("black"),
             labels_cex = 0.9,
             main = NULL,
             quiet = T, )
    
  })
  
  output$timePie <- renderPlot({
    PieChart(sample,
             rows = (col %in% input$Subtype & ageGp %in% input$pickAge & sample %in% input$pickTime & sampleType %in% input$pickType),
             data = dge$samples,
             #data = dbGetQuery(LBank,"SELECT * FROM Sample_Anno WHERE [exclude] != 'Yes';"),
             hole = 0.5,
             values = "input",
             values_position = "out",
             values_color = c("black"),
             labels_cex = 0.9,
             main = NULL,
             quiet = T, )
    
  })
  
  output$sampPie <- renderPlot({
    PieChart(sampleType,
             rows = (col %in% input$Subtype & ageGp %in% input$pickAge & sample %in% input$pickTime & sampleType %in% input$pickType),
             data = dge$samples,
             #data = dbGetQuery(LBank,"SELECT * FROM Sample_Anno WHERE [exclude] != 'Yes';"),
             hole = 0.5,
             values = "input",
             values_position = "out",
             values_color = c("black"),
             labels_cex = 0.9,
             main = NULL,
             quiet = T, )
    
  })
  
  ### Browse and Download ###   
  
  # Download buttons: by patient ID
  
  # Patient Meta
  output$downloadPats <- downloadHandler(
    filename = function() {
      paste0("patientList_Meta.csv")
    },
    content = function(file) {
      patMeta <- reactive({dge$samples %>% filter(LB_ID %in% input$patSelect)
      })
      write.table(patMeta(), file, sep = ",", col.names = TRUE, quote = FALSE)
    }
  )
  
  # Patient Counts
  output$downloadPatCounts <- downloadHandler(
    filename = function() {
      paste0("patientList_Counts.csv")
    },
    content = function(file) {
      patCounts <- reactive({dge$counts[,colnames(dge$counts) %in% input$patSelect]
        
      })
      write.table(patCounts(), file, sep = ",", col.names = TRUE, quote = FALSE)
    }
  )
  
  # Patient SNVs
  output$downloadPatSNVs <- downloadHandler(
    filename = function() {
      paste0("patientList_SNVs.csv")
    },
    content = function(file) {
      
      patSamples <- reactive({dge$samples %>% filter(col %in% input$patSelect) %>% pull(variant_file_ID)%>% 
          pull(variant_file_ID) %>% 
          gsub("^", "'", .) %>% 
          gsub("$", "'", .) %>%
          paste(., collapse = ", ")})
      
      patSNVs <- reactive({dbGetQuery(LBank, paste0("SELECT * FROM Filtered_SNVs WHERE [Sample] IN (", patSamples(), ");")
      )})
      write.table(patSNVs(), file, sep = ",", col.names = TRUE, quote = FALSE)
    }
  )
  
  
  # Download buttons: by subtype   
  
  # Subtype Meta
  output$downloadMeta <- downloadHandler(
    filename = function() {
      paste0("subtype_Meta.csv")
    },
    content = function(file) {
      subMeta <- reactive({dge$samples %>% filter(col %in% input$subSelect)
      })
      write.table(subMeta(), file, sep = ",", col.names = TRUE, quote = FALSE)
    }
  )
  # Subtype Counts
  output$downloadSubCounts <- downloadHandler(
    filename = function() {
      paste0("subtype_Counts.csv")
    },
    content = function(file) {
      subSamples <- reactive({dge$samples %>% filter(col %in% input$subSelect) %>% pull(LB_ID)
      })
      subCounts <- dge$counts[,colnames(dge$counts) %in% subSamples()]
      
      write.table(subCounts, file, sep = ",", col.names = TRUE, quote = FALSE)
    }
  )
  
  # Subtype SNVs
  output$downloadSubSNVs <- downloadHandler(
    filename = function() {
      paste0("subtype_SNVs.csv")
    },
    content = function(file) {
      
      subSamples <- reactive({dge$samples %>% filter(col %in% input$subSelect) %>% pull(variant_file_ID)%>% 
          pull(variant_file_ID) %>% 
          gsub("^", "'", .) %>% 
          gsub("$", "'", .) %>%
          paste(., collapse = ", ")})
      
      subSNVs <- reactive({dbGetQuery(LBank, paste0("SELECT * FROM Filtered_SNVs WHERE [Sample] IN (", subSamples(), ");")
      )})
      write.table(subSNVs(), file, sep = ",", col.names = TRUE, quote = FALSE)
    }
  )
  
  # Table Output 
  
  output$patTable <- shiny::renderDataTable(dge$samples %>% 
                                              dplyr::select(patientID,LB_ID,sample,sex,ageGp,sampleType,key_alt,
                                                            primary.subtype,col,exclude,reason_excl) %>% 
                                              filter(LB_ID %in% input$patSelect) %>%
                                              set_colnames(c("Patient ID", "Sample ID", "Timepoint", "Sex","Age Group", "Sample Type", "Key Alteration", "Primary Subtype", "Subtype", "Exclude", "Reason")),
                                            options = list(
                                              searching = F,
                                              pageLength = 10,
                                              lengthChange = F))
  
  output$subTable <- shiny::renderDataTable(test<-dge$samples %>% 
                                              dplyr::select(patientID,LB_ID,sample,sex,ageGp,sampleType,key_alt,
                                                            primary.subtype,col,exclude,reason_excl) %>% 
                                              filter(col %in% input$subSelect & ageGp %in% input$subAgeSelect & sample %in% input$subTimeSelect) %>%
                                              set_colnames(c("Patient ID", "Sample ID", "Timepoint", "Sex","Age Group", "Sample Type", "Key Alteration", "Primary Subtype", "Subtype", "Exclude", "Reason")),
                                            options = list(
                                              searching = F,
                                              pageLength = 10,
                                              lengthChange = F))
  
  
  
  ### TSNE ###  
  # Function for tSNE data
  tsne_plot <- plotTSNE()
  
  # Set Jacqui's colour palette
  palette_tsne <- c("BCL2-MYC" = "#9E0142", 
                    "Ph" = "#E31A1C", 
                    "IL1B" = "black", 
                    "DUX4" = "#B79F00", 
                    "Ph-like" = "#33A02C", 
                    "iAMP21" = "#619CFF",
                    "PAX5.P80R" = "#6A3D9A",
                    # "HLF" = "#c2b280",
                    "HLF" = "#9BB306",
                    "IKZF1.N159Y" = "#1F78B4",
                    "KMT2A" = "#00BFC4",
                    "CRLF2" = "#FDBF6F", 
                    "PAX5alt" = "#A6CEE3",
                    "ETV6-RUNX1" = "#AE87FF",
                    "T-ALL" = "#08306B",
                    "TCF3-PBX1" = "#F0027F", 
                    "ZNF384" = "#FF61C3",
                    "Hyperdiploid" = "#276419",
                    "Hypodiploid" = "#F8766D",
                    "iAMP21" = "#8db600",
                    "MEF2D" = "#f38400",
                    "MLLT10" = "#654522")
  
  tsnePlot <- ggplot(tsne_plot, 
                     aes(x,y,color=subtype, 
                         label = LB_ID,
                         text = paste("patientID: ", patientID, "\n",
                                      "key_alt: ", key_alt, sep = ""))) + 
    geom_point(size = 3) + 
    theme_bw() + 
    labs(title = "tSNE: Good quality diagnostic samples - top 2000 variable genes",        
         colour = "Key alteration") + 
    scale_colour_manual(values = palette_tsne, na.value = "grey50") + 
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank())
  
  
  output$tsne <- renderPlotly({ggplotly(tsnePlot)
  })
  
  
  ### Genes ###
  
  # Get Ensembl IDs
  geneEnsembl <- reactive({dge$genes %>%
      filter(external_gene_name %in% input$geneNames) %>%
      pull(ensembl_gene_id)})
  
  # geneEnsembl <- reactive({dbGetQuery(LBank, 
  #                                     paste0("SELECT [ensembl_gene_id] FROM Gene_Anno WHERE [external_gene_id] = '", input$geneNames,"';"))
  # })
  
  # Extract Counts
  geneData <- reactive({removeBatchEffect(cpm(dge, log = TRUE), 
                                          batch=as.factor(c(dge$samples$ref)), 
                                          batch2 = as.factor(c(dge$samples$sex))) %>%
      as.data.frame %>%
      rownames_to_column("ensembl_gene_id") %>%
      dplyr::filter(ensembl_gene_id %in% geneEnsembl()) %>%
      melt() %>%
      set_colnames(c("ensembl_gene_id", "LB_ID", "logCPM")) %>%
      left_join(dge$samples, by = "LB_ID")  %>%
      left_join(dge$genes, by = "ensembl_gene_id")})
  
  # Reactive Plot Button
  geneOut <- eventReactive(input$plotClick, {ggplot(geneData(), aes(x=col, y=logCPM, fill=col)) +  
      geom_point(aes(color=col), position = position_jitterdodge()) +
      theme_classic() +
      theme(legend.position = "none", axis.text.x = element_text(angle=45, vjust = 1, hjust= 1, size=9))+
      labs(x="Subtype") +
      #facet_wrap(.~external_gene_name, scales='free_y') +
      ggtitle("Expression of selected gene across subtypes")})
  
  output$plotGene <- renderPlot({geneOut()
  })
  
  
  ### Explore SNVs ###
  
  groupSamples <- reactive({dge$samples %>% 
      filter(col %in% input$snvSubtype) %>%
      as.data.frame()})
  
  varIDs <- reactive({groupSamples() %>% 
      pull(variant_file_ID) %>% 
      gsub("^", "'", .) %>% 
      gsub("$", "'", .) %>%
      paste(., collapse = ", ")})
  
  # Filter SNV data for only these samples
  sampleSNVs <- reactive({dbGetQuery(LBank, paste0("SELECT [Gene.refGene], [Sample] FROM Filtered_SNVs WHERE [Sample] IN (", varIDs(), ");"))})
  
  # Filter for gene of interest
  varData <- reactive({sampleSNVs() %>% dplyr::filter(Gene.refGene %in% input$snvNames) %>%
      dplyr::right_join(groupSamples(), by = c("Sample" = "variant_file_ID")) %>%
      mutate(Gene.refGene = if_else(is.na(Gene.refGene), "None", Gene.refGene)) %>%
      dplyr::count(Gene.refGene, primary.subtype) %>%
      mutate(Total = sum(n)) %>%
      dplyr::group_by(primary.subtype) %>%
      mutate(Percentage = n/sum(n)) %>% 
      mutate(varTotal = sum(n)) %>%
      as.data.frame()})
  
  varPal <- c("#E69F00", "#0072B2")
  groupFreq <- reactive({varData()[varData()$Gene.refGene %in% input$snvNames,] %>% mutate(Freq = (sum(n)/Total)*100) %>% pull(Freq)})
  
  percPlot <- eventReactive(input$plotClickSNVs, {ggplot(varData(), aes(x=primary.subtype, y=Percentage, fill=Gene.refGene)) +
      geom_col(position="stack") +
      geom_text(aes(label=paste(paste0(Gene.refGene, " = ", round(Percentage*100, digits=0), "% (", n, ")"))), 
                position=position_stack(vjust = 0.5),
                colour="white",
                size=2.5,
                fontface="bold") +
      theme_classic() + 
      scale_fill_manual(values = varPal) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      labs(x = "Primary Subtype",
           y = "Percentage of Samples",
           fill = "Mutation") +
      ggtitle(paste0("Overall Group Frequency = ", groupFreq(), "%"))})
  
  output$percPlotSNVs <- renderPlot({percPlot()})
  
  
  
  #  session$onSessionEnded(function() {
  #   DBI::dbDisconnect(LBank)
  # }
}



# Launch ------------------------------------------------------------------


shinyApp(ui, server)

DBI::dbDisconnect(LBank)
