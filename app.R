library(shiny)
library(eulerr)
library(scales)


pathway_name <- 'Pathway of Interest'

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Calculating enrichment"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(width = 5, 
         
         h4("Genes"),
         
         
         fluidRow(
            column(width = 9, 
                   sliderInput("n_test_genes", "Differential genes:", value = 300, min = 1, max = 5000, round=TRUE)
                   ),
            column(width = 3,
               fluidRow(numericInput(inputId = "n_test_genes_box", "", value=NA, min=1,  max=5000)),
               fluidRow(align = "centre",  actionButton("setButton_n_test_genes", "Set")) 
               )
         ),

         fluidRow(
            column(width = 9, 
                   sliderInput("n_background", "Total Genes (background):",
                               value = 16000, min=5000, max=30000, round=TRUE)
                   ),
            column(width = 3,
                   fluidRow(numericInput(inputId = "n_background_box","", value=NA, min=5000, max=30000)),
                   fluidRow(align = "centre",  actionButton("setButton_n_background", "Set")) 
            )

         ),

         
         hr(),
         h4("Terms"),
      
         
         fluidRow(
            
            column(width = 9, sliderInput("n_term_genes", "Genes with Pathway Term:",
                                          value = 500, min=1, max=2000, round=TRUE)),
            column(width = 3,
                   fluidRow(numericInput(inputId = "n_term_genes_box","", value=NA, min=1, max=2000)),
                   fluidRow(align = "centre",  actionButton("setButton_n_term_genes", "Set")) 
            )

         ),
         

         
         
         
         fluidRow(
            
            column(width = 9, sliderInput("n_terms", "Number of Terms:",
                                          value = 10000, min=1, max=100000, round=TRUE)),
            column(width = 3,
                   fluidRow(numericInput(inputId = "n_terms_box","",value=NA,min=1, max=100000)),
                   fluidRow(align = "centre",  actionButton("setButton_n_terms", "Set")) 
            )
         ),
         
         hr(),

         selectInput("test", "Statistical test: ",
                     choices = c("Hypergeometric enrichment", "Fishers exact"))
         
      ),
      
      
      
      
      
      
      # Show a plot of the generated distribution
      mainPanel( width= 7,
         fluidRow( 
            column(width = 9 , align="center", 
                   sliderInput("n_hit", paste0("Number of genes hitting term '",pathway_name,"':"), 
                               min=0,  
                               round=TRUE, step=1,
                               max   = 300,
                               value = 30, 
                   )),
            column(width = 2 , align="left", #style = "margin-top: 20px;",
                   fluidRow(  numericInput("n_hit_box", "", 
                                min   = 0, 
                                max   = 300,
                                value = NA) 
                   ),
                   fluidRow(actionButton("setButton", "Set")) 
            )
            
         ), 
         
         
         br(),
         fluidRow(
            
            column(5,
               h5("Uncorrected p value = "),
               verbatimTextOutput("pval")
            ),
            column(5,
               h5("Corrected p value = "),
               verbatimTextOutput("pval_corr")
            )
         ),
         br(),
         plotOutput("vennPlot"),
         hr(),
         htmlOutput("prettyhtmlsummary")

      )
   ),
   
   # ( removes the up/down adjust arrows from the numeric boxes)
   tags$style(HTML("
        input[type=number] {
              -moz-appearance:textfield;
        }
        input[type=number]::{
              -moz-appearance:textfield;
        }
        input[type=number]::-webkit-outer-spin-button,
        input[type=number]::-webkit-inner-spin-button {
              -webkit-appearance: none;
              margin: 0;
        }
    "))
)


server <- function(input, output, session) {
   
   

   #----------------------------------------------------------------------------
   # Event handling.

   #Use obseve event to update the range of the n hit genes box
   # That box doesn't have reactive input
   observeEvent({input$n_test_genes | input$n_term_genes}, {
      max_allowed <- min(input$n_test_genes, input$n_term_genes)
      
      updateSliderInput(session, "n_hit", value=input$n_hit,   max = max_allowed)
      updateNumericInput(session,"n_hit_box", 
                         #also need to explicity change the range on input box:
                         value = min(input$n_hit_box, max_allowed), 
                         max   = max_allowed)
      
   })

   

   observeEvent(input$setButton, {
      if (! is.na(input$n_hit_box)) {
         updateNumericInput(session,"n_hit", value = input$n_hit_box)
         isolate(updateNumericInput(session,"n_hit_box", value=NA)) #rset
      }
   })   
   #https://stackoverflow.com/questions/45007730/update-goes-into-a-loop-shiny

   
   observeEvent(input$setButton_n_test_genes, {
      if (! is.na(input$n_test_genes_box)) {
         updateNumericInput(session,"n_test_genes", value = input$n_test_genes_box)
         isolate(updateNumericInput(session,"n_test_genes_box", value=NA)) #rset
      }
   }) 
   
   observeEvent(input$setButton_n_background, {
      if (! is.na(input$n_background_box)) {
         updateNumericInput(session,"n_background", value = input$n_background_box)
         isolate(updateNumericInput(session,"n_background_box", value=NA)) #rset
      }
   })  
   
   observeEvent(input$setButton_n_term_genes, {
      if (! is.na(input$n_term_genes_box)) {
         updateNumericInput(session,"n_term_genes", value = input$n_term_genes_box)
         isolate(updateNumericInput(session,"n_term_genes_box", value=NA)) #rset
      }
   })  
   
   observeEvent(input$setButton_n_terms, {
      if (! is.na(input$n_terms_box)) {
         updateNumericInput(session,"n_terms", value = input$n_terms_box)
         isolate(updateNumericInput(session,"n_terms_box", value=NA)) #rset
      }
   })  
   
   
   #----------------------------------------------------------------------------
   # Outputs
   output$prettyhtmlsummary<- renderUI({
      
      the_text<- paste0( 
         "There were ",strong(input$n_test_genes), " differentially expressed genes, ",
         "out of ", strong(input$n_background)," total background genes. ",
         br(),
         "Of those differentially expressed genes, ",
         strong(input$n_hit), 
         " (",em(percent(input$n_hit/input$n_test_genes)),") ",
         " were annotated with term '",pathway_name,"' (e.g. 'neurite outgrowth', 'TGFB signalling'). ",
         " The term has ",strong(input$n_term_genes)," genes in total; ",
         em(percent(input$n_term_genes/input$n_background))," of the background genes. ",
         br(),
         "The p-val calculated by ",
         strong(paste0(input$test, " test")),
         " is ", 
         em(formatC(calculated_pval(), format = "e", digits = 2)),
         ", and after multiple hypothesis correction for the ", strong(input$n_terms), " terms tested, ",
         "the corrected p-val is ", 
         em(formatC(p.adjust(calculated_pval(), n = input$n_terms ), format = "e", digits = 2))
      )
      HTML(the_text)
   })
   
   
   output$vennPlot <- renderPlot({
      vennfit <- euler(c(
         Query      = 0,
         Term       = 0,
         Background = input$n_background - input$n_test_genes - input$n_term_genes + input$n_hit,
         "Query&Term&Background" = input$n_hit, 
         "Query&Background" = input$n_test_genes - input$n_hit,
         "Term&Background"  = input$n_term_genes - input$n_hit
      ))
         
      plot(vennfit,
           fills = c("hotpink", "darkgoldenrod1", "grey90"),
           edges = FALSE,
           fontsize = 20,
           quantities = list(fontsize = 15))
      
   })
   
   calculated_pval <- reactive({
      
      #updateNumericInput(session,"n_hit_box",   value = input$n_hit)
      
      pval <- NA
      num_non_term_genes = input$n_background - input$n_term_genes
      
      if (input$test == "Hypergeometric enrichment" ) {
         pval = 1 - phyper(
            q = input$n_hit, 
            m = input$n_term_genes, 
            n = num_non_term_genes,
            k = input$n_test_genes,
            lower.tail = TRUE
         )
      }
      else if (input$test == "Fishers exact")  {

         n_not_test_or_hit <- input$n_background - input$n_test_genes - (input$n_term_genes - input$n_hit)
         contingency <- 
            matrix(c(input$n_hit, input$n_term_genes - input$n_hit,
                     input$n_test_genes - input$n_hit, n_not_test_or_hit ), 
                   nrow=2)
         
         pval=fisher.test(contingency, alternative = 'greater')$p.value
         
      }
      return(pval)
   })
   

   output$pval <- renderText({ 
      return(calculated_pval())
   })
   
   output$pval_corr <- renderText ({
      p.adjust(calculated_pval(), n = input$n_terms )
   })
   

   

   
   
   

}

# Run the application 
shinyApp(ui = ui, server = server)



