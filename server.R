library(shiny)
library(tidyverse) #for data wrangling
library(plotly)    #for interactive plots
library(limma)     #for making the model

default_in_file <-
  data.frame(
    name = "template.csv",
    size = 20618L,
    type = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    datapath = "www/template.csv",
    stringsAsFactors = FALSE
  )

instr <- "Select the conditions to compare. Gene groups can be switched on and 
off by clicking on them in the right-side panel. <br><br> Hover over points to 
view values or over the graph for more options or to save the plot."

shinyServer(function(input, output) {
  v <- reactiveValues(inFile = NULL)
  observeEvent(input$settempl, {v$inFile <- default_in_file})
  observeEvent(input$fo, {v$inFile <- input$fo})
  dat <- reactive({ #reacts to file being uploaded
    if (is.null(v$inFile)) NULL
    else{ #Read in data and extract relevant info
      dat_in <- read.csv(v$inFile$datapath, na = "0")
      values <- dat_in %>%  'rownames<-'(.[, 1]) %>% .[-1,-c(1,2)] %>% 
       data.matrix() %>%  '/'(colSums(., na.rm = TRUE)[col(.)]) * 100 %>% log2()  
      cond_groups <- factor(unlist(dat_in[1, -c(1,2)]))
      gene_groups <- factor(dat_in[-1, "gene_group"])
      list(values = values, cond_groups = cond_groups, gene_groups= gene_groups)
      }
    })

  drops <- reactive({ #reacts to dropdown values
    d1 <- input$c1
    d2 <- input$c2
    if (is.null(d1) || is.null(d2) || d1 == d2) NULL
    else TRUE  #TRUE is 2 different values are selected in dropdowns
    })
  
  #dropdowns are shown after data is uploaded
  output$dropdown1 <- renderUI({ 
    if (is.null(dat())) NULL
    else{cond_groups <- dat()[["cond_groups"]]
      selectInput(
      "c1", "compare these 2 conditions", choices = c(levels(cond_groups)), 
        selected = c(levels(cond_groups))[1])
        }
    })
  
  output$dropdown2 <- renderUI({
    if (is.null(dat())) NULL
    else{cond_groups <- dat()[["cond_groups"]]
         selectInput("c2", "", choices = c(levels(cond_groups)), 
           selected = c(levels(cond_groups))[2])
        }
    })
  
  #further instructions are shown after data is uploaded and conditions chosen
  output$instruction <- renderText({
    if (is.null(drops())) NULL
    else instr
  })
  
  output$volcplot <- renderPlotly({ #when drops==TRUE returns volcano plot
    if (is.null(drops())) return(NULL)
    else{
      #read in relevant data
      values <- dat()[["values"]]
      cond_groups_orig <- dat()[["cond_groups"]]  #copy of original group names
      gene_groups <- dat()[["gene_groups"]]
      groups <- cond_groups_orig %>%  #make syntactically valid names
          mapvalues(levels(.), make.names(levels(.), unique = TRUE)) 
      my_c1 <- input$c1
      my_c2 <- input$c2
      c1 <- groups[match(my_c1, cond_groups_orig)]  #match syntact valid names
      c2 <- groups[match(my_c2, cond_groups_orig)]
      
      #create linear fit and topTable on values
      ctrsts <- paste(c2, c1, sep = "-")   #define LIMMA contrasts 
      design <- model.matrix( ~ 0 + groups)
      colnames(design) <- levels(groups)
      contrast.matrix <- makeContrasts(contrasts = ctrsts, levels = levels(groups))
      fit <- lmFit(values, design) %>% contrasts.fit(contrast.matrix) %>% eBayes()
      plot_df <- topTable(fit, coef = 1, number = nrow(values), adjust = "BH") %>%
      mutate(gene = rownames(.), gene_groups = gene_groups)
      maxFC <- max(abs(plot_df$logFC), na.rm = T)
      maxPval <- max(-log10(plot_df$P.Value), na.rm = T)
        
      #return volcano plot from topTable
      gg <- ggplot(data = plot_df, aes(x = logFC, y = -log10(P.Value),
                                  color = gene_groups, text = plot_df$gene)) +
            labs(
          x = paste("up in", my_c1, "<---", "||logFC||", "---> up in", my_c2),
          y = "-10log P-value",
          title = "Volcano plot comparing the 2 selected conditions"
        ) +
        lims(x = c(-maxFC, maxFC), y = c(0, maxPval)) +
        geom_vline(xintercept = 0, linetype = 2) +
        theme(  #! theme minimal
          legend.key = element_rect(fill = "white"),
          aspect.ratio = 1, panel.background = element_rect(fill = "white", colour = "black")
          ) +
        scale_colour_brewer(palette = "Set1") +
        geom_point(size = 2, alpha = 0.7)
      print(gg)
      p <- ggplotly(gg)
      p
      }
    })
})