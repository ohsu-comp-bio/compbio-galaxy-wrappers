library(shiny)

load("BeatAML_DD_All.RData")
load("BeatAML_DG_All.RData")
load("BeatAML_DD_AML.RData")
load("BeatAML_DG_AML.RData")

shinyUI(
  navbarPage(title = "BeatAML Explorer", fluid = FALSE, position = "static-top",
    tabPanel(title = "Drug-Gene Correlation", 
      sidebarLayout(
        sidebarPanel(width = 3,
          selectInput(inputId = 'dg_corr_disease', label = 'Disease', choices = c("AML", "All Diseases")),
          selectInput(inputId = 'dg_corr_sensitivity', label = 'Sensitivity', choices = c("AUC", "IC50")),
          selectInput(inputId = 'dg_corr_correlation', label = 'Function', choices = c("Pearson", "Spearman"))
        ),
        mainPanel(
          helpText('Higher AUC or IC50 values mean more resistance.'),
          helpText('If the correlation between drug sensitivity and gene expression values is positive, it means that resistant patients have relatively higher expression values of the gene.'),
          helpText('If the correlation between drug sensitivity and gene expression values is negative, it means that sensitive patients have relatively higher expression values of the gene.'),
          tabsetPanel(type = "tabs", 
                      tabPanel(title = "Plot", plotOutput(outputId = "dg_corr_legend", height = 50, width = 800), plotOutput(outputId = "dg_corr_plot")), 
                      tabPanel(title = "Table", dataTableOutput(outputId = "dg_corr_table"))
          )
        )
      )
    ),
    tabPanel(title = "Drug-Gene Heatmap",
      sidebarLayout(
        sidebarPanel(width = 3,
                     selectInput(inputId = 'dg_heat_disease', label = 'Disease', choices = c("AML", "All Diseases")),
                     selectInput(inputId = 'dg_heat_drug', label = 'Drug', choices = colnames(BeatAML_DG_AML$IC50)),
                     selectInput(inputId = 'dg_heat_sensitivity', label = 'Sensitivity', choices = c("AUC", "IC50")),
                     selectInput(inputId = 'dg_heat_gene', label = 'Gene', choices = colnames(BeatAML_DG_AML$RNASeq))
        ),
        mainPanel(
          helpText('Higher AUC or IC50 values mean more resistance.'),
          plotOutput(outputId = "dg_heat_plot")
        ) 
      )
    ),
    tabPanel(title = "Drug-Drug Correlation", 
             sidebarLayout(
               sidebarPanel(width = 3,
                            selectInput(inputId = 'dd_corr_disease', label = 'Disease', choices = c("AML", "All Diseases")),
                            selectInput(inputId = 'dd_corr_sensitivity', label = 'Sensitivity', choices = c("AUC", "IC50")),
                            selectInput(inputId = 'dd_corr_correlation', label = 'Function', choices = c("Pearson", "Spearman"))
               ),
               mainPanel(
                 helpText('If the correlation between drug sensitivity values of two drugs is positive, it means that patients have similar responses against these two drugs.'),
                 helpText('If the correlation between drug sensitivity values of two drugs is negative, it means that patients have different responses against these two drugs.'),
                 tabsetPanel(type = "tabs", 
                             tabPanel(title = "Plot", plotOutput(outputId = "dd_corr_legend", height = 50, width = 800), plotOutput(outputId = "dd_corr_plot")), 
                             tabPanel(title = "Table", dataTableOutput(outputId = "dd_corr_table"))
                 )
               )
             )
    ),
    tabPanel(title = "Drug-Drug Heatmap",
      sidebarLayout(
       sidebarPanel(width = 3,
                    selectInput(inputId = 'dd_heat_disease', label = 'Disease', choices = c("AML", "All Diseases")),
                    selectInput(inputId = 'dd_heat_drug1', label = 'Drug #1', choices = colnames(BeatAML_DD_AML$IC50)),
                    selectInput(inputId = 'dd_heat_drug2', label = 'Drug #2', choices = colnames(BeatAML_DD_AML$IC50)),
                    selectInput(inputId = 'dd_heat_sensitivity', label = 'Sensitivity', choices = c("AUC", "IC50"))
       ),
       mainPanel(
         helpText('Higher AUC or IC50 values mean more resistance.'),
         plotOutput(outputId = "dd_heat_plot")
       ) 
      )
    )
  )
)