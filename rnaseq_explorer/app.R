library(shiny)
library(corrplot)
library(pheatmap)

# absolute paths should always be used when inputing data
#setwd("/Volumes/Data/academic/svmlab/ml/cancer_projects/aml/OHSU/data/")
drug_screen <- read.csv("drug.csv", stringsAsFactors = FALSE, sep=",", check.names = FALSE)
AUC_All <- tapply(drug_screen[,"Inhibitor Interpreted Result: Area under the curve"], list(r = drug_screen[,"Specimen: Lab ID"], c = drug_screen[,"Inhibitor Interpreted Result: Drug"]), mean)
IC50_All <- tapply(drug_screen[,"Inhibitor Interpreted Result: IC50"], list(r = drug_screen[,"Specimen: Lab ID"], c = drug_screen[,"Inhibitor Interpreted Result: Drug"]), mean)

drug_screen <- drug_screen[which(drug_screen[,"Heme Malignancy: Diagnosis"] == "ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS"),]
AUC_AML <- tapply(drug_screen[,"Inhibitor Interpreted Result: Area under the curve"], list(r = drug_screen[,"Specimen: Lab ID"], c = drug_screen[,"Inhibitor Interpreted Result: Drug"]), mean)
IC50_AML <- tapply(drug_screen[,"Inhibitor Interpreted Result: IC50"], list(r = drug_screen[,"Specimen: Lab ID"], c = drug_screen[,"Inhibitor Interpreted Result: Drug"]), mean)

BeatAML_DD_All <- list()
BeatAML_DD_All$AUC <- AUC_All
BeatAML_DD_All$IC50 <- IC50_All
BeatAML_DD_All$AUC_AUC_Pearson <- cor(AUC_All, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DD_All$AUC_AUC_Spearman <- cor(AUC_All, method = "spearman", use = "pairwise.complete.obs")
BeatAML_DD_All$IC50_IC50_Pearson <- cor(IC50_All, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DD_All$IC50_IC50_Spearman <- cor(IC50_All, method = "spearman", use = "pairwise.complete.obs")
# save(BeatAML_DD_All, file = "BeatAML_DD_All.RData")

BeatAML_DD_AML <- list()
BeatAML_DD_AML$AUC <- AUC_AML
BeatAML_DD_AML$IC50 <- IC50_AML
BeatAML_DD_AML$AUC_AUC_Pearson <- cor(AUC_AML, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DD_AML$AUC_AUC_Spearman <- cor(AUC_AML, method = "spearman", use = "pairwise.complete.obs")
BeatAML_DD_AML$IC50_IC50_Pearson <- cor(IC50_AML, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DD_AML$IC50_IC50_Spearman <- cor(IC50_AML, method = "spearman", use = "pairwise.complete.obs")
# save(BeatAML_DD_AML, file = "BeatAML_DD_AML.RData")

RNASeq <- read.csv("exp.csv", stringsAsFactors = FALSE, sep=",", row.names = "gene", check.names = FALSE)
RNASeq <- log2(t(RNASeq) + 1)

common_patients_All <- intersect(rownames(IC50_All), rownames(RNASeq))
common_patients_AML <- intersect(rownames(IC50_AML), rownames(RNASeq))

gene_names <- c("CTLA-4", "Galectin-9", "PD-L1", "PD-L2", "PD1", "TIM-3", "VISTA", "LAG3", "TNFRSF4", "CD3G", "ARG1", "ARG2")
gene_ids <- c("CTLA4", "LGALS9", "CD274", "PDCD1LG2", "PDCD1", "HAVCR2", "C10orf54", "LAG3", "TNFRSF4", "CD3G", "ARG1", "ARG2")
RNASeq_All <- RNASeq[common_patients_All, gene_ids]
RNASeq_AML <- RNASeq[common_patients_AML, gene_ids]
colnames(RNASeq_All) <- gene_names
colnames(RNASeq_AML) <- gene_names

AUC_All <- AUC_All[common_patients_All,]
IC50_All <- IC50_All[common_patients_All,]
#selected_drugs <- which(colSums(!is.na(AUC_All)) >= 6)
#AUC_All <- AUC_All[,selected_drugs]
#IC50_All <- IC50_All[,selected_drugs]
AUC_AML <- AUC_AML[common_patients_AML,]
IC50_AML <- IC50_AML[common_patients_AML,]
#selected_drugs <- which(colSums(!is.na(AUC_AML)) >= 6)
#AUC_AML <- AUC_AML[,selected_drugs]
#IC50_AML <- IC50_AML[,selected_drugs]

BeatAML_DG_All <- list()
BeatAML_DG_All$AUC <- AUC_All
BeatAML_DG_All$IC50 <- IC50_All
BeatAML_DG_All$RNASeq <- RNASeq_All
BeatAML_DG_All$AUC_RNASeq_Pearson <- cor(AUC_All, RNASeq_All, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DG_All$AUC_RNASeq_Spearman <- cor(AUC_All, RNASeq_All, method = "spearman", use = "pairwise.complete.obs")
BeatAML_DG_All$IC50_RNASeq_Pearson <- cor(IC50_All, RNASeq_All, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DG_All$IC50_RNASeq_Spearman <- cor(IC50_All, RNASeq_All, method = "spearman", use = "pairwise.complete.obs")
# save(BeatAML_DG_All, file = "BeatAML_DG_All.RData")

BeatAML_DG_AML <- list()
BeatAML_DG_AML$AUC <- AUC_AML
BeatAML_DG_AML$IC50 <- IC50_AML
BeatAML_DG_AML$RNASeq <- RNASeq_AML
BeatAML_DG_AML$AUC_RNASeq_Pearson <- cor(AUC_AML, RNASeq_AML, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DG_AML$AUC_RNASeq_Spearman <- cor(AUC_AML, RNASeq_AML, method = "spearman", use = "pairwise.complete.obs")
BeatAML_DG_AML$IC50_RNASeq_Pearson <- cor(IC50_AML, RNASeq_AML, method = "pearson", use = "pairwise.complete.obs")
BeatAML_DG_AML$IC50_RNASeq_Spearman <- cor(IC50_AML, RNASeq_AML, method = "spearman", use = "pairwise.complete.obs")
# save(BeatAML_DG_AML, file = "BeatAML_DG_AML.RData")

server <- function(input, output){
	output$dg_corr_legend <- renderPlot({
    legend(x = "center", legend = c("-1.00", "-0.50", "+0.50", "+1.00"), col = c(corr_colors[1], corr_colors[51], corr_colors[151], corr_colors[200]), pch = c(19, 19), pt.cex = c(8, 5, 5, 8), cex = 2, horiz = TRUE, box.lwd = 0)
  }) 
  output$dg_corr_plot <- renderPlot({
    if (input$dg_corr_disease == "AML") {
      BeatAML <- BeatAML_DG_AML
    } else {
      BeatAML <- BeatAML_DG_All
    }
    if (input$dg_corr_sensitivity == "AUC" && input$dg_corr_correlation == "Pearson") {
      X <- BeatAML$AUC_RNASeq_Pearson[, input_genes, drop = FALSE]
    }
    if (input$dg_corr_sensitivity == "AUC" && input$dg_corr_correlation == "Spearman") {
      X <- BeatAML$AUC_RNASeq_Spearman[, input_genes, drop = FALSE]
    }
    if (input$dg_corr_sensitivity == "IC50" && input$dg_corr_correlation == "Pearson") {
      X <- BeatAML$IC50_RNASeq_Pearson[, input_genes, drop = FALSE]
    }
    if (input$dg_corr_sensitivity == "IC50" && input$dg_corr_correlation == "Spearman") {
      X <- BeatAML$IC50_RNASeq_Spearman[, input_genes, drop = FALSE]
    }    
    row_order <- sort(rowSums(X^2, na.rm = TRUE), decreasing = TRUE, index.return = TRUE)$ix
    col_order <- sort(colSums(X^2, na.rm = TRUE), decreasing = TRUE, index.return = TRUE)$ix
    X <- X[row_order, col_order, drop = FALSE]
    X <- X[which(rowMeans(is.na(X)) != 1),]
    X[is.na(X)] <- 0
    corrplot(X, tl.col = "black", tl.cex = 1.25, method = "circle", cl.pos = "n", col = corr_colors)
  }, height = 1240 * 5, width = 800)
  output$dg_corr_table <- renderDataTable({
    if (input$dg_corr_disease == "AML") {
      BeatAML <- BeatAML_DG_AML
    } else {
      BeatAML <- BeatAML_DG_All
    }
    if (input$dg_corr_sensitivity == "AUC" && input$dg_corr_correlation == "Pearson") {
      X <- BeatAML$AUC_RNASeq_Pearson[, input_genes, drop = FALSE]
    }
    if (input$dg_corr_sensitivity == "AUC" && input$dg_corr_correlation == "Spearman") {
      X <- BeatAML$AUC_RNASeq_Spearman[, input_genes, drop = FALSE]
    }
    if (input$dg_corr_sensitivity == "IC50" && input$dg_corr_correlation == "Pearson") {
      X <- BeatAML$IC50_RNASeq_Pearson[, input_genes, drop = FALSE]
    }
    if (input$dg_corr_sensitivity == "IC50" && input$dg_corr_correlation == "Spearman") {
      X <- BeatAML$IC50_RNASeq_Spearman[, input_genes, drop = FALSE]
    }
    X <- X[which(rowMeans(is.na(X)) != 1),]
    X <- cbind.data.frame(as.vector(matrix(rownames(X), nrow = nrow(X), ncol = ncol(X), byrow = FALSE)), as.vector(matrix(colnames(X), nrow = nrow(X), ncol = ncol(X), byrow = TRUE)), as.vector(X))
    colnames(X) <- c("Drug", "Gene", "Correlation")
    X[is.na(X)] <- 0
    X <- X[sort(abs(X[,"Correlation"]), decreasing = TRUE, index.return = TRUE, na.last = NA)$ix,]
    return (X)
  }, options = list(searching = FALSE, paging = FALSE))
  output$dg_heat_plot <- renderPlot({
    if (input$dg_heat_disease == "AML") {
      BeatAML <- BeatAML_DG_AML
    } else {
      BeatAML <- BeatAML_DG_All
    }
    if (input$dg_heat_sensitivity == "AUC") {
      X_drug <- BeatAML$AUC[, input$dg_heat_drug, drop = FALSE]
    }
    if (input$dg_heat_sensitivity == "IC50") {
      X_drug <- BeatAML$IC50[, input$dg_heat_drug, drop = FALSE]
    }
    X_expression <- BeatAML$RNASeq[, input$dg_heat_gene, drop = FALSE]
    X <- cbind(X_drug, X_expression)
    
    X_top <- X[which(!is.na(X[,1])),]
    X_top <- X_top[sort(X_top[,1], decreasing = FALSE, index.return = TRUE)$ix,]
    X_bottom <- X[which(is.na(X[,1])),]
    X_bottom <- X_bottom[sort(X_bottom[,2], decreasing = FALSE, index.return = TRUE)$ix,]
    
    pheatmap(scale(X_top), fontsize = 16, fontsize_col = 20, fontsize_row = 16, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE)
  }, height = 1200)
  output$dd_corr_legend <- renderPlot({
    legend(x = "center", legend = c("-1.00", "-0.50", "+0.50", "+1.00"), col = c(corr_colors[1], corr_colors[51], corr_colors[151], corr_colors[200]), pch = c(19, 19), pt.cex = c(2.5, 1.5, 1.5, 2.5), cex = 2, horiz = TRUE, box.lwd = 0)
  })
  output$dd_corr_plot <- renderPlot({
    if (input$dd_corr_disease == "AML") {
      BeatAML <- BeatAML_DD_AML
    } else {
      BeatAML <- BeatAML_DD_All
    }
    if (input$dd_corr_sensitivity == "AUC" && input$dd_corr_correlation == "Pearson") {
      X <- BeatAML$AUC_AUC_Pearson
    }
    if (input$dd_corr_sensitivity == "AUC" && input$dd_corr_correlation == "Spearman") {
      X <- BeatAML$AUC_AUC_Spearman
    }
    if (input$dd_corr_sensitivity == "IC50" && input$dd_corr_correlation == "Pearson") {
      X <- BeatAML$IC50_IC50_Pearson
    }
    if (input$dd_corr_sensitivity == "IC50" && input$dd_corr_correlation == "Spearman") {
      X <- BeatAML$IC50_IC50_Spearman
    }
    X[is.na(X)] <- 0
    corrplot(X, tl.col = "black", tl.cex = 1.25, method = "circle", cl.pos = "n", order = "hclust", col = corr_colors)
  }, height = 1240 * 2, width = 1240 * 2)
  output$dd_corr_table <- renderDataTable({
    if (input$dd_corr_disease == "AML") {
      BeatAML <- BeatAML_DD_AML
    } else {
      BeatAML <- BeatAML_DD_All
    }
    if (input$dd_corr_sensitivity == "AUC" && input$dd_corr_correlation == "Pearson") {
      X <- BeatAML$AUC_AUC_Pearson
    }
    if (input$dd_corr_sensitivity == "AUC" && input$dd_corr_correlation == "Spearman") {
      X <- BeatAML$AUC_AUC_Spearman
    }
    if (input$dd_corr_sensitivity == "IC50" && input$dd_corr_correlation == "Pearson") {
      X <- BeatAML_DD$IC50_IC50_Pearson
    }
    if (input$dd_corr_sensitivity == "IC50" && input$dd_corr_correlation == "Spearman") {
      X <- BeatAML_DD$IC50_IC50_Spearman
    }
    X[lower.tri(X)] <- NA
    X <- cbind.data.frame(as.vector(matrix(rownames(X), nrow = nrow(X), ncol = ncol(X), byrow = FALSE)), as.vector(matrix(colnames(X), nrow = nrow(X), ncol = ncol(X), byrow = TRUE)), as.vector(X))
    colnames(X) <- c("Drug #1", "Drug #2", "Correlation")
    X <- X[which(!is.na(X[, 3])), ]
    X <- X[which(X[,1] != X[,2]),]
    X <- X[sort(abs(X[,"Correlation"]), decreasing = TRUE, index.return = TRUE)$ix,]
    return (X)
  }, options = list(searching = FALSE, paging = TRUE, pageLength = 1000, lengthMenu = list(c(1000, 5000, -1), c('1K', '5K', 'All'))))
  output$dd_heat_plot <- renderPlot({
    if (input$dd_heat_disease == "AML") {
      BeatAML <- BeatAML_DD_AML
    } else {
      BeatAML <- BeatAML_DD_All
    }
    if (input$dd_heat_sensitivity == "AUC") {
      X_drug1 <- BeatAML$AUC[, input$dd_heat_drug1, drop = FALSE]
    }
    if (input$dd_heat_sensitivity == "IC50") {
      X_drug1 <- BeatAML$IC50[, input$dd_heat_drug1, drop = FALSE]
    }
    if (input$dd_heat_sensitivity == "AUC") {
      X_drug2 <- BeatAML$AUC[, input$dd_heat_drug2, drop = FALSE]
    }
    if (input$dd_heat_sensitivity == "IC50") {
      X_drug2 <- BeatAML$IC50[, input$dd_heat_drug2, drop = FALSE]
    }
    X <- cbind(X_drug1, X_drug2)
    
    X_top <- X[which(!is.na(X[,1])),]
    X_top <- X_top[sort(X_top[,1], decreasing = FALSE, index.return = TRUE)$ix,]
    X_bottom <- X[which(is.na(X[,1])),]
    X_bottom <- X_bottom[sort(X_bottom[,2], decreasing = FALSE, index.return = TRUE)$ix,]
    
    pheatmap(scale(X_top), fontsize = 16, fontsize_col = 20, fontsize_row = 16, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE)
  }, height = 1200 * 2)
}

ui <- navbarPage(title = "BeatAML Explorer", fluid = FALSE, position = "static-top",
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

shinyApp(ui = ui, server = server)