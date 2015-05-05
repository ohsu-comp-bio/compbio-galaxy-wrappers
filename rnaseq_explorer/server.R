library(shiny)
library(corrplot)
library(pheatmap)

load("BeatAML_DD_All.RData")
load("BeatAML_DG_All.RData")
load("BeatAML_DD_AML.RData")
load("BeatAML_DG_AML.RData")
input_genes <- c("CTLA-4", "Galectin-9", "PD-L1", "PD-L2", "PD1", "TIM-3", "VISTA", "LAG3", "TNFRSF4", "CD3G", "ARG1", "ARG2")

colors <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
corr_colors <- rev(colors(201))

shinyServer(function(input, output) {
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
})