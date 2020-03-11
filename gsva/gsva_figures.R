library(ggplot2)
library(pheatmap)


# Rename the columns of a melted table
rename.melted <- function(input, feature=NULL, sample=NULL, value=NULL, in.place=FALSE) {
    # Set up the data for use in-or-out-of place
    if (in.place) {
        data <- input
    } else {
        data <- data.frame(input)
    }
    cols <- colnames(data)
    # Set the feature column
    if (!is.null(feature)) {
        if (!(feature %in% names(data))) {
            stop(paste0('Error in table_visualization::preprocess - no column ',
                        feature, ' in data'))
        }
        cols[cols==feature] <- 'Feature'
    }
    # Set the sample column
    if (!is.null(sample)) {
        if (!(sample %in% names(data))) {
            stop(paste0('Error in table_visualization::preprocess - no column ',
                        sample, ' in data'))
        }
        cols[cols==feature] <- 'Sample'
    }
    # Set the values column
    if (!is.null(value)) {
        if (!(value %in% names(data))) {
            stop(paste0('Error in table_visualization::preprocess - no column ',
                        value, ' in data'))
        }
        cols[cols==feature] <- 'Value'
    }
    # Assign the new column names
    colnames(data) <- cols
    return(data)
}

# Produce a box and whisker plot for the cohort
make.box <- function(data, sample, feature.name=NULL, sample.name=NULL, value.name=NULL) {
    df = data.frame(data)
    cols <- colnames(df)
    # Establish feature column name
    if (is.null(feature.name)) {
        if(!('Feature' %in% colnames(data))){
            stop(paste0('Error in table_visualization::make.box - no Feature column in data with columns ',colnames(data)))
        }
        feature.name <- "Feature"
    }
    # Establish value column name
    if (is.null(value.name)) {
        if(!('Value' %in% colnames(data))){
            stop(paste0('Error in table_visualization::make.box - no Value column in data with columns ',colnames(data)))
        }
        value.name <- "Value"
    }
    # Establish sample column name
    if (is.null(sample.name)) {
        if(!('Sample' %in% colnames(data))){
            stop(paste0('Error in table_visualization::make.box - no Sample column in data with columns ',colnames(data)))
        }
        sample.name <- "Sample"
    }
    colnames(df)[cols==sample.name] <- "Sample"
    colnames(df)[cols==value.name] <- "Value"
    colnames(df)[cols==feature.name] <- "Feature"
    
    # Subset the sample from the dataframe
    sample.df <- df[df$Sample==sample,]
    #stop(order(sample.df$Value))

    # Order the factors
    fn <- sample.df[order(sample.df$Value),]
    df$Feature <- factor(df$Feature, fn$Feature)
    sample.df$Feature <- factor(sample.df$Feature, fn$Feature)

    # Make the plot
    cols = c("#F8766D")
    g <- ggplot(df) +
      geom_boxplot(ggplot2::aes(y=Value, x=Feature)) +
      stat_boxplot(data=sample.df,aes(y=Value, x=Feature, col=Sample)) +
      coord_flip() +
      theme_minimal() + 
      scale_colour_manual(name="Sample", values=cols) +
      labs(x=value.name, y=feature.name)

    # Return the plot
    return(g)
}

# Produce a clustered heatmap for the data
make.heatmap <- function(data, colors=colorRampPalette(c("#2166ac", "#FFFFFF", "#b2182b")), breaks=seq(-1,1,by=0.02), size=4) {
    hm <- pheatmap(data, 
                   color=colors(length(breaks)),
                   breaks=breaks, 
                   cellwidth=size, 
                   cellheight=size, 
                   fontsize=size,
                   cluster_cols=TRUE, 
                   border_color=NA)
    return(hm)
}
