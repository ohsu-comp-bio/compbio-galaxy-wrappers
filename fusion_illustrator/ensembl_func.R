# custom functions --------------------------------------------------------
#functions to pull out ENST ID for use in code
#Sometimes an ENST is not provided for the most abundant breakpoint discovered by STAR-fusion
#this breaks the code
#instead scan down the STAR-Fusion list to try to find something to plot.
#If that does not work pick the first ENSTID on the Ensembl list
getLeftENST_function <- function(fusions_df) {
  if(all(fusions_df$ENST_LEFT_ID %in% "." == "TRUE")){
    #if there is no ENSTID for this gene pair
    #then just pick the first ENSTID from ensembl db
    #I'm guessing that will be the longest transcript

    gtf %>%
      dplyr::filter(gene_name %in% left_gene_name) %>%
      dplyr::pull(tx_id) %>%
      dplyr::first()

    #Maybe I should add an asterisk for this case?


  } else {
    #otherwise, get the ENST for the most abundant breakpoint that has one

    left_ENST <- fusions_df %>%
      dplyr::filter(Left.Gene %in% left_gene_name) %>% #probably redundant, but is good check?
      dplyr::filter(!ENST_LEFT_ID %in% ".") %>%
      dplyr::pull(ENST_LEFT_ID) %>%
      dplyr::first()

    left_ENST <- strsplit(left_ENST, ".", fixed = TRUE)[[1]][1]
  }
}
getRightENST_function <- function(fusions_df) {
  if(all(fusions_df$ENST_RIGHT_ID %in% "." == "TRUE")){
    #if there is no ENSTID for this gene pair
    #then just pick the first ENSTID from ensembl db
    #I'm guessing that will be the longest transcript

    gtf %>%
      dplyr::filter(gene_name %in% right_gene_name) %>%
      dplyr::pull(tx_id) %>%
      dplyr::first()


  } else { #otherwise, get the ENST for the most abundant breakpoint that has one

    right_ENST <- fusions_df %>%
      dplyr::filter(Right.Gene %in% right_gene_name) %>% #probably redundant, but is good check?
      dplyr::filter(!ENST_RIGHT_ID %in% ".") %>%
      dplyr::pull(ENST_RIGHT_ID) %>%
      dplyr::first()

    right_ENST <- strsplit(right_ENST, ".", fixed = TRUE)[[1]][1]
  }
}