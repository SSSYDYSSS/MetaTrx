#' Process Core Enrichment Genes and Count Positive/Negative Correlations
#'
#' This function processes core enrichment genes from a given string, identifies unique and non-unique genes,
#' and counts the number of positive and negative correlations for both unique and non-unique gene sets.
#'
#' @param core_enrichment A string of core enrichment genes separated by "/".
#' @param correlation_data A data frame containing gene correlations with columns "gene2" and "cor".
#' @examples
#' \dontrun{
#' core_enrichment <- "gene1/gene2/gene3/gene1"
#' correlation_data <- data.frame(gene2 = c("gene1", "gene2", "gene3", "gene4"),
#'                                cor = c(0.5, -0.3, 0.7, 0.2))
#' process_core_enrichment(core_enrichment, correlation_data)
#' }
#' @export
process_core_enrichment <- function(core_enrichment, correlation_data) {
  # Separate all genes and merge them into a vector
  core_enrichment_list <- unlist(strsplit(core_enrichment, "/"))
  core_enrichment_list_unique_genes <- unique(core_enrichment_list)
  core_enrichment_list_notunique <- length(core_enrichment_list)
  core_enrichment_list_unique <- length(core_enrichment_list_unique_genes)

  cat("#=================================================================================#", "\n")
  cat("core_enrichment_list_unique_genes:", "\n")
  print(core_enrichment_list_unique_genes)

  # Unique
  # Filter out genes in the gene2 column that match core_enrichment_list_unique_genes
  matching_genes_unique <- correlation_data[correlation_data$gene2 %in% core_enrichment_list_unique_genes, ]
  # Count the number of genes with positive and negative correlations
  positive_cor_genes_unique <- sum(matching_genes_unique$cor > 0)
  negative_cor_genes_unique <- sum(matching_genes_unique$cor < 0)

  # Not unique
  # Create a data frame to count the occurrences of each gene
  gene_count <- table(core_enrichment_list)
  # Initialize an empty data frame for storing the matching genes
  matching_genes_notunique <- data.frame()
  # Loop over each gene and its count to extract matching genes
  for (gene in names(gene_count)) {
    count <- gene_count[gene]
    matches <- correlation_data[correlation_data$gene2 == gene, ]
    # Repeat the rows according to the count and add to the matching_genes data frame
    matching_genes_notunique <- rbind(matching_genes_notunique, do.call("rbind", replicate(count, matches, simplify = FALSE)))
  }
  # Count the number of genes with positive and negative correlations
  positive_cor_genes_notunique <- sum(matching_genes_notunique$cor > 0)
  negative_cor_genes_notunique <- sum(matching_genes_notunique$cor < 0)


  # Print the results
  cat("#=================================================================================#", "\n")
  cat("All Number of genes with core enrichment (notunique):", core_enrichment_list_notunique, "\n")
  cat("All Number of genes with core enrichment (unique):", core_enrichment_list_unique, "\n")
  cat("#=================================================================================#", "\n")
  cat("Number of genes with positive correlation (notunique):", positive_cor_genes_notunique, "\n")
  cat("Number of genes with negative correlation (notunique):", negative_cor_genes_notunique, "\n")
  cat("#=================================================================================#", "\n")
  cat("Number of genes with positive correlation (unique):", positive_cor_genes_unique, "\n")
  cat("Number of genes with negative correlation (unique):", negative_cor_genes_unique, "\n")
}


