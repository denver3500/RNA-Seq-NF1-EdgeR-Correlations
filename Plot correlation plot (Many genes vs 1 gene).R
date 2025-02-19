library(dplyr)
library(tidyr)
library(SummarizedExperiment)
library(edgeR)


summarized_experiment <- readRDS("combined_se.rds")
edgeRlist <- SE2DGEList(summarized_experiment)

plotCorrelationPlot <- function(edgeRlist, tumor_type, tissue, gene_names, gene_to_compare) {
  
  # Create sample filter based on tumor type and tissue.
  sample_filter <- edgeRlist$samples$tumorType == tumor_type &
                   edgeRlist$samples$tissue == tissue
  
  # Subset samples and counts
  filtered_samples <- edgeRlist$samples[sample_filter, ]
  filtered_counts <- edgeRlist$counts[, sample_filter]
  filtered_edgeRlist <- edgeRlist
  filtered_edgeRlist$samples <- filtered_samples
  filtered_edgeRlist$counts <- filtered_counts
  
  # Compute normalization factors and CPM
  filtered_edgeRlist <- calcNormFactors(filtered_edgeRlist)
  filtered_cpm <- cpm(filtered_edgeRlist)
  
  # Ensure metadata has gene_name column.
  if (!"gene_name" %in% colnames(edgeRlist$genes)) {
    stop("The column 'gene_name' was not found in edgeRlist$genes")
  }
  
  # Retrieve gene IDs corresponding to the gene names.
  gene_ids <- rownames(edgeRlist$genes)[edgeRlist$genes$gene_name %in% gene_names]
  
  # Look up gene ID
  lookup_id <- rownames(edgeRlist$genes)[edgeRlist$genes$gene_name == gene_to_compare]
  if (length(lookup_id) == 0) {
    stop(paste("The gene_to_compare", gene_to_compare, "was not found in the metadata"))
  }
  gene_to_compare_id <- lookup_id[1]
  
  # Ensure gene_to_compare_id is included among the gene_ids.
  if (!(gene_to_compare_id %in% gene_ids)) {
    gene_ids <- c(gene_to_compare_id, gene_ids)
  }

  # Subset the filtered CPM matrix using the retrieved gene IDs.
  subset_cpm <- filtered_cpm[gene_ids, ]
  
  # Check that gene_to_compare_id is in the subset.
  if (!(gene_to_compare_id %in% rownames(subset_cpm))) {
    stop(paste("The gene", gene_to_compare, "was not found in the subset CPM data"))
  }
  
  # Extract expression profile for gene_to_compare.
  gene_expr <- subset_cpm[gene_to_compare_id, ]
  
  # Compute Pearson correlations of gene_to_compare with each gene.
  correlations <- apply(subset_cpm, 1, function(x) cor(gene_expr, x, method = "pearson"))
  correlations <- correlations[names(correlations) != gene_to_compare_id]
  print(correlations)
  
  # Create data frame of correlations.
  df_cor <- data.frame(
    gene_id = names(correlations),
    correlation = as.numeric(correlations)
  )
  
  # Change Gene IDs to Gene Names
  df_cor$gene <- edgeRlist$genes$gene_name[match(df_cor$gene_id, rownames(edgeRlist$genes))]
  
  # Define desired order from gene_names excluding gene_to_compare.
  desired_order <- rev(gene_names[gene_names != gene_to_compare])
  
  # Convert gene column to a factor
  df_cor$gene <- factor(df_cor$gene, levels = desired_order)
  
  # Load ggplot2 and generate the bar plot.
  library(ggplot2)
  p <- ggplot(df_cor, aes(x = gene, y = correlation)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    scale_y_continuous(limits = c(-1, 1)) +
    labs(title = paste("Correlation of", gene_to_compare, "in", tumor_type),
         x = "Gene",
         y = "Pearson Correlation") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # Build filename using tumor_type and gene_to_compare.
  file_name <- paste0("correlation_plot_", gsub(" ", "_", tumor_type), "_", gene_to_compare, ".png")
  
  # Save the plot.
  ggsave(file_name, plot = p, width = 8, height = 6, dpi = 300)
  
  return(p)
}

# Arguments for PN
tumor_type <- "Plexiform Neurofibroma"
tissue <- "primary tumor"
gene_names <- c("COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6",
                "CDH1", "TJP1", "DSP", "SNAI1", "SNAI2",
                "TWIST1", "LEF1", "CDH2", "VIM", "CTNNB1", "ACTA2", "TRAP1")
gene_to_compare <- "TRAP1"

# Function to plot for PN
plotCorrelationPlot(edgeRlist, tumor_type, tissue, gene_names, gene_to_compare)

# Arguments for MPNST
tumor_type <- "Malignant Peripheral Nerve Sheath Tumor"
tissue <- "primary tumor"
gene_names <- c("COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6",
                "CDH1", "TJP1", "DSP", "SNAI1", "SNAI2",
                "TWIST1", "LEF1", "CDH2", "VIM", "CTNNB1", "ACTA2", "TRAP1")
gene_to_compare <- "TRAP1"

#Function to plot for MPNST
plotCorrelationPlot(edgeRlist, tumor_type, tissue, gene_names, gene_to_compare)


