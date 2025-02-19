library(dplyr)
library(tidyr)
library(SummarizedExperiment)
library(edgeR)
library(corrplot)

summarized_experiment <- readRDS("combined_se.rds")
edgeRlist <- SE2DGEList(summarized_experiment)
plotSquareCorrPlot <- function(edgeRlist, tumor_type, tissue, gene_names) {
  
  # Filter samples
  sample_filter <- edgeRlist$samples$tumorType == tumor_type &
                   edgeRlist$samples$tissue == tissue
  filtered_samples <- edgeRlist$samples[sample_filter, ]
  filtered_counts <- edgeRlist$counts[, sample_filter]
  filtered_edgeRlist <- edgeRlist
  filtered_edgeRlist$samples <- filtered_samples
  filtered_edgeRlist$counts <- filtered_counts

  # Normalize and compute CPM.
  filtered_edgeRlist <- calcNormFactors(filtered_edgeRlist)
  filtered_cpm <- cpm(filtered_edgeRlist)
  
  # Check that the metadata has gene_name
  if(!"gene_name" %in% colnames(edgeRlist$genes)) {
    stop("The column 'gene_name' was not found in edgeRlist$genes")
  }
  
  # Retrieve gene IDs corresponding to the gene_names vector.
  gene_ids <- rownames(edgeRlist$genes)[edgeRlist$genes$gene_name %in% gene_names]
  
  # Subset the filtered CPM matrix based on gene_ids.
  subset_cpm <- filtered_cpm[gene_ids, ]
  
  # Compute the pairwise Pearson correlation matrix.
  cor_matrix <- cor(t(subset_cpm), method = "pearson")
  
  # Change Gene IDs to Gene Names.
  gene_labels_all <- edgeRlist$genes$gene_name[match(rownames(cor_matrix), rownames(edgeRlist$genes))]
  rownames(cor_matrix) <- gene_labels_all
  colnames(cor_matrix) <- gene_labels_all
  
  
  # Create filename. Replace spaces with underscores.
  file_name <- paste0("corrplot_", gsub(" ", "_", tumor_type), ".png")
  png(filename = file_name, width = 1000, height = 800)
  
  corrplot(cor_matrix, method = "color", 
           title = paste("Correlation Comparison for", tumor_type),
           mar = c(0, 0, 1, 0))
  
  dev.off()
  
  return(cor_matrix)
}
#Arguments for PN
tumor_type <- "Plexiform Neurofibroma"
tissue <- "primary tumor"
gene_names <- c("COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6",
                "CDH1", "TJP1", "DSP", "SNAI1", "SNAI2", "TWIST1", "LEF1",
                "CDH2", "VIM", "CTNNB1", "ACTA2", "TRAP1")
#Function call to PN
cor_matrix <- plotSquareCorrPlot(edgeRlist, tumor_type, tissue, gene_names)

#Arguments for MPNST
tumor_type <- "Malignant Peripheral Nerve Sheath Tumor"
tissue <- "primary tumor"
gene_names <- c("COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6",
                "CDH1", "TJP1", "DSP", "SNAI1", "SNAI2", "TWIST1", "LEF1",
                "CDH2", "VIM", "CTNNB1", "ACTA2", "TRAP1")
#Function call to MPNST
cor_matrix <- plotSquareCorrPlot(edgeRlist, tumor_type, tissue, gene_names)