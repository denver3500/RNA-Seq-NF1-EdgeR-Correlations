library(SummarizedExperiment)
library(dplyr)


files <- c("JH_batch1_salmon.merged.gene_counts.tsv", "JHU_batch2_salmon.merged.gene_counts.tsv", "JHU_batch2_addl_salmon.merged.gene_counts.tsv", "JHU_batch3_salmon.merged.gene_counts.tsv")
annotation <- read.csv("Simple_sample_map_10062023.csv")
merged_data <- read.table(files[1], header = TRUE, sep = "\t")

# Iterate over files and merge them
for (file in files[-1]) {
  data <- read.table(file, header = TRUE, sep = "\t")
  merged_data <- cbind(merged_data, data[,-c(1,2)])
}

#Filter metadata to leave only Bulk RNA sequencing samples
filtered_annotation <- annotation %>%
  filter(assay == "Bulk RNA sequencing") %>%
  mutate(merged_id = paste0(gsub("-", ".", specimenID), ".", aliquotID))

#Filter data to leave only samples present in metadata. In this case, it's all the samples. 
merged_ids <- filtered_annotation$merged_id
columns_to_keep <- c("gene_id", "gene_name", intersect(merged_ids, colnames(merged_data)))
filtered_data <- merged_data[, columns_to_keep, drop = FALSE]

#Building summarize experiment
expression_matrix <- filtered_data
gene_metadata <- expression_matrix %>% select(gene_id, gene_name)
expression_matrix <- expression_matrix %>% select(-gene_name)
rownames(expression_matrix) <- expression_matrix$gene_id
expression_matrix <- expression_matrix %>% select(-gene_id)
samples_metadata <- filtered_annotation %>% select(merged_id, everything())
samples_metadata <- samples_metadata %>% filter(merged_id %in% colnames(expression_matrix))
stopifnot(rownames(expression_matrix) == gene_metadata  $gene_id)
stopifnot(colnames(expression_matrix) == samples_metadata$merged_id)

se <- SummarizedExperiment(assays = list(counts = as.matrix(expression_matrix)), colData = samples_metadata, rowData = gene_metadata)

saveRDS(se, "combined_se.rds")
