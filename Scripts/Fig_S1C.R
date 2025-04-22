library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(tidyverse)
library(readxl)
library(stringi)
library(tidyr)
library(openxlsx)
library(fgsea)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(biomaRt)


# load working directory
setwd("set_working_directory")

set.seed(1507)

#First load the gene list
gene_table <- read.xlsx("MAGeCK_output.xlsx") #(See Table S1 for MAGeCK output)
colnames(gene_table)[1]="symbol_id"

#gene_table = gene_table[gene_table$neg_fdr_essentiality > 0.01,] #filter on non-significant genes for essentiality for enriched gene analysis

top_100_genes <- gene_table %>%
  arrange(Score_dox) %>%  # Sort in descending order (ascending if interested in enriched genes)
  slice_head(n = 100)             # Select the top 100 rows

# Convert gene symbols to ENTREZ IDs
top_100_genes$entrez_id <- mapIds(
  org.Hs.eg.db,
  keys = top_100_genes$symbol_id,  # Input gene symbols
  column = "ENTREZID",            # Target ID type
  keytype = "SYMBOL",             # Input ID type
  multiVals = "first"             # Handle multiple matches
)

gene_list <- top_100_genes$entrez_id[!is.na(top_100_genes$entrez_id)]
#gene_list <- top_100_genes$symbol_id[!is.na(top_100_genes$symbol_id)]


# Load gene sets (e.g., MSigDB hallmark sets)
gene_sets <- msigdbr(species = "Homo sapiens", category = "C2") %>%
  dplyr::select(gs_name, entrez_gene)

# Perform enrichment analysis
overlap_result <- enricher(
  gene = gene_list,
  TERM2GENE = gene_sets,
  pvalueCutoff = 0.05)

############################# Calculating normalized gene ratio & plotting #####################################################

#changing "X/Y" to ratio for GeneRatio
overlap_result@result$GeneRatio <- sapply(overlap_result@result$GeneRatio, function(x) {
  # Split the string by "/"
  split_values <- unlist(strsplit(x, "/"))
  
  # Convert to numeric and calculate the ratio
  num1 <- as.numeric(split_values[1])  # Numerator
  num2 <- as.numeric(split_values[2])  # Denominator
  
  # Return the ratio, ensuring no division by zero or NA values
  if (!is.na(num1) && !is.na(num2) && num2 != 0) {
    return(num1 / num2)
  } else {
    return(NA)  # Return NA if any part is invalid
  }
})

#changing "X/Y" to ratio for BgRatio
overlap_result@result$BgRatio <- sapply(overlap_result@result$BgRatio, function(x) {
  # Split the string by "/"
  split_values <- unlist(strsplit(x, "/"))
  
  # Convert to numeric and calculate the ratio
  num1 <- as.numeric(split_values[1])  # Numerator
  num2 <- as.numeric(split_values[2])  # Denominator
  
  # Return the ratio, ensuring no division by zero or NA values
  if (!is.na(num1) && !is.na(num2) && num2 != 0) {
    return(num1 / num2)
  } else {
    return(NA)  # Return NA if any part is invalid
  }
})

overlap_result@result$NormalizedGeneRatio <- overlap_result@result$GeneRatio / overlap_result@result$BgRatio

# Filter for significant pathways (adjusted p-value < 0.05)
filtered_result <- subset(overlap_result@result, p.adjust < 0.05)

top_10 <- filtered_result %>% 
  arrange(p.adjust) %>% 
  slice_head(n = 10)

# Plot the top 10 pathways
ggplot(top_10, aes(x = NormalizedGeneRatio, y = reorder(Description, NormalizedGeneRatio), color = p.adjust)) +
  geom_point(size = 5) +  # Larger dots
  scale_color_gradientn(colors = c("red", "blue"), limits = c(0.01, 0.05)) +
  theme_minimal() +
  scale_y_discrete(
    breaks = top_10$Description,  # Make sure y-axis ticks are from the full top_10 dataset
    expand = expansion(mult = c(0.05, 0.05))  # Keeps spacing consistent
  ) + 
  labs(
    x = "Normalized Gene Ratio",
    y = NULL,  # Removes y-axis label
    title = "Top 10 Significant Pathways (Adjusted p-value < 0.05)",
    color = "Adjusted p-value"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  )
