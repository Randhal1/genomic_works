# Read the CSV file into a data frame
df <- read.csv("condition_Atrophic_vs_Control.csv")

# Filter for significant genes (adj-p < 0.05)
significant_genes <- df[df$padj < 0.05, ]

# Filter the data frame for genes with padj < 0.05 and count the rows
count_significant_genes <- sum(df$padj < 0.05, na.rm = TRUE)

# Print the result
print(count_significant_genes)

# Count upregulated genes (log2FoldChange > 0)
upregulated_genes <- sum(significant_genes$log2FoldChange > 0, na.rm = TRUE)

# Count downregulated genes (log2FoldChange < 0)
downregulated_genes <- sum(significant_genes$log2FoldChange < 0, na.rm = TRUE)

# Print the results
print(paste("Number of upregulated genes:", upregulated_genes))
print(paste("Number of downregulated genes:", downregulated_genes))

# Count genes with a log2FoldChange > 2
count_log2_gt_2 <- sum(significant_genes$log2FoldChange > 2, na.rm = TRUE)

# Count genes with a log2FoldChange < -2
count_log2_lt_neg_2 <- sum(significant_genes$log2FoldChange < -2, na.rm = TRUE)

# Print the results
print(paste("Number of genes with log2FoldChange > 2:", count_log2_gt_2))
print(paste("Number of genes with log2FoldChange < -2:", count_log2_lt_neg_2))

