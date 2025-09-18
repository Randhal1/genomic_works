#Load into R
library(tximport)

samples <- read.table("rsem_results/samples.txt", header = TRUE)
samples

files <- file.path("rsem_results", paste0(samples$sample, ".genes.results"))
names(files) <- samples$sample
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)
write.csv(as.data.frame(txi.rsem$counts), file="raw_counts.csv")

zero_length_and_unexpressed = (apply(txi.rsem$counts, 1, max) > 0) &
                              (apply(txi.rsem$length, 1, min) > 0)

txi.rsem$length = txi.rsem$length[zero_length_and_unexpressed,]
txi.rsem$abundance = txi.rsem$abundance[zero_length_and_unexpressed,]
txi.rsem$counts = txi.rsem$counts[zero_length_and_unexpressed,]


library('DESeq2')
sampleTable <- data.frame(condition = samples$pop)
sampleTable$condition <- relevel(factor(sampleTable$condition), ref = "ctl")

rownames(sampleTable) <- colnames(txi.rsem$counts)
dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)

dds <- DESeq(dds)
resultsNames(dds)

baseMeanPerLvl <- sapply( levels(dds$condition), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$condition== lvl] ) )
write.csv(baseMeanPerLvl, file="normalized_counts.csv")

library(tidyverse)
trans_cts_corr <-counts(dds,normalized=TRUE) %>% cor(method = "spearman")
  # we use Spearman's correlation, a non-parametric metric based on ranks

trans_cts_corr[1:5, 1:5]
library(corrr)

rplot(trans_cts_corr) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
write.csv(trans_cts_corr, file="correlation.csv")


res<-results(dds,name="condition_atr_vs_ctl")
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file="condition_Atrophic_vs_Control.csv")

