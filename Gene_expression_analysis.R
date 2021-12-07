## Differential expression analysis (R 4.0.4).

## Libraries used
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
library(edgeR)
library(goseq)
library(enrichplot)
library(WGCNA)
library(MDSeq)

## Import data resulting from Salmon into R.
samples <- list.files(path = "./data", full.names = T)
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "./data/", "")
tx2gene <- read.delim("annotation_file.txt") # Load genome annotation file.
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[,c("gene", "genefull")], countsFromAbundance = "lengthScaledTPM")
data <- txi$counts %>% 
  round() %>% 
  data.frame() # Write counts to object.

## Create metadata.
sampletype <- factor(c("Source","Core","Intermediate",…))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))

## Filtering.
data.unfiltered <- DGEList(counts = data, group = sampletype)
keep <- filterByExpr(data.unfiltered, min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7)
data.prefiltered <- data.unfiltered[keep,,keep.lib.sizes=FALSE]
data.prefiltered <- data.prefiltered$counts %>% 
  data.frame()
keep2 <- apply(data.prefiltered, 1, function(x) sum(x >= 10) >= 10)
data.filtered <- data.prefiltered[keep2, ]

## Differential expression analysis.
dds <- DESeqDataSetFromMatrix(countData = data.filtered, colData = meta, design = ~ sampletype)
dds <- DESeq(dds) # Wald test.

## Define contrasts for Core vs Source. Repeat for all other pairwise comparisons.
contrast_C_S <- c("sampletype", "Core", "Source")
res_table_C_S <- results(dds, contrast = contrast_C_S, alpha = 0.05)
write.table(res_table_C_S, file = "res_table_C_S.txt", sep = "\t", quote = F, col.names = NA)
padj.cutoff <- 0.05 #Set thresholds.
res_table_C_S_tb <- res_table_C_S %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>% 
  as_tibble() #Create tibble.
sig_C_S <- res_table_C_S_tb %>%
  filter(padj < padj.cutoff) # Keep only significant genes.
write.table(sig_C_S, file = "sig_C_S.txt", sep = "\t", quote = F, col.names = NA)

## Data representation.
normalized_counts <- counts(dds, normalized = T) %>% 
  data.frame() %>%
  rownames_to_column(var = "gene")
gene_names <- tx2gene %>% 
  dplyr::select(gene, genefull) %>% 
  dplyr::distinct()
colnames(gene_names) <- c("genefull", "gene")
normalized_counts_genes <- merge(normalized_counts, gene_names, by.x = "gene")
normalized_counts_genes <- normalized_counts_genes %>%
  as_tibble()
norm_C_S_sig <- normalized_counts_genes[,c(22,24:29,36:53)] %>% 
  filter(normalized_counts_genes$gene %in% sig_C_S$gene) # Extract normalised expression for significant genes from the Core and Source samples. Repeat for all other pairwise comparisons.
heat_colors <- brewer.pal(6, "RdYlBu") # Set colour palette.
pheatmap(norm_C_S_sig[1:25], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = meta, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20) # Plot heatmap of DEGs.
res_table_C_S_tb <- res_table_C_S_tb %>% 
  mutate(threshold_ = padj < 0.05 & abs(log2FoldChange) >= (log2(1)))
ggplot(res_table_C_S_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_)) +
  ggtitle("Core v Source") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) # Volcano plot.

## Functional analysis.
gene_len <- read.table("/path/gene_length.txt", header = TRUE) # Load data frame with transcript length of each gene.
len <- gene_len[order(gene_len$gene),]
gene_cat <- read.table("/path/gene_cat.txt", header = TRUE) # Load data frame with mapping between genes and GO categories.
gene_cat$GO <- substr(gene_cat$GO, 0, 10)
ALL <- res_table_C_S$gene
DEG <- sig_C_S$gene
DEG.vector <- c(t(DEG))
ALL.vector <- c(t(ALL))
gene.vector <- as.integer(ALL.vector %in% DEG.vector)
names(gene.vector) <- ALL.vector
gene.vector <- sort(gene.vector)
len_ALL <- subset(len, len$gene %in% ALL)
len.vector <- len_ALL$length
pwf <- nullp(gene.vector, bias.data = len.vector, plot.fit = TRUE)
GO <- goseq(pwf, gene2cat = gene_cat, method = "Wallenius")
GO$padj <- p.adjust(GO$over_represented_pvalue, method = "fdr")
resGO <- subset(GO, GO$padj < 0.05)

## Likelihood ratio test analysis.
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
res_LRT <- results(dds_lrt)
res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>% 
  as_tibble()
sig_LRT <- res_LRT_tb %>% 
  filter(padj < padj.cutoff)
clustering_sig_genes <- sig_LRT %>%
  arrange(padj)
cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
clusters <- degPatterns(cluster_rlog, metadata = meta, time = "sampletype", col = NULL) # Show gene clusters across sample groups.

## Differential variability analysis.
dat.normalised <- normalize.counts(data.filtered, group = sampletype, method = "TMM")
groups <- factor(sampletype, labels = c("Source", "Australia"))
groups_contrasts <- get.model.matrix(groups, contrast.type = "contr.treatment")
fit <- MDSeq(dat.normalised, contrast = groups_contrasts, mc.cores = 4)
results_S_A <- extract.ZIMD(fit, compare = list(A = "Source", B = "Australia"))
padj.cutoff <- 0.05
results_S_A_tb <- results_S_A %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>% 
  as_tibble()
sig_S_A <- results_S_A_tb %>%
  filter(FDR.dispersion < padj.cutoff)


