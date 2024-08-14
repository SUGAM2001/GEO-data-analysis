#script to perform differential gene expression amalysis using DESeq2 packages

setwd("C:/Users/comp_/OneDrive/Documents/IOB/tc_samples")

# load libraries
library(DESeq2)
library(tidyverse)
library(airway)

# Step 1: Preparing count data

## read in counts data
counts_data <- read.csv("GSE94464_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t")

# # Install and load GEOquery
if (!require(GEOquery)) {
  install.packages("GEOquery")
}
library(GEOquery)

# Download the dataset
gse <- getGEO("GSE94464", GSEMatrix = TRUE)

# Extract the first ExpressionSet object
eset <- gse[[1]]

# Extract sample information
sample_info <- pData(eset)

# View the sample information
head(sample_info)

write.csv(sample_info, file = "GSE94464_sample_info.csv")

# select desired columns
sample_info.subset <- select(sample_info,c(2,10))
sample_info.modified<-sample_info %>%
  select(2,10) %>%
  # rename columns
  rename(cell_line = characteristics_ch1)%>%
  # remove characters that are not required
  mutate(cell_line=gsub("cell line: ","", cell_line))%>%
  mutate(cell_line=gsub("Anaplastic thyroid cancer cell line","Anaplastic thyroid cancer", cell_line))%>%
  mutate(cell_line=gsub("papillary thyroid cancer cell line","Papillary thyroid cancer", cell_line))

col_data <- sample_info.modified

# Set GeneID column as index
rownames(counts_data) <- counts_data$GeneID

# Remove GeneID column
counts_data <- counts_data[, -which(names(counts_data) == "GeneID")]

# Make sure the row names in col_data matches col_names in count_data

all(colnames(counts_data)%in% rownames(col_data))

# Are they in same order
all(colnames(counts_data) == rownames(col_data))

# Construct a DESeqDataSet object -----------
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = col_data, design= ~ cell_line)

# pre-filtering: removing rows with low gene counts
# keeping rows: that have at least 10 reads total
keep= rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# set the factor level
dds$cell_line <- relevel(dds$cell_line, ref = "Papillary thyroid cancer")

dds$cell_line
dds <- DESeq(dds)


# Explore Results----

res <- results(dds)
res

# Explore Results
summary(res)

# reset adjusted p-value cutoff
res0.01 <- results(dds, alpha=0.01)
summary(res0.01)

resultsNames(dds)

# Visualisation
# MA plot
plotMA(res)

volcano_data <- as.data.frame(res)
volcano_data$GeneID <- rownames(volcano_data)
#getting gene names for gene ids
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)
volcano_data$GeneID <- as.character(volcano_data$GeneID)
volcano_data$GeneName <- mapIds(org.Hs.eg.db, keys = volcano_data$GeneID, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
head(volcano_data)

vol_plot <- volcano_data %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) + 
  geom_point() 

# Visualise simple volcano plot
vol_plot + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")


volcano_data<-volcano_data %>%
  mutate(gene_type = case_when(log2FoldChange >= 2 & padj <= 0.05 ~ "up",
                               log2FoldChange <= 0.5 & padj <= 0.05 ~ "down",
                               TRUE ~ "ns"))   
# Load necessary libraries
library(ggplot2)
library(ggrepel)
# Define criteria for top regulated genes
top_upregulated_genes <- volcano_data %>%
  dplyr::filter(padj < 0.05, log2FoldChange > 0) %>%  # Filter for significance and upregulation
  dplyr::arrange(desc(log2FoldChange)) %>%  # Sort by fold change in descending order
  dplyr::slice_head(n = 10)  # Select top 10 upregulated genes

downregulated_genes <- volcano_data %>%
  dplyr::filter(padj < 0.05, log2FoldChange < 0) %>%  # Filter for significance and downregulation
  dplyr::arrange(desc(abs(log2FoldChange))) %>%  # Sort by absolute fold change
  dplyr::slice_head(n = 10)  # Select top 10 downregulated genes

# Obtain gene_type counts ------------------------------------------------------           
volcano_data %>%
  count(gene_type)
# Check gene_type categories ---------------------------------------------------
volcano_data %>%
  distinct(gene_type)%>%
  pull()

# Add colour, size and alpha (transparency) to volcano plot -
cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

volcano_data %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  geom_label_repel(data = top_upregulated_genes,  
                   aes(label = GeneName),
                   force = 2,
                   nudge_y = 1)+
  geom_label_repel(data = downregulated_genes,  
                   aes(label = GeneName),
                   force = 2,
                   nudge_y = 1)+
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-10, 10, 2)),       
                     limits = c(-10, 10)) 
