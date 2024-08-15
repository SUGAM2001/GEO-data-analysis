# Set Working Directory
setwd("C:/Users/comp_/OneDrive/Documents/IOB/tc_samples")

# Install packages if required
install.packages("dplyr")  

#load libraries
library(dplyr)
library(dplyr)
library(tidyverse)
library(GEOquery)

#Upload data
data <- read.csv("GSE94464_norm_counts_FPKM_GRCh38.p13_NCBI.tsv", sep = "\t")
dim(data)
head(data)

#get metadata
gse<-getGEO(GEO="GSE94464", GSEMatrix=TRUE)
gse

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# data processing
# select desired columns
metadata.subset <- select(metadata,c(1,10,11,17))
metadata.modified<-metadata %>%
  select(1,10,11,17) %>%
# rename columns
  rename(tissue=characteristics_ch1)%>%
  rename(metastasis=characteristics_ch1.1)%>%
# remove characters that are not required
  mutate(tissue=gsub("cell line: ","", tissue))%>%
  mutate(tissue=gsub("Anaplastic thyroid cancer cell line","Anaplastic thyroid cancer", tissue))%>%
  mutate(tissue=gsub("papillary thyroid cancer cell line","Papillary thyroid cancer", tissue))%>%
  mutate(metastasis=gsub("tissue: ","", metastasis))

# reshaping data
data.long<- data %>%
  gather(key="samples",value="FPKM",-GeneID)

# convert metadata index to column so that sample ids can be accessed
metadata.modified$index <- rownames(metadata.modified)

# join data frames = data.long +metadata.modified
data.long <- data.long %>%
  left_join(.,metadata.modified,by=c("samples"="index"))

# explore data
head(data.long)

#getting gene names for gene ids
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)
data.long$GeneID <- as.character(data.long$GeneID)
data.long$GeneName <- mapIds(org.Hs.eg.db, keys = data.long$GeneID, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
head(data.long)

# data analysis
data.long %>%
  filter(GeneName=="BRAF"|GeneName=="TP53") %>%
  group_by(GeneName, tissue) %>%
  summarize(mean_FPKM = mean(FPKM),median_FPKM=median(FPKM))%>%
  arrange(mean_FPKM)

#data visualization

# 1.Bar Plot
ba<-data.long %>%
  filter(GeneName == "BRAF") %>%
  ggplot(.,aes(x=samples,y=FPKM, fill=tissue)) + 
  geom_col() + theme(legend.position = "bottom")+
  ggtitle("Expression of BRAF Gene Across Samples")
ggsave(ba,filename="barplot.png", width=12,height=8)

# 2.Density plot
d<-data.long %>%
  filter(GeneName == "BRAF")%>%
  ggplot(.,aes(x=FPKM,fill=tissue))+
  ggtitle("Distribution of BRAF Gene Expression Across Different Tissues") +
  geom_density(alpha=0.5)
ggsave(d,filename="density.png", width=10,height=8)

# 3. Box Plot
b<-data.long %>%
  filter(GeneName =="BRAF")%>%
  ggplot(.,aes(x=tissue,y=FPKM))+
  geom_boxplot()+ theme(legend.position = "bottom")+
  ggtitle("Variation in BRAF Gene Expression by Tissue Type")
ggsave(b,filename="boxplot.png", width=12,height=8)

# 4. violin plot
v<-data.long %>%
  filter(GeneName =="BRAF")%>%
  ggplot(.,aes(x=tissue,y=FPKM,fill=tissue))+
  geom_violin()+ theme(legend.position = "bottom")+
  ggtitle("Distribution and Density of BRAF Gene Expression Across Tissues")
ggsave(v,filename="violin.png", width=10,height=8)

# remove Gene Id column to avoid getting null values during scatter plotting
data.long <- data.long[, !(names(data.long) %in% "GeneID")]

# 5. Scatter Plot
s<-data.long %>%
  filter(GeneName =="BRAF" | GeneName == "NRAS")%>%
  spread(key=GeneName, value=FPKM)%>%
  ggplot(., aes(x=BRAF, y=NRAS, color=tissue))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE)+
  ggtitle("Correlation Between BRAF and NRAS Gene Expression Across Tissues")
ggsave(s,filename="scatterplot.png", width=10,height=8)

# 6. HeatMap
genes.of.interest <- c("BRAF","TP53","NRAS","TERT","HRAS")
p <- data.long %>%
  filter(GeneName %in% genes.of.interest)%>%
  ggplot(.,aes(x=samples,y=GeneName, fill=FPKM))+
  geom_tile()+
  scale_fill_gradient(low = "white",high="red")+
  ggtitle("Heatmap of BRAF,TP53,NRAS,TERT,HRAS")
ggsave(p,filename="heatmap.png", width=12,height=8)


