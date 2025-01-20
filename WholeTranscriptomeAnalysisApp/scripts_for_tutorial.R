#####RNA-Seq data analysis and functional enrichment analysis 
#### I have used DESEq2 package for normalization
#### Downstream functional analysis was done using Pathview and ReactomePA packages
getwd()
setwd("/Users/aigulsharip/Library/CloudStorage/GoogleDrive-Aigul.Sharip@nu.edu.kz/My Drive/Bioinformatics/R/RShinyProjects/RShinyProjects/WholeTranscriptomeAnalysisApp/")

############################################
#Installation and uploading necessary libraries
################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("biomaRt")
install.packages("pheatmap")
BiocManager::install("airway")
BiocManager::install("pathview")
BiocManager::install("gage")
BiocManager::install("gageData")

#Load packages
library(DESeq2)
library(readxl)
library(biomaRt)
library("pheatmap")
library(readxl)

#This file contains reads per gene
xlsx = "count_matrix.xlsx"
df = as.data.frame(read_excel(xlsx))
head(df)


############################################################
# prepare dataframe
############################################################
colnames(df)[1] <- "genes"  
row.names(df) <- df$genes   
df <- df[, -1]




tmp = colnames(df)
tmp[1] = "genes"
colnames(df) = tmp
row.names(df) = df[,1]
df = df[, tmp[2:length(tmp)]]

############################################################
# prepare colData design matrix
############################################################


samples <- colnames(df)
libType <- rep("paired-end", length(samples))
condition <- c(rep("normal", 11), rep("tumor", 22))
tnm <- c(rep("normal", 11), rep(c("T3", "T3", "T3", "T3", "T3", "T3", "T3", "T3", "T3", "T4", "T2", "T3", "T1", "T3", "T3", "T2", "T3", "T3", "T3", "T4", "T3", "T3"), each=1))
replica <- ifelse(seq_along(samples) <= 11, "normal", samples[12:length(samples)])

colData <- data.frame(samples, condition, libType, tnm, replica, row.names = samples)
View(colData)


############################################################
# DESeq Tumor vs Normal
############################################################
filtered = log10(rowSums(df))>3
counts = df[filtered,]
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds); 
dds
res <-results(dds)
res
#write.csv(as.data.frame(res), file=results_tumor_vs_normal.csv")
summary(res)
res <- results(dds, alpha=.05, lfcThreshold = 1)
table(res$padj < .05)




#plotting the results
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("condition"))
plotMA(res, ylim=c(-2,2), main="DESeq2")

#Gene clustering
#install.packages("pheatmap")
library(pheatmap)
rld <- rlog(dds)
head(assay(rld), 3)

select <- order(rowVars(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)


#ReactomePA

BiocManager::install("ReactomePA")
library(ReactomePA)
geneList<- res$entrez
head(geneList)
x<-enrichPathway(gene=geneList, pvalueCutoff = 0.05, readable = TRUE)
head(as.data.frame(x))
barplot(x, showCategory = 8)
dotplot(x, showCategory = 8)
emapplot(x)
cnetplot(x, categorySize="pvalue", foldChange = geneList, showCategory = 5)
















