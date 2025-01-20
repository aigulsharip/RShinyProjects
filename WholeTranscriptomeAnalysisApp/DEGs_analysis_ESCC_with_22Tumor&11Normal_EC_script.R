#####RNA-Seq data analysis and functional enrichment analysis of ESCC patients
#### I have used DESEq2 package for normalization
#### Downstream functional analysis was done using Pathview and ReactomePA packages

############################################
#Installation and uploading necessary libraries
################################################

# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
BiocManager::install("DESeq2")    # for DESeq2
BiocManager::install("biomaRt")   # for repalce gene names
install.packages("readxl")        # for reading data from excel file 
install.packages("pheatmap")      # for heatmaps
#BiocManager::install("biomaRt", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("biomaRt")

install.packages("pheatmap")

BiocManager::install("airway")
BiocManager::install("pathview")
BiocManager::install("gage")
BiocManager::install("gageData")

library(DESeq2)
library(readxl)
library(biomaRt)
library("pheatmap")

#This file contains reads per gene for 22 ESCC tumor tissue and 11 normal esophegeal tissue
xlsx = "Merge_STARReadsPerGene_22Tumor_Samples&11Normal_Tissues.xlsx"
df = as.data.frame(read_excel(xlsx))
head(df)


############################################################
# prepare dataframe
############################################################
tmp = colnames(df)
tmp
tmp[1] = "genes"
colnames(df) = tmp
row.names(df) = df[,1]
df = df[, tmp[2:length(tmp)]]
head(df)
############################################################
# prepare colData design matrix
############################################################
samples = colnames(df)
head(samples)
libType = rep("paired-end", length(samples)); libType
condition = c(rep("normal", times=11), rep("tumor", times=22))
condition
tnm = c(rep("normal", times =11), "T3", "T3", "T3", "T3", "T3", "T3", "T3", "T3", "T3", "T4", "T2", "T3", "T1", "T3", "T3", "T2", "T3", "T3", "T3", "T4", "T3", "T3")
tnm
replica = colnames(df); 
replica[1:11] = rep("normal", times=11)
colData = cbind(samples, condition, libType, tnm, replica)
rownames(colData) = colnames(df)
View(colData)
############################################################
# DESeq Tumor vs Normal
############################################################
project_name = "DESeq_CONDITION_filtered"
filtered = log10(rowSums(df))>3
counts = df[filtered,]
dds <- DESeqDataSetFromMatrix(counts, colData, formula(~ condition)); 
dds <- DESeq(dds); resultsNames(dds)
resTumorvsNormal <-results(dds)
write.csv(as.data.frame(resTumorvsNormal), file="DESeq_filtered_all_22TUMOR_vs_11NORMAL.CSV")
summary(resTumorvsNormal)

#TUMOR_vs_NORMAL_adding gene names
resTumorvsNormal$ensembl <- sapply(strsplit(rownames(resTumorvsNormal),split="\\+"),"[",1)
ensembl=useMart("ensembl",dataset = "hsapiens_gene_ensembl")
write.csv(as.data.frame(resTumorvsNormal), file="DESeq__filtered_all_22TUMOR_vs_11NORMAL_withgenename.CSV")

#selecting the the most DEG by p-adjusted <0.05 and ordering
resTumorvsNormal<-resTumorvsNormal[which(resTumorvsNormal$padj<0.05),]
resTumorvsNormal <- resTumorvsNormal[order(resTumorvsNormal$log2FoldChange),]
summary(resTumorvsNormal)
#We can also filter based on log2FoldChange
resTumorvsNormal<-resTumorvsNormal[abs(resTumorvsNormal$log2FoldChange) >= 1.0,]
write.csv(as.data.frame(resTumorvsNormal), file="DESeq_22TUMOR_vs_11NORMAL_significant_threshold_padj&log2fc_withgenename.csv")
summary(resTumorvsNormal)

######sort based on p-adj
resTumorvsNormal <- resTumorvsNormal[order(resTumorvsNormal$padj),]
summary(resTumorvsNormal)
write.csv(as.data.frame(resTumorvsNormal), file="DESeq_22TUMOR_vs_11NORMAL_significant_threshold_padj&log2fc_withgenename_sorted_padj.csv")



#plotting the results
topGene <- rownames(resTumorvsNormal)[which.min(resTumorvsNormal$padj)]
plotMA(resTumorvsNormal, ylim=c(-13,13), main="DESeq2")


############################################################
# We can find DEGs for each tumor stage (T1-T4)
# DESeq by tnm group
# create results tables for each TNM groups_filtered
counts = df
filtered = log10(rowSums(df))>3 ; counts = df[filtered,]
dds <- DESeqDataSetFromMatrix(counts, colData, formula(~ tnm)); dds <- DESeq(dds); resultsNames(dds)
resT1_N <- results(dds, contrast = c("tnm", "T1","000"))
resT2_N <- results(dds, contrast = c("tnm", "T2","000"))
resT3_N <- results(dds, contrast = c("tnm", "T3","000"))
resT4_N <- results(dds, contrast = c("tnm", "T4","000"))

resT1_N<-resT1_N[which(resT1_N$padj<0.05),]
resT1_N <- resT1_N[order(resT1_N$log2FoldChange),]
resT1_N<-resT1_N[abs(resT1_N$log2FoldChange) >= 1.0,]

resT2_N<-resT2_N[which(resT2_N$padj<0.05),]
resT2_N <- resT2_N[order(resT2_N$log2FoldChange),]
resT2_N<-resT2_N[abs(resT2_N$log2FoldChange) >= 1.0,]


resT3_N<-resT3_N[which(resT3_N$padj<0.05),]
resT3_N <- resT3_N[order(resT3_N$log2FoldChange),]
resT3_N<-resT3_N[abs(resT3_N$log2FoldChange) >= 1.0,]

resT4_N<-resT4_N[which(resT4_N$padj<0.05),]
resT4_N <- resT4_N[order(resT4_N$log2FoldChange),]
resT4_N<-resT4_N[abs(resT4_N$log2FoldChange) >= 1.0,]

#exporting DEGS to file
write.csv(as.data.frame(resT1_N), file="DESeq_significant_threshold_padj&log2fc_T1_vs_NORMAL.csv")
write.csv(as.data.frame(resT2_N), file="DESeq_significant_threshold_padj&log2fc_T2_vs_NORMAL.csv")
write.csv(as.data.frame(resT3_N), file="DESeq_significant_threshold_padj&log2fc_T3_vs_NORMAL.csv")
write.csv(as.data.frame(resT4_N), file="DESeq_significant_threshold_padj&log2fc_T4_vs_NORMAL.csv")

#summary stats
summary(resTumorvsNormal)



#######Sample-sample distances
rld <- rlog(dds, blind=FALSE)
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$tnm, rld$samples, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#Heatmap on individual samples for all DEGs
select <- order(rowVars(counts(dds,normalized=TRUE)),
                decreasing=TRUE)
pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE)
mat1 = assay(rld)[select,]
q050 = add_genes_and_codes(mat1); pheatmap(q050) 
write.csv(as.data.frame(mat1), file="heatmap_matrix_individual_samples_all_DEGs.csv")
write.csv(as.data.frame(q050), file="heatmap_matrix_individual_samples_all_DEGs_q050.csv")


#Heatmap on individual samples for top 100 genes
select <- order(rowVars(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE)
mat1 = assay(rld)[select,]
q050 = add_genes_and_codes(mat1); pheatmap(q050) 
write.csv(as.data.frame(mat1), file="heatmap_matrix_individual_samples_mat1.csv")

write.csv(as.data.frame(q050), file="heatmap_matrix_individual_samples_q050.csv")
colnames(q050) <- c("N31","N32","N33","N34","N35","N36","N90","N97","N54","N55","N56","001T","002T","003T","004T","005T","007T","008T","009T","010T","011T","012T","013T","014T","015T","016T","017T","019T","021T","023T","025T","029T","033T")
## changes the genes names
rownames(q050) <-c("MT-CO1","IGHG1","MT-RNR2","MT-ND4","MT-CO3","RN7SL2","MALAT1","NEAT1","KRT4","KRT13","RNA45S5","MT-CO2","MT-CYB","RN7SL1","MT-ATP6","TTN","7SK","MT-ND2","S100A8","MT-ND1","DES","SPRR3","MYH11","MT-ND5","RF00009","AHNAK","ACTG2","ACTB","KRT5","COL3A1","KRT6A","ANXA1","MT-ND3","FTL","EEF1A1","DSP","TPM2","CSTB","COL1A2","IGHG4","KRT16","IGHA1","RPL37","COL1A1","S100A9","EMP1","H19","FN1","CRNN","KRT14","B2M","MT-RNR1","MT-TY","RMRP","S100A7","IGHG3","TAGLN","MMP7","S100A2","CSTA","MUC5B","SPINK5","MT-ND6","IGHG2","MYL6","MUC21","PIGR","S100A6","HLA-B","RPS27","RPLP0","ACTG1","SPARC","RPS18","PPL","TMSB4X","UBC","FLNA","RHCG","CSRP1","AQP3","MT-ND4L","PPIA","RPLP1","NEB","MTATP6P1","RPL10","KRT6B","RPL7","RPS2","MYLK","MMP1","TGM3","TPT1","RPL9","KRT19","LCN2","MAL","IL1RN","PERP")


#Heatmap on individual samples, version 2 
select <- order(rowVars(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:50]
df_tnm <- as.data.frame(colData[,c("tnm")])
mat_ind = assay(rld)[select,]
pheatmap(mat_ind)
pheatmap(mat_ind, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col = df_tnm)
q_ind = add_genes_and_codes(mat_ind); 
pheatmap(q_ind, annotation_col = df_tnm)
View(mat_ind)
View(q_ind)

colnames(q_ind) <- c("T1","T2", "T3","T4")
rownames(q_ind) <- c("MT-CO1","MT-CO3","IGHG1","MT-RNR2","MT-ND4","MT-CO2","GC00P7E0035","RN7SL2 ","MALAT1","NEAT1","KRT4","KRT13","MT-CYB","MT-ND2","MT-ATP6","MT-ND1","RN7SL1 ","TTN","RN7SK","S100A8","SPRR3","DES","MYH11","MT-ND5","PRKAB2","KRT5","AHNAK","FTL","ACTG2","ACTB","KRT14","MT-RNR1","KRT6A","MT-ND3","COL3A1","EEF1A1","ANXA1","IGHG4","DSP","COL1A1","TPM2","CSTB","COL1A2","S100A9","KRT16","IGHG3","RPL37","EMP1","IGHA1","S100A2","CRNN","FN1","H19","B2M","RMRP","MT-TY","IGHG2","S100A7","TAGLN","CSTA","S100A6","SPINK5","MT-ND6","MUC5B","MMP7","MTATP6P1","HLA-B","MYL6","MUC21","ACTG1","RPS18","RPLP0","PIGR","SPARC","RPS27","MT-ND4L","PPL","FLNA","UBC","PPIA","CD74","TMSB4X","RPLP1","RHCG","S100A11","RPL7","HLA-A","AQP3","CSRP1","RPL8","RPL19","RPS2","NEB","KRT8","KRT6B","TGM3","RPS20","RPL9","RPL10","RPL13A")
pheatmap(q_ind)
pheatmap(q_ind, annotation_col = df_tnm)
View(q_ind)



############################################################
# Contructing heatmaps
############################################################
#Function to add genes names to heatmap 
add_genes = function(res) {
  list_ensembl = sapply(strsplit(rownames(res),split="\\+"),"[",1)
  ensembl = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
  genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = list_ensembl, mart = ensembl )
  idx <- match(list_ensembl, genemap$ensembl_gene_id )
  list_entrez <- genemap$entrezgene[ idx ]
  list_hgnc_symbol <- genemap$hgnc_symbol[ idx ]
  return (list_hgnc_symbol)
}
add_genes_and_codes = function(m) {
  q = add_genes(m)
  w = row.names(m)
  q[is.na(q)] = w[is.na(q)]
  q[q == ""] = w[q == ""]
  row.names(m) = q
  return (m)
}

# Option 1: Top DEGs
betas <- coef(dds)
colnames(betas)
topGenes <- head(order(resTumorvsNormal$padj),50)
mat <- betas[topGenes, -c(1,1)]
pheatmap(mat)
View(mat)
write.csv(as.data.frame(mat), file="heatmap_matrix_res.csv")
#write.csv(as.data.frame(topGenes), file="topGenes_50_real.csv")
View(topGenes)
pheatmap(mat,display_numbers = TRUE)
mat050 = head(mat, 050); pheatmap(mat050)   
q050 = add_genes_and_codes(mat); pheatmap(q050) 
View(q050)

setwd("G:/My Drive/LBSB_Work/Kaz_ESCC/Bioinformatic analysis/KAZ_ESCC_DeSeq_22Tumor&11Normal/Outputs/")
write.csv(as.data.frame(q050), file="heatmap genes with genes for modification.csv")

###Gene names

genes=c("WASH7P","RP11-34P13.9","MTND1P23","MTND2P28","MTCO1P12","MTCO2P12","MTATP6P1","MTCO3P12","RP11-206L10.6","LINC01409","LINC00115","LINC01128","LINC02593","SAMD11","NOC2L","KLHL17","PLEKHN1","PERM1","AL645608.7","HES4","ISG15","AGRN","AL390719.47","RNF223","C1orf159","RP11-465B22.8","TNFRSF18","TNFRSF4","SDF4","B3GALT6","C1QTNF12","UBE2J2","SCNN1D","ACAP3","PUSL1","INTS11","RP5-890O3.9","CPTP","TAS1R3","DVL1","MXRA8","AURKAIP1","CCNL2","MRPL20-AS1","MRPL20","MRPL20-DT","ANKRD65","LINC01770","VWA1","ATAD3C")
colnames(q050) <- c("I","II", "III","IV")
rownames(q050) = genes
View(q050)
pheatmap(q050)


#Option 2: rowVars(assay(rld)
rld <- rlog(dds)
head(assay(rld), 3)
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),50)
mat <- betas[topVarGenes, -c(1,1)]
mat = mat[complete.cases(mat),]
pheatmap(mat, display_numbers = TRUE)
pheatmap(mat)
q050_row = add_genes_and_codes(mat); pheatmap(q050_row) 
colnames(q050_row) <- c("T1","T2", "T3","T4")
## changes the genes names
rownames(q050_row) <-c("AL772337.2","REG4","MTND1P23","MYH1","MT-TY","MIR3687-2","MAGEC1","PPDPFL","AL355075.4 ","RMRP","SCARNA10","CREB3L3","MIR3648-2","AL162581.1 ","GSTM1","RN7SKP227","AL356488.2 ","RN7SKP9","MAGEA4","MYOG","MYL1","RN7SL2","HIST1H4F","RNU1-11P","SNORA54","MUC5B","XIST","MUC21","RN7SL38P","RPS4Y1","MYL2","MYH7","ACTA1","CDR1","MMP1","CHRND","BPIFB2","PGC","HIST1H3B ","SCARNA7","SNORA49","RN7SL1","SNORA37","IGHG4","RN7SL4P","IGHGP","PIGR","HIST1H3I","NRAP","RN7SL3 ","ANKRD20A11P","KRT4","SCARNA5","AL355336.1 ","MAGEA3","H4C3","CRISP3","CRNN","LINC01287","RNU6ATAC","SPINK1","MMP13","H1-5","BPIFB1","AMTN","H1-1","MAGEA10","SCARNA13","TERC","MYH2","KDM5D","S100A7","DDX3Y","NTS","C5orf66-AS2","MAGEA6","MT-TC","MUCL3","MAGEA1","CAV3","RN7SL128P","OLFM4","H1-4","ANKRD20A5P","MMP10","COL18A1-AS1","DMBT1","HIST1H4B ","LINC00626","MMP3","SCGB3A1","MYBPC1","AC012613.2 ","SPP1","TXLNGY","CST1","HIST1H4D ","XIRP2","COL11A1","HIST1H3F ")
View(q050_row)
pheatmap(q050_row)


#### Heatmap for top 100 genes
# Option 1: Top DEGs
betas <- coef(dds)
colnames(betas)
topGenes <- head(order(resTumorvsNormal$padj),100)
mat <- betas[topGenes, -c(1,1)]
pheatmap(mat)
View(mat)
write.csv(as.data.frame(mat), file="heatmap_matrix_res_top100.csv")

 
View(topGenes)
pheatmap(mat,display_numbers = TRUE)
mat050 = head(mat, 050); pheatmap(mat050)   
q050 = add_genes_and_codes(mat); pheatmap(q050) 
View(q050)

### Reading top100 Genes from CSV file
top100Genes = read.csv("Top100_DEGs_for_heatmap.csv", header = T, sep = ",", row.names = 1)
pheatmap(top100Genes)
head(top100Genes)

colnames(top100Genes) <- c("I","II", "III","IV")
pheatmap(top100Genes)

#####STRING PIP analysis
top500Genes <- resTumorvsNormal[!is.na(resTumorvsNormal$entrez), ]
summary(top500Genes)
top500Genes <- top500Genes[order(top500Genes$padj),]
top500Genes <- head (top500Genes,500)

top500Genes <- top500Genes[order(top500Genes$log2FoldChange),]

write.csv(as.data.frame(top500Genes), file="TopGenes_ordered_TumorvsNormal_500.csv")






############################################################
# Functional enrichment analysis
############################################################
#Pathview package for KEGG pathways

#install.packages("Pathview")
#BiocManager::install("pathview")

library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs)
summary(resTumorvsNormal)
foldchanges = resTumorvsNormal$log2FoldChange
names(foldchanges) =resTumorvsNormal$entrez
head(foldchanges, 10)
head(rownames(foldchanges))


keggres<-gage(foldchanges, gsets =kegg.sets.hs, same.dir = TRUE)
lapply(keggres,head, 15)
lapply(keggres[1:3], head)
head(keggres$greater)
head(keggres$less, 30)
write.csv(as.data.frame(lapply(keggres,head, 15)), file="Kegg_pathways.CSV")
write.csv(as.data.frame(lapply(keggres)), file="Kegg_pathways.CSV")

write.csv(as.data.frame(keggres), file="Kegg_pathways.CSV")

library(org.Hs.eg.db)
mapped <- mappedkeys(org.Hs.egPATH2EG)
L <- as.list(org.Hs.egPATH2EG[mapped])
Kegg_ID <- names(L)
Gene_IDs <- sapply(L, paste, collapse=",")
write.table(cbind(Kegg_ID, Gene_IDs), file="KEGG to Genes.txt", sep="\t", row.names=FALSE, col.names=FALSE)


library(org.Hs.eg.db)
mapped <- mappedkeys(org.Hs.egPATH2EG)
L <- as.list(org.Hs.egPATH2EG[mapped])
Kegg_ID <- names(L)
Gene_IDs <- sapply(L, paste, collapse=",")
write.table(cbind(Kegg_ID, Gene_IDs), file="KEGG to Genes.csv", sep="\t", row.names=FALSE, col.names=FALSE)
#Changing the directory
getwd()
setwd("G:/My Drive/LBSB_Work/Kaz_ESCC/Bioinformatic analysis/KAZ_ESCC_DeSeq_22Tumor&11Normal/Outputs/KEGG")
# Drawing pathview diagrams 
install.packages("dplyr")
library(dplyr)
#Pathview for upregulated ($greater)pathways
remove.packages("pathview")
devtools::install_github("javadnoorb/pathview")

keggrespathways <-data.frame(id = rownames(keggres$greater), keggres$greater) %>%
tbl_df() %>%
filter(row_number()<=8) %>%
.$id %>%
as.character()
keggrespathways
keggresids = substr(keggrespathways,start= 1, stop =8)
keggresids
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))
rm(korg)


require(clusterProfiler)
require(pathview)
data(geneList)
hsa04110 <- pathview(gene.data = geneList, 
                       pathway.id = "hsa04110",  species    = "hsa", 
                       limit = list(gene=max(abs(geneList)), cpd=1))

#Pathview for downregulated ($less)pathways
keggrespathways <-data.frame(id = rownames(keggres$less), keggres$less) %>%
  tbl_df() %>%
  filter(row_number()<=2) %>%
  .$id %>%
  as.character()
keggrespathways
keggresids = substr(keggrespathways,start= 1, stop =8)
keggresids
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))

############################################################
# Extract the list of top genes
############################################################
getwd()
setwd("G:/My Drive/LBSB_Work/Kaz_ESCC/Bioinformatic analysis/KAZ_ESCC_DeSeq_22Tumor&11Normal/Outputs/Reactome_GO_ggplots")
#TUMORvsNORMAL
#just ordered by logfdc
topGenes_TumorvsNormal_50 <- head (resTumorvsNormal,50)
write.csv(as.data.frame(topGenes_TumorvsNormal_50), file="TopGenes_ordered_TumorvsNormal_50.csv")

topGenes_TumorvsNormal_100 <- head (resTumorvsNormal,100)
write.csv(as.data.frame(topGenes_TumorvsNormal_100), file="TopGenes_ordered_TumorvsNormal_100.csv")

signDownDEGS <-resTumorvsNormal[which(resTumorvsNormal$log2FoldChange<1),]
write.csv(as.data.frame(signDownDEGS), file="sign_Down_DEGS.csv")

signUPDEGS <-resTumorvsNormal[which(resTumorvsNormal$log2FoldChange>1),]
write.csv(as.data.frame(signUPDEGS), file="sign_Up_DEGS.csv")

topDownDEGS100 = head (signDownDEGS,100)
topUpDEGS100 = head (signUPDEGS,100)
write.csv(as.data.frame(signDownDEGS), file="sign_Down_DEGS_top100.csv")
write.csv(as.data.frame(topUpDEGS100), file="sign_Up_DEGS_top100.csv")




#################
#ReactomePA package for REACTOME packages
#reactome analysis for downregulated tumorvsnormal genesreactome pathways
#BiocManager::install("ReactomePA")

library(ReactomePA)
setwd("G:/My Drive/LBSB_Work/Kaz_ESCC/Bioinformatic analysis/KAZ_ESCC_DeSeq_22Tumor&11Normal/Outputs/Reactome_GO_ggplots")

dataset = "DAVID_TumorVsNormal.xlsx"

## Reactome PA plots for all top 100 DEGS genes
reactome_top100_all = as.data.frame(read_excel(dataset, sheet = "TopGenes_100_All"))
head(reactome_top100_all)
geneList_down<-reactome_top100_all$entrez
head(geneList_down)
x_down<-enrichPathway(gene=geneList_down, pvalueCutoff = 0.05, readable = TRUE)
head(as.data.frame(x_down))
barplot(x_down, showCategory = 10)
dotplot(x_down, showCategory = 10)
BiocManager::install("enrichplot")
library(enrichplot)


x2 <- pairwise_termsim(x_down) 
emapplot(x2, max.overlaps = Inf)

emapplot(x_down)
cnetplot(x_down, categorySize="pvalue", foldChange = geneList_down)




######All plots for top 100 Down DEGS
reactome_down = as.data.frame(read_excel(dataset, sheet = "TopGenes_Down_100"))
head(reactome_down)
geneList_down<-reactome_down$entrez
head(geneList_down)
x_down<-enrichPathway(gene=geneList_down, pvalueCutoff = 0.05, readable = TRUE)
head(as.data.frame(x_down))
barplot(x_down, showCategory = 10)
dotplot(x_down, showCategory = 10)
BiocManager::install("enrichplot")
library(enrichplot)


x2 <- pairwise_termsim(x_down) 
emapplot(x2)
install.packages("ggnewscale")
library(ggnewscale)

emapplot(x_down)
cnetplot(x_down, categorySize="pvalue", foldChange = geneList_down)

########All the plots for UP
reactome_down = as.data.frame(read_excel(dataset, sheet = "TopGenes_Up_100"))
head(reactome_down)
geneList_down<-reactome_down$entrez
head(geneList_down)
x_down<-enrichPathway(gene=geneList_down, pvalueCutoff = 0.05, readable = TRUE)
head(as.data.frame(x_down))
barplot(x_down, showCategory = 10)
dotplot(x_down, showCategory = 10)
x2 <- pairwise_termsim(x_down) 
emapplot(x2)

emapplot(x_down)
cnetplot(x_down, categorySize="pvalue", foldChange = geneList_down)


# Drawing bubbleplot using ggplot2 package 
setwd("G:/Мой диск/LBSB_Aigul_Sharip/ESCC_DESEQ_New_Normal_ESCCsamples_April2020/DeSEQ_with11Normal_25Tumor sample/DAVID")
library(ggplot2) 
dataset = "DAVID_TumorVsNormal.xlsx"
par(mfrow=c(2,2))

#KR_UP
dataset = as.data.frame(read_excel(dataset, sheet = "TumorvsNormal_UP_KR_20"))
head(dataset)
ggplot(dataset, aes(x=Category, y=Term, size = Count, colour = Category)) + 
  geom_point()

#GO_UP
dataset = as.data.frame(read_excel(dataset, sheet = "TumorvsNormal_UP_GO_20"))
head(dataset)
ggplot(dataset, aes(x=Category, y=Term, size = Count, colour = Category)) + 
  geom_point()

#KR_DOWN
dataset = as.data.frame(read_excel(dataset, sheet = "TumorvsNormal_Down_KR_20"))
ggplot( dataset, aes(x=Category, y=Term, size = Count, colour = Category)) + 
  geom_point()

#GO_DOWN
dataset = as.data.frame(read_excel(dataset, sheet = "TumorvsNormal_Down_GO_20"))
ggplot(dataset, aes(x=Category, y=Term, size = Count, colour = Category)) + 
  geom_point()


############################################################
#Top 100 genes for Table S5 and Table S6
topGenes <- head(order(resTumorvsNormal$padj),100)
mat <- betas[topGenes, -c(1,1)]
pheatmap(mat)
View(mat)
write.csv(as.data.frame(mat), file="heatmap_matrix_res.csv")
#write.csv(as.data.frame(topGenes), file="topGenes_50_real.csv")
View(topGenes)
pheatmap(mat,display_numbers = TRUE)
mat050 = head(mat, 050); pheatmap(mat050)   
q050 = add_genes_and_codes(mat); pheatmap(q050) 
View(q050)

setwd("G:/My Drive/LBSB_Work/Kaz_ESCC/Bioinformatic analysis/KAZ_ESCC_DeSeq_22Tumor&11Normal/Outputs/")
write.csv(as.data.frame(q050), file="top100 genes for heatmap with genes names.csv")


############################################################
# Functions to filter the genes and adding gene names 
############################################################
filterRes <- function(res) {
  res <- res[which(res$padj<0.05),]
  res <- res[order(res$log2FoldChange),]
  res <-res[abs(res$log2FoldChange) >= 1.0,]
  return (res)
}
add_genes_names = function(res) {
  res$ensembl <- sapply(strsplit(rownames(res),split="\\+"),"[",1)
  ensembl=useMart("ensembl",dataset = "hsapiens_gene_ensembl")
  genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = res$ensembl, mart = ensembl )
  idx <- match( res$ensembl, genemap$ensembl_gene_id )
  res$entrez <- genemap$entrezgene[ idx ]
  res$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
  return (res)
}
#We can generate the list DEGs for each tumor stage using these functions
summary(resT1_N)
resT1_N = filterRes(resT1_N)
#Add gene names
resT1_N$ensembl <- sapply(strsplit(rownames(resT1_N),split="\\+"),"[",1)
ensembl=useMart("ensembl",dataset = "hsapiens_gene_ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = resT1_N$ensembl, mart = ensembl )
idx <- match( resT1_N$ensembl, genemap$ensembl_gene_id )
resT1_N$entrez <- genemap$entrezgene[ idx ]
resT1_N$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
summary(resT1_N)
head(resT1_N)
write.csv(as.data.frame(resT1_N), file="sign_genes_ordered_T1_vs_N_filteredwithgenename.CSV")

#For resT2_N
#We can generate the list DEGs for each tumor stage using these functions
summary(resT2_N)
resT2_N = filterRes(resT2_N)

#Add gene names
resT2_N$ensembl <- sapply(strsplit(rownames(resT2_N),split="\\+"),"[",1)
ensembl=useMart("ensembl",dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = resT2_N$ensembl, mart = ensembl)
idx <- match( resT2_N$ensembl, genemap$ensembl_gene_id )
resT2_N$entrez <- genemap$entrezgene[ idx ]
resT2_N$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
summary(resT2_N)
head(resT2_N)
write.csv(as.data.frame(resT2_N), file="sign_genes_ordered_T2_vs_N_filteredwithgenename.CSV")

#For resT3_N
#We can generate the list DEGs for each tumor stage using these functions
summary(resT3_N)
resT3_N = filterRes(resT3_N)
#Add gene names
resT3_N$ensembl <- sapply(strsplit(rownames(resT3_N),split="\\+"),"[",1)
ensembl=useMart("ensembl",dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = resT3_N$ensembl, mart = ensembl)
idx <- match(resT3_N$ensembl, genemap$ensembl_gene_id )
resT3_N$entrez <- genemap$entrezgene[ idx ]
resT3_N$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
summary(resT3_N)
head(resT3_N)
write.csv(as.data.frame(resT3_N), file="sign_genes_ordered_T3_vs_N_filteredwithgenename.CSV")

#For resT4_N
#We can generate the list DEGs for each tumor stage using these functions
summary(resT4_N)
resT4_N = filterRes(resT4_N)
#Add gene names
resT4_N$ensembl <- sapply(strsplit(rownames(resT4_N),split="\\+"),"[",1)
ensembl=useMart("ensembl",dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = resT4_N$ensembl, mart = ensembl)
idx <- match(resT4_N$ensembl, genemap$ensembl_gene_id )
resT4_N$entrez <- genemap$entrezgene[ idx ]
resT4_N$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
summary(resT4_N)
head(resT4_N)
write.csv(as.data.frame(resT4_N), file="sign_genes_ordered_T4_vs_N_filteredwithgenename.CSV")


###############Functional enrichment analysis with clusterProfiler################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MatrixGenerics")
BiocManager::install("SummarizedExperiment")
remove.packages(c("matrixStats", "DelayedArray", "MatrixGenerics", "SummarizedExperiment"))
BiocManager::install(c("matrixStats", "DelayedArray", "MatrixGenerics", "SummarizedExperiment"))


library(BiocManager)

BiocManager::install("clusterProfiler")
library (clusterProfiler)


# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install the clusterProfiler package
BiocManager::install("clusterProfiler") 

# Load the library
library(clusterProfiler)

# Load the DESeq2 package
library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts, colData, formula(~ condition)); 

# For this hypothetical example, let's assume you have a DESeq2
# result object named res that contains your differentially expressed genes
ddsEnrichment <- DESeq(dds) 


# res <- DESeq(dds) is the output that contains diff exp genes.
resEnrichment <-results(ddsEnrichment)

# Extract significant DEGs
sig_genes <- subset(resEnrichment, padj < 0.05)

# Extract gene symbols
gene_symbols <- row.names(sig_genes)

# Convert gene symbol to Entrez gene ID for KEGG analysis
summary(resTumorvsNormal)
select <- AnnotationDbi::select
entrez_ids <- select(org.Hs.eg.db, keys=rownames(resTumorvsNormal), columns="ENTREZID", keytype="ENSEMBL")
entrez_ids <- na.omit(entrez_ids)

# Perform KEGG pathway enrichment analysis
kk <- enrichKEGG(gene             = entrez_ids$ENTREZID,
                 organism         = 'hsa',
                 pvalueCutoff     = 0.05, 
                 pAdjustMethod    = "BH")

# Visualize the results
dotplot(kk)
barplot(kk)
cnetplot(kk)
enrichMap(kk)



######Draw heatmap for common consistently aberrantly expressed 1002 genes
### Data preparation
interest_genes_df <- read_excel("1002_common_genes.xlsx")
interest_genes <- interest_genes_df$EnsemblID
head(interest_genes)
colnames(resT1_N)
head(resT1_N)
# Create a subset of genes_data for the 1002 genes
extracted_data_T1 <- resT1_N[rownames(resT1_N)%in% interest_genes, ]
extracted_data_T2 <- resT2_N[rownames(resT2_N) %in% interest_genes, ]
extracted_data_T3 <- resT3_N[rownames(resT3_N) %in% interest_genes, ]
extracted_data_T4 <- resT4_N[rownames(resT4_N) %in% interest_genes, ]

# Print the extracted data
print(extracted_data)
write.csv(as.data.frame(extracted_data_T1), file="1002_common_genes_T1.csv")
write.csv(as.data.frame(extracted_data_T2), file="1002_common_genes_T2.csv")
write.csv(as.data.frame(extracted_data_T3), file="1002_common_genes_T3.csv")
write.csv(as.data.frame(extracted_data_T4), file="1002_common_genes_T4.csv")


##########Draw
install.packages("pheatmap")
install.packages("readxl")

library(pheatmap)
library(readxl)
data_df <- read_excel('Heatmap_common_1002 genes.xlsx', sheet = 'all')
rownames(data_df) <- data_df[[1]]  # Assuming the first column contains Ensembl IDs
data_df <- data_df[, -1]  # Remove the first column if it's used for row names
pheatmap(data_df, 
         cluster_rows = TRUE,  # Cluster rows
         cluster_cols = TRUE,  # Cluster columns
         main = "Heatmap Title",  # Title of the heatmap
         fontsize = 8,  # Font size for labels
         cellwidth = 20,  # Width of each cell
         cellheight = 10,  # Height of each cell
         scale = "row",  # Scale rows (other options include "column" or "none")
         border_color = NA  # No border color
)
pheatmap(data_df, cluster_cols=TRUE, show_rownames=TRUE)
