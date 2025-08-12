rm(list = ls())

library(Seurat)
library(SeuratData)
library(ggplot2)
library(sctransform)
library(dplyr)
library(patchwork)
library(SeuratDisk)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DT)
library(biomaRt)
library(Scissor)

# memory.limit(512000)
# options(future.globals.maxSize = 256000*1024^2)

output.dir = "D:/KP/EKP_RKP/Scissor"
dir.create(output.dir, showWarnings = F)

# 1. prepare mouse gene symbol annotated TCGA-STAD RNA-seq gene matrix
# download TCGA-STAD gene matrix including 'Primary Tumor' and 'Solid Tissue Normal'
query <- GDCquery(project = "TCGA-STAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  access= 'open', data.format='TXT', experimental.strategy='RNA-Seq')
GDCdownload(query)
RNA_seq <- GDCprepare(query=query, save=TRUE, save.filename='TCGA_STAD_Exp.rda')
GeneMatrix <- assay(RNA_seq)
write.csv(GeneMatrix, file=paste0(output.dir, "/TCGA_STAD_GeneMatrix.csv"), sep = ",", row.names=TRUE, col.names = TRUE)
# GeneMatrix = read.csv(paste0(output.dir, "/TCGA_STAD_GeneMatrix.csv"), sep = ",", header = T, row.names = 1)
# colnames(GeneMatrix) = gsub("\\.", "-", colnames(GeneMatrix))

# convert human ensembl id to mouse gene symbol
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host='https://dec2021.archive.ensembl.org/')
mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", host='https://dec2021.archive.ensembl.org/')
# Use follows if makes error
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
####

human_gene <- rownames(GeneMatrix)
mouse_gene <- getLDS(attributes = c("ensembl_gene_id_version"), filters = "ensembl_gene_id_version", values = human_gene, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=TRUE)
head(mouse_gene)

mouse_gene <- mouse_gene[!(mouse_gene$MGI.symbol==""),]
mouse_gene <- mouse_gene[!duplicated(mouse_gene$Gene.stable.ID.version),]
mouse_gene <- mouse_gene[!duplicated(mouse_gene$MGI.symbol),]

# make the mouse gene symbol annotated gene matrix
class(GeneMatrix)
GeneMatrix_df<- as.data.frame(GeneMatrix)
class(GeneMatrix_df)

GeneMatrix_MGI <- GeneMatrix_df
GeneMatrix_MGI$genes <- mouse_gene$MGI.symbol[match(rownames(GeneMatrix_df), mouse_gene$Gene.stable.ID.version)]
GeneMatrix_MGI <- na.omit(GeneMatrix_MGI)
rownames(GeneMatrix_MGI) <- GeneMatrix_MGI$genes
GeneMatrix_MGI = GeneMatrix_MGI[,!(colnames(GeneMatrix_MGI) %in% 'genes')]
write.csv(GeneMatrix_MGI, file=paste0(output.dir, "/GeneMatrix_MGI.csv"), sep = ",", row.names=TRUE, col.names = TRUE)
#GeneMatrix_MGI = read.csv('D:/KP/EKP_RKP/Scissor/GeneMatrix_MGI.csv', row.names=1)



# 2. prepare phenotype list
# download TCGA-STAD clinical index data
clinical_index <- GDCquery_clinic(project = "TCGA-STAD", type = "Clinical", save.csv = TRUE)

# 2. get clinical info

query <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Clinical",
  data.format = "bcr xml"
)
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")
datatable(clinical, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)
write.csv(clinical, file = paste0(output.dir, "/clinical_STAD.csv"), sep = ",", row.names=TRUE, col.names = TRUE)
#clinical = read.csv("D:/KP/EKP_RKP/Scissor/clinical_STAD.csv")
clinical$bcr_patient_barcode <- gsub("-", ".", clinical$bcr_patient_barcode)


#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################

# Subsetting Signet ring cell carcinoma

SRCC_clinical = subset(clinical, histological_type %in% "Stomach Adenocarcinoma, Signet Ring Type")
GeneMatrix_SRCC <- subset(GeneMatrix_MGI, select = substr(colnames(GeneMatrix_MGI),1,12) %in% SRCC_clinical$bcr_patient_barcode)
dim(GeneMatrix_SRCC)

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################

# Normal
GeneMatrix_Normal <- subset(GeneMatrix_MGI, select = substr(colnames(GeneMatrix_MGI),14,15) %in% '11')
write.csv(GeneMatrix_Normal, file="D:/KP/EKP_RKP/Scissor/GeneMatrix_Normal.csv", sep = ",", row.names=TRUE, col.names = TRUE)
Normal_phenotype <- as.data.frame(colnames(GeneMatrix_Normal))
Normal_phenotype$Disease <- 'Normal'
Normal_phenotype$Status <- 0
colnames(Normal_phenotype)<-  c("Patient", 'Disease', "Status")

# SRCC
SRCC_phenotype <- as.data.frame(colnames(GeneMatrix_SRCC))
SRCC_phenotype$Disease <- 'SRCC'
SRCC_phenotype$Status <- 1
colnames(SRCC_phenotype)<-  c("Patient", 'Disease', "Status")
Phenotype_normal_SRCC <- rbind(Normal_phenotype, SRCC_phenotype)
write.table(Phenotype_normal_SRCC, file="D:/KP/EKP_RKP/Scissor/Phenotype_normal_SRCC.csv", sep = ",", row.names=FALSE)

# 3. match gene matrix and phenotype list
# only include patients that also exist in the phenotype list
GeneMatrix_normal_SRCC <- subset(GeneMatrix_MGI, select = colnames(GeneMatrix_MGI) %in% Phenotype_normal_SRCC$Patient)
patient = as.data.frame(colnames(GeneMatrix_normal_SRCC))

duplicated(Phenotype_normal_SRCC$Patient)
Phenotype_normal_SRCC_1 <- Phenotype_normal_SRCC[!duplicated(Phenotype_normal_SRCC$Patient),]
Phenotype_normal_SRCC <- Phenotype_normal_SRCC_1

# make gene matrix and phenotype list the same order
GeneMatrix_normal_SRCC <- GeneMatrix_normal_SRCC[ , order(match(colnames(GeneMatrix_normal_SRCC), Phenotype_normal_SRCC$Patient))]
all(colnames(GeneMatrix_normal_SRCC) == Phenotype_normal_SRCC$Patient)

# make the sample id annotated phenotype list
rownames(Phenotype_normal_SRCC) <- Phenotype_normal_SRCC$Patient
Phenotype_normal_SRCC_input = Phenotype_normal_SRCC[,!(colnames(Phenotype_normal_SRCC) %in% c('Patient','Disease'))]
Phenotype_normal_SRCC_input
names(Phenotype_normal_SRCC_input) = rownames(Phenotype_normal_SRCC)
Phenotype_normal_SRCC_input
table(Phenotype_normal_SRCC_input)

# 4. scRNA-seq data preprocess before
# Data prepared from '8_Data_prep_for_Phate_082524.ipynb'
library(data.table)
rna.df = as.data.frame(fread("D:/KP/EKP_RKP/EKP_RKP_raw_counts.csv", header = T))
rna.df[,1] = NULL


gene.df = read.table("D:/KP/EKP_RKP/EKP_RKP_genes.csv", header = T, as.is = T, sep = ",")
meta.df = read.table("D:/KP/EKP_RKP/EKP_RKP_metadata.csv", header = T, as.is = T, sep = ",")
rownames(meta.df) = meta.df$X
colnames(rna.df) = gene.df[,2]
rownames(rna.df) = rownames(meta.df)
t.rna.df = as.data.frame(t(rna.df))




# 5. run Scissor
# Seurat v5 does not work. Downgrade to version 4.4.0
# EKP vs RKP
phenotype <- Phenotype_normal_SRCC_input
tag <- c('Normal', 'SRCC')
GeneMatrix_normal_SRCC_input  <- as.matrix(GeneMatrix_normal_SRCC)

infos <- Scissor(bulk_dataset=GeneMatrix_normal_SRCC_input, sc_dataset=t.rna.df, phenotype=phenotype, tag = tag, alpha = 0.32, family = "binomial", Save_file = "D:/KP/EKP_RKP/Scissor/normal_SRCC.RData")
Scissor_select <- rep(0, ncol(t.rna.df))
names(Scissor_select) <- colnames(t.rna.df)
Scissor_select[infos$Scissor_pos] <- 1
Scissor_select[infos$Scissor_neg] <- 2
# sc_dataset <- AddMetaData(t.rna.df, metadata = Scissor_select, col.name = "scissor_0.32")
# Scissor_organoid <- DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor_0.32', cols = c('grey','red','blue'), order = c(2,1))
write.table(as.data.frame(Scissor_select), "D:/KP/EKP_RKP/Scissor/Scissor_normal_SRCC_a0.32.csv", sep = ",", row.names = T, col.names = T, quote = F)    # export meta.data of sc_dataset, you can check Scissor information

infos <- Scissor(bulk_dataset=GeneMatrix_normal_SRCC_input, sc_dataset=t.rna.df, phenotype=phenotype, tag = tag, alpha = 0.05, family = "binomial", Save_file = "D:/KP/EKP_RKP/Scissor/normal_SRCC_a0.05.RData")
Scissor_select <- rep(0, ncol(t.rna.df))
names(Scissor_select) <- colnames(t.rna.df)
Scissor_select[infos$Scissor_pos] <- 1
Scissor_select[infos$Scissor_neg] <- 2
write.table(as.data.frame(Scissor_select), "D:/KP/EKP_RKP/Scissor/Scissor_normal_SRCC_a0.05.csv", sep = ",", row.names = T, col.names = T, quote = F)    # export meta.data of sc_dataset, you can check Scissor information

infos <- Scissor(bulk_dataset=GeneMatrix_normal_SRCC_input, sc_dataset=t.rna.df, phenotype=phenotype, tag = tag, alpha = 0.01, family = "binomial", Save_file = "D:/KP/EKP_RKP/Scissor/normal_SRCC_a0.01.RData")
Scissor_select <- rep(0, ncol(t.rna.df))
names(Scissor_select) <- colnames(t.rna.df)
Scissor_select[infos$Scissor_pos] <- 1
Scissor_select[infos$Scissor_neg] <- 2
write.table(as.data.frame(Scissor_select), "D:/KP/EKP_RKP/Scissor/Scissor_normal_SRCC_a0.01.csv", sep = ",", row.names = T, col.names = T, quote = F)    # export meta.data of sc_dataset, you can check Scissor information


# Make UMAP in scanpy again
