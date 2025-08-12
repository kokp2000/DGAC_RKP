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


# Loading patients data (prepared by Lana before)
GeneMatrix_DGAC = read.csv("D:/KP/EKP_RKP/Scissor/19_pts_raw_Pseudobulk_allcells.csv")
DGAC12_phenotype= read.csv('D:/KP/EKP_RKP/Scissor/19_pts_raw_Pseudobulk_allcells_phenotype.csv')

rownames(DGAC12_phenotype) = DGAC12_phenotype$Patient
Phenotype_DGAC12_input = DGAC12_phenotype[,!(colnames(DGAC12_phenotype) %in% c('Patient','Disease'))]
Phenotype_DGAC12_input
names(Phenotype_DGAC12_input) = rownames(DGAC12_phenotype)
Phenotype_DGAC12_input
table(Phenotype_DGAC12_input)
#Phenotype_DGAC12_input = read.csv('D:/KP/EKP_RKP/Scissor/19_pts_raw_Pseudobulk_allcells_phenotype_input.csv')



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


# 4. scRNA-seq data preprocess before
# Data prepared from '8_Data_prep_for_Phate_082524.ipynb'
library(data.table)
rna.df = as.data.frame(fread("D:/KP/EKP_RKP/Immune_EKP_RKP_raw_counts.csv", header = T))
rna.df[,1] = NULL


gene.df = read.table("D:/KP/EKP_RKP/Immune_EKP_RKP_genes.csv", header = T, as.is = T, sep = ",")
meta.df = read.table("D:/KP/EKP_RKP/Immune_EKP_RKP_metadata.csv", header = T, as.is = T, sep = ",")
rownames(meta.df) = meta.df$X
colnames(rna.df) = gene.df[,2]
rownames(rna.df) = rownames(meta.df)
t.rna.df = as.data.frame(t(rna.df))

####
#srna <- CreateSeuratObject(t.rna.df, meta = meta.df)
#srna@assays$RNA$data <- srna@assays$RNA$counts
#srna <- NormalizeData(srna)
###


# 5. run Scissor
# EKP vs RKP
phenotype <- Phenotype_DGAC12_input
tag <- c('DGAC1', 'DGAC2')

rownames(GeneMatrix_DGAC) = GeneMatrix_DGAC$genes
GeneMatrix_DGAC <- GeneMatrix_DGAC[,-1]
GeneMatrix_DGAC_input  <- as.matrix(GeneMatrix_DGAC)


infos <- Scissor(bulk_dataset=GeneMatrix_DGAC_input, sc_dataset=t.rna.df, phenotype=phenotype, tag = tag, alpha = 0.32, family = "binomial", Save_file = "D:/KP/EKP_RKP/Scissor/DGAC12.RData")
Scissor_select <- rep(0, ncol(t.rna.df))
names(Scissor_select) <- colnames(t.rna.df)
Scissor_select[infos$Scissor_pos] <- 1
Scissor_select[infos$Scissor_neg] <- 2
# sc_dataset <- AddMetaData(t.rna.df, metadata = Scissor_select, col.name = "scissor_0.32")
# Scissor_organoid <- DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor_0.32', cols = c('grey','red','blue'), order = c(2,1))
write.table(as.data.frame(Scissor_select), "D:/KP/EKP_RKP/Scissor/Scissor_DGAC12_immune_a0.32.csv", sep = ",", row.names = T, col.names = T, quote = F)    # export meta.data of sc_dataset, you can check Scissor information

infos <- Scissor(bulk_dataset=GeneMatrix_DGAC_input, sc_dataset=t.rna.df, phenotype=phenotype, tag = tag, alpha = 0.05, family = "binomial", Save_file = "D:/KP/EKP_RKP/Scissor/normal_DGAC_a0.05.RData")
Scissor_select <- rep(0, ncol(t.rna.df))
names(Scissor_select) <- colnames(t.rna.df)
Scissor_select[infos$Scissor_pos] <- 1
Scissor_select[infos$Scissor_neg] <- 2
write.table(as.data.frame(Scissor_select), "D:/KP/EKP_RKP/Scissor/Scissor_DGAC12_immune_a0.05.csv", sep = ",", row.names = T, col.names = T, quote = F)    # export meta.data of sc_dataset, you can check Scissor information

infos <- Scissor(bulk_dataset=GeneMatrix_DGAC_input, sc_dataset=t.rna.df, phenotype=phenotype, tag = tag, alpha = 0.01, family = "binomial", Save_file = "D:/KP/EKP_RKP/Scissor/normal_DGAC_a0.01.RData")
Scissor_select <- rep(0, ncol(t.rna.df))
names(Scissor_select) <- colnames(t.rna.df)
Scissor_select[infos$Scissor_pos] <- 1
Scissor_select[infos$Scissor_neg] <- 2
write.table(as.data.frame(Scissor_select), "D:/KP/EKP_RKP/Scissor/Scissor_DGAC12_immune_a0.01.csv", sep = ",", row.names = T, col.names = T, quote = F)    # export meta.data of sc_dataset, you can check Scissor information


# Make UMAP in scanpy again
