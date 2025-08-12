install.packages(
  "https://cran.r-project.org/src/contrib/Archive/Rmagic/Rmagic_2.0.3.tar.gz",
  repos = NULL,
  type = "source"
)

if (!require(viridis)) install.packages("viridis")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(readr)) install.packages("readr")
if (!require(Rmagic)) install.packages("Rmagic")

system("pip install wrapt")
system("pip install future")
system("pip install matplotlib")
system("pip install pandas")



rm(list = ls())

library(data.table)
library(phateR)
library(ggplot2)
library(readr)
library(viridis)
library(Rmagic)
library(reticulate)

reticulate::py_config()

#use_condaenv('KP_py38') 
use_condaenv('phate_KP')

reticulate::py_config()
reticulate::py_install("phate", pip=TRUE)

# Prepare input data
count.df = as.data.frame(fread("D:/KP/EKP_RKP/EKP_RKP_raw_counts.csv"))

count.df_1 = count.df
rownames(count.df_1) = count.df_1[,1]
count.df_1 = count.df_1[, -1]
count.df = count.df_1

gene.df = read.table("D:/KP/EKP_RKP/EKP_RKP_genes.csv", header = T, as.is = T, sep = ",")
meta.df = read.table("D:/KP/EKP_RKP/EKP_RKP_metadata.csv", header = T, as.is = T, sep = ",")
rownames(meta.df) = meta.df[,1]
rownames(gene.df) = gene.df[,2]

output.dir = "D:/KP/EKP_RKP/phate"

#input.df = as.data.frame(t(gene.df)) 
# Row as Cell, Column as Gene
input.df = as.data.frame(count.df)

# Keep genes expressed in at least 10 cells
keep_cols <- colSums(input.df > 0) > 10
input.df <- input.df[,keep_cols]

# Keep cells with at least 1000 UMIs
keep_rows <- rowSums(input.df) > 1000
input.df <- input.df[keep_rows,]

input.df <- library.size.normalize(input.df)
input.df <- sqrt(input.df)

input_PCA <- as.data.frame(prcomp(input.df)$x)

# Run PHATE
input_PHATE <- phate(input_PCA)
output.df = as.data.frame(input_PHATE[[1]])
write.table(output.df, paste0(output.dir, "/PHATE.output.txt"), sep = "\t", row.names = T, col.names = T, quote = F)
