rm(list = ls())
install.packages("ade4")

install.packages('NMF')
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")
BiocManager::install("Biobase")
BiocManager::install("BiocNeighbors")
devtools::install_github("sqjin/CellChat")
##### CellChat download blocked due to the MD Anderson firewall
## download in the mdaguest wifi


library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(BiocNeighbors)
library(ade4)
library(CellChat)

setwd('D:/KP/rkp_RKP')



###################################################################

###################################################################

rkp <- readRDS('rkp_cellchat_by_celltype.rds')
levels(rkp@idents) # show factor levels of the cell labels
table(rkp@idents)

groupSize <- as.numeric(table(rkp@idents)) # number of cells in each cell group
#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
rkp@DB <- CellChatDB.use

#Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
rkp <- subsetData(rkp) # This step is necessary even if using the whole database
future::plan("multisession", workers = 8) # do parallel
options(future.globals.maxSize= 943718400)

#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
rkp <- identifyOverExpressedGenes(rkp)
rkp <- identifyOverExpressedInteractions(rkp)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
#rkp <- projectData(rkp, PPI.mouse)
#Part II: Inference of cell-cell communication network
#Compute the communication probability and infer cellular communication network

unique(rkp@idents)
is.na(rkp@idents)
#rkp@meta$celltype = droplevels(rkp@meta$celltype, exclude = setdiff(levels(rkp@meta$celltype), unique(rkp@meta$celltype)))
rkp <- computeCommunProb(rkp)

#rkp <- computeCommunProb(
#  rkp,
#  type = c("triMean", "truncatedMean", "median"),
#  trim = NULL,
#  LR.use = NULL,
#  raw.use = TRUE,
#  population.size = FALSE,
#  do.fast = TRUE,
#  nboot = 100,
#  seed.use = 1L,
#  Kh = 0.5,
#  n = 1
#)


# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
rkp <- filterCommunication(rkp, min.cells = 10)

#Extract the inferred cellular communication network as a data frame
table(rkp@meta$celltype)


df.net <- subsetCommunication(rkp, sources.use = c(1,13,16,17,19), targets.use = c(2,3,4,5,6,7,8,9,10,11,12,14,15,18,20,21,22)) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

df.net <- subsetCommunication(rkp, signaling = c("CCL")) #gives the inferred cell-cell communications mediated by 
#Infer the cell-cell communication at a signaling pathway level
rkp <- computeCommunProbPathway(rkp)
#Calculate the aggregated cell-cell communication network
rkp <- aggregateNet(rkp)

groupSize <- as.numeric(table(rkp@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(rkp@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(rkp@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


mat <- rkp@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
rkp@netP$pathways
#    [1] "CCL"         "MHC-I"       "CD80"        "MIF"         "LAMININ"     "SPP1"       
#[7] "GALECTIN"    "CXCL"        "ALCAM"       "CD6"         "CD22"        "CD45"       
#[13] "PD-L1"       "FN1"         "CD52"        "THY1"        "ICAM"        "LCK"        
#[19] "ITGAL-ITGB2" "SELPLG"      "TNF"         "CD86"        "MHC-II"      "CD23"       
#[25] "PECAM1"      "COMPLEMENT"  "TGFb"        "PDL2"        "SEMA4"       "VISFATIN"   
#[31] "ICOS"        "FASLG"       "IFN-II"      "ANNEXIN"     "IL16"  

pathways.show <- c("PD-L1")

par(mfrow=c(1,1))
netVisual_aggregate(rkp, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(rkp, signaling = pathways.show, layout = "chord")




levels(rkp@idents) # show factor levels of the cell labels

#  [1] "B cell"          "Cd4 T cell"      "Cd8 T cell"      "M1_Macrophage"  
#  [5] "M2_Macrophage"   "Macrophage"      "T cell"          "epithelial cell"

# show all the interactions between epithelial cell and T cell 
#netVisual_chord_gene(rkp, sources.use = c(12), targets.use = c(2,8,9,11), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions between epithelial cell and Macrophage 
#netVisual_chord_gene(rkp, sources.use = c(12), targets.use = c(5,6), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions between epithelial cells and [T cell and all immune]
#netVisual_chord_gene(rkp, sources.use = c(12), targets.use = c(1,2,3,4,5,6,7,8,9,10,11,13), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)


# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(rkp, sources.use = c(2), targets.use = c(8),  signaling = c("PD-L1"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_bubble(rkp, sources.use = c(3), targets.use = c(8), signaling = c("PD-L1"), remove.isolate = FALSE)

netVisual_chord_gene(rkp, signaling = c("PD-L1"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_chord_gene(rkp, signaling = c("CD80"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_chord_gene(rkp, signaling = c("MHC-II"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_chord_gene(rkp, signaling = c("MHC-I"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_chord_gene(rkp, signaling = c("PDL2"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)

netVisual_aggregate(rkp, signaling = c("CCL"), layout = "circle")
netVisual_aggregate(rkp, signaling = c("CXCL"), layout = "circle")



#netVisual_chord_gene(rkp, sources.use = c(12), targets.use = c(1,2,3,4,5,6,7,8,9,11),  signaling = c("MIF"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
#netVisual_bubble(rkp, sources.use = c(12), targets.use = c(1,2,3,4,5,6,7,8,9,11), signaling = c("MIF"), remove.isolate = FALSE)
#netVisual_chord_cell(rkp, sources.use = c(12), targets.use = c(1,2,3,4,5,6,7,8,9,11),  signaling = c("MIF"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)


netVisual_chord_gene(rkp, sources.use = c(11), targets.use = c(2,3,5,6,7,9,10,11,12), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_bubble(rkp, sources.use = c(11), targets.use = c(2,3,5,6,7,9,10,11,13), remove.isolate = FALSE)
netVisual_bubble(rkp, sources.use = c(5,6,7,9), targets.use = c(2,3,10,11,12), remove.isolate = FALSE)
netVisual_bubble(rkp, sources.use = c(2,3,10,12), targets.use = c(5,6,7,9,11), remove.isolate = FALSE)


#> Comparing communications on a single object

saveRDS(rkp, file = 'rkp_after_cellchat.rds')
rkp <- readRDS('rkp_after_cellchat.rds')

