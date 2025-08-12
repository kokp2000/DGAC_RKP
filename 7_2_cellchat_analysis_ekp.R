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

setwd('D:/KP/EKP_RKP')



###################################################################

###################################################################

ekp <- readRDS('ekp_cellchat_by_celltype.rds')
levels(ekp@idents) # show factor levels of the cell labels
table(ekp@idents)

groupSize <- as.numeric(table(ekp@idents)) # number of cells in each cell group
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
ekp@DB <- CellChatDB.use

#Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
ekp <- subsetData(ekp) # This step is necessary even if using the whole database
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
ekp <- identifyOverExpressedGenes(ekp)
ekp <- identifyOverExpressedInteractions(ekp)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
#ekp <- projectData(ekp, PPI.mouse)
#Part II: Inference of cell-cell communication network
#Compute the communication probability and infer cellular communication network

unique(ekp@idents)
is.na(ekp@idents)
#ekp@meta$celltype = droplevels(ekp@meta$celltype, exclude = setdiff(levels(ekp@meta$celltype), unique(ekp@meta$celltype)))
ekp <- computeCommunProb(ekp)

#ekp <- computeCommunProb(
#  ekp,
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
ekp <- filterCommunication(ekp, min.cells = 10)

#Extract the inferred cellular communication network as a data frame
table(ekp@meta$celltype)


df.net <- subsetCommunication(ekp, sources.use = c(1,13,16,17,19), targets.use = c(2,3,4,5,6,7,8,9,10,11,12,14,15,18,20,21,22)) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

df.net <- subsetCommunication(ekp, signaling = c("CCL")) #gives the inferred cell-cell communications mediated by 
#Infer the cell-cell communication at a signaling pathway level
ekp <- computeCommunProbPathway(ekp)
#Calculate the aggregated cell-cell communication network
ekp <- aggregateNet(ekp)

groupSize <- as.numeric(table(ekp@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(ekp@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(ekp@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


mat <- ekp@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
ekp@netP$pathways
#    [1] "APP"         "CCL"         "FN1"         "MIF"         "SPP1"        "THBS"       
#[7] "CD80"        "LAMININ"     "CXCL"        "GALECTIN"    "TGFb"        "COMPLEMENT" 
#[13] "VISFATIN"    "CD45"        "CD23"        "CD22"        "CSF"         "ICAM"       
#[19] "TNF"         "CD39"        "CD52"        "PECAM1"      "SEMA4"       "ANNEXIN"    
#[25] "OSM"         "NOTCH"       "ITGAL-ITGB2" "JAM"         "CD48"        "NRG"        
#[31] "APRIL"       "BST2"        "CDH"         "LAIR1"  

pathways.show <- c("TGFb")

par(mfrow=c(1,1))
netVisual_aggregate(ekp, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(ekp, signaling = pathways.show, layout = "chord")




levels(ekp@idents) # show factor levels of the cell labels

#  [1] "B cell"          "Cd4 T cell"      "Cd8 T cell"      "M1_Macrophage"  
#  [5] "M2_Macrophage"   "Macrophage"      "epithelial cell"

# show all the interactions between epithelial cell and T cell 
#netVisual_chord_gene(ekp, sources.use = c(12), targets.use = c(2,8,9,11), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions between epithelial cell and Macrophage 
#netVisual_chord_gene(ekp, sources.use = c(12), targets.use = c(5,6), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions between epithelial cells and [T cell and all immune]
#netVisual_chord_gene(ekp, sources.use = c(12), targets.use = c(1,2,3,4,5,6,7,8,9,10,11,13), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)


# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(ekp, sources.use = c(2), targets.use = c(7),  signaling = c("TGFb"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_bubble(ekp, sources.use = c(3), targets.use = c(7), signaling = c("TGFb"), remove.isolate = FALSE)

netVisual_chord_gene(ekp, signaling = c("TGFb"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_chord_gene(ekp, signaling = c("NECTIN"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_chord_gene(ekp, signaling = c("CXCL"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_chord_gene(ekp, signaling = c("CCL"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_chord_gene(ekp, signaling = c("PD-L1"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)



#netVisual_chord_gene(ekp, sources.use = c(12), targets.use = c(1,2,3,4,5,6,7,8,9,11),  signaling = c("MIF"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
#netVisual_bubble(ekp, sources.use = c(12), targets.use = c(1,2,3,4,5,6,7,8,9,11), signaling = c("MIF"), remove.isolate = FALSE)
#netVisual_chord_cell(ekp, sources.use = c(12), targets.use = c(1,2,3,4,5,6,7,8,9,11),  signaling = c("MIF"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)


netVisual_chord_gene(ekp, sources.use = c(11), targets.use = c(2,3,5,6,7,9,10,11,12), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_bubble(ekp, sources.use = c(11), targets.use = c(2,3,5,6,7,9,10,11,13), remove.isolate = FALSE)
netVisual_bubble(ekp, sources.use = c(5,6,7,9), targets.use = c(2,3,10,11,12), remove.isolate = FALSE)
netVisual_bubble(ekp, sources.use = c(2,3,10,12), targets.use = c(5,6,7,9,11), remove.isolate = FALSE)


#> Comparing communications on a single object

saveRDS(ekp, file = 'ekp_after_cellchat.rds')
ekp <- readRDS('ekp_after_cellchat.rds')

