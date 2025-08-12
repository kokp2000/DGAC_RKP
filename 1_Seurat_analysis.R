library(Seurat)
library(future)
plan("multisession", workers = 10)
library(ggplot2)
library(sf)

options(future.globals.maxSize = 8000 * 1024^2)


setwd("D:/KP/EKP_RKP/Xenium/data_from_R/")

# Load the RKP Xenium data
path <- "S:/Dept. Exp Rad Oncology/Lab Files/Lab Park/JPark/Xenium_data/Xenium_06-21-24/output-XETG00320__0038464__Region_3__20240619__215126"
xenium_RKP.obj <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
xenium_RKP.obj <- subset(xenium_RKP.obj, subset = nCount_Xenium > 0)
VlnPlot(xenium_RKP.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(xenium_RKP.obj, fov = "fov", molecules = c("Muc1", "Muc4", "Aqp5", "Tff2"), nmols = 20000)


# Load the EKP Xenium data
path <- "S:/Dept. Exp Rad Oncology/Lab Files/Lab Park/JPark/Xenium_data/Xenium_06-21-24/output-XETG00320__0038464__Region_4__20240619__215126"
xenium_EKP.obj <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
xenium_EKP.obj <- subset(xenium_EKP.obj, subset = nCount_Xenium > 0)
VlnPlot(xenium_EKP.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(xenium_EKP.obj, fov = "fov", molecules = c("Muc1", "Muc4", "Aqp5", "Tff2"), nmols = 20000)


# Expression level
ImageFeaturePlot(xenium_RKP.obj, fov="fov", features = c("Muc1", "Muc4", "Aqp5", "Tff2"), max.cutoff = c(25,35, 12, 5), size = 1.5, cols = c("white", "red"))                                                                                          
ImageFeaturePlot(xenium_EKP.obj, fov="fov", features = c("Muc1", "Muc4", "Aqp5", "Tff2"), max.cutoff = c(25,35, 12, 5), size = 1.5, cols = c("white", "red"))




# Zoomed image by crop-RKP
cropped.coords <- Crop(xenium_RKP.obj[["fov"]], x = c(1200, 2900), y = c(3750, 4550), coords = "plot")
xenium_RKP.obj[["zoom"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium_RKP.obj[["zoom"]]) <- "segmentation"
ImageDimPlot(xenium_RKP.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = FALSE, molecules = c("Muc1", "Muc4"), nmols = 10000)+
  theme(panel.grid = element_blank())

# Zoomed image by crop-EKP
cropped.coords_EKP <- Crop(xenium_EKP.obj[["fov"]], x = c(1200, 2900), y = c(3750, 4550), coords = "plot")
xenium_EKP.obj[["zoom"]] <- cropped.coords_EKP
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium_EKP.obj[["zoom"]]) <- "segmentation"
ImageDimPlot(xenium_EKP.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = FALSE, molecules = c("Muc1", "Muc4"), nmols = 10000)+
  theme(panel.grid = element_blank())



# Normalization and clustering - RKP
xenium_RKP.obj <- SCTransform(xenium_RKP.obj, assay = "Xenium")
xenium_RKP.obj <- RunPCA(xenium_RKP.obj, npcs = 30, features = rownames(xenium_RKP.obj))
xenium_RKP.obj <- RunUMAP(xenium_RKP.obj, dims = 1:30)
xenium_RKP.obj <- FindNeighbors(xenium_RKP.obj, reduction = "pca", dims = 1:30)
xenium_RKP.obj <- FindClusters(xenium_RKP.obj, resolution = 0.3)
DimPlot(xenium_RKP.obj)


# Normalization and clustering - EKP
xenium_EKP.obj <- SCTransform(xenium_EKP.obj, assay = "Xenium")
xenium_EKP.obj <- RunPCA(xenium_EKP.obj, npcs = 30, features = rownames(xenium_EKP.obj))
xenium_EKP.obj <- RunUMAP(xenium_EKP.obj, dims = 1:30)
xenium_EKP.obj <- FindNeighbors(xenium_EKP.obj, reduction = "pca", dims = 1:30)
xenium_EKP.obj <- FindClusters(xenium_EKP.obj, resolution = 0.3)
DimPlot(xenium_EKP.obj)



ImageDimPlot(xenium_RKP.obj, fov = "fov",border.size = 0, axes = TRUE, , cols = "polychrome",
             coord.fixed = FALSE, nmols = 10000, group.by = 'seurat_clusters')+theme(panel.grid = element_blank())
ImageDimPlot(xenium_EKP.obj, fov = "fov",border.size = 0, axes = TRUE, , cols = "polychrome",
             coord.fixed = FALSE, nmols = 10000, group.by = 'seurat_clusters')+theme(panel.grid = element_blank())


FeaturePlot(xenium_RKP.obj, features = c("Muc1", "Muc4"),max.cutoff = 2)
FeaturePlot(xenium_EKP.obj, features = c("Muc1", "Muc4"),max.cutoff = 2)


# Importing the cluster information used in Xenium Explorer
cluster_EKP = read.csv("S:/Dept. Exp Rad Oncology/Lab Files/Lab Park/JPark/Xenium_data/Xenium_06-21-24/output-XETG00320__0038464__Region_4__20240619__215126/analysis/clustering/gene_expression_kmeans_10_clusters/clusters.csv")
head(cluster_EKP)
xenium_EKP.obj@meta.data[cluster_EKP$Barcode, "kmean_cluster10_EKP"] = cluster_EKP$Cluster
head(xenium_EKP.obj@meta.data)
DimPlot(xenium_EKP.obj, group.by = "kmean_cluster10_EKP")
# rename clusters by as it was annotated in Xenium explorer
Idents(xenium_EKP.obj) = "kmean_cluster10_EKP"
table(Idents(xenium_EKP.obj))
xenium_EKP.obj <- RenameIdents(xenium_EKP.obj, `1` = "effector_T_cell", `2` = "Tumor_1", `3` = "Myeloid_cell_1", 
                                    `4` = "Tumor_3/Lgr5+", `5` = "exhausted_T_cell", `6` = "Myeloid_cell_3", 
                               `7` = "Tumor_2", `8` = "Fibroblast_like", 
                                    `9` = "Endothelial_cell", `10` = "Myeloide_cell_2" )
xenium_EKP.obj$celltype <- Idents(xenium_EKP.obj)

saveRDS(xenium_EKP.obj, "EKP_clustered.rds")
#xenium_EKP.obj = readRDS("EKP_clustered.rds")



# Importing the cluster information used in Xenium Explorer
cluster_RKP = read.csv("S:/Dept. Exp Rad Oncology/Lab Files/Lab Park/JPark/Xenium_data/Xenium_06-21-24/output-XETG00320__0038464__Region_3__20240619__215126/analysis/clustering/gene_expression_kmeans_10_clusters/clusters.csv")
head(cluster_RKP)
xenium_RKP.obj@meta.data[cluster_RKP$Barcode, "kmean_cluster10_RKP"] = cluster_RKP$Cluster
head(xenium_RKP.obj@meta.data)
str(xenium_RKP.obj@meta.data)
DimPlot(xenium_RKP.obj, group.by = "kmean_cluster10_RKP")
# rename clusters by as it was annotated in Xenium explorer
length(Idents(xenium_RKP.obj))
Idents(xenium_RKP.obj) = "kmean_cluster10_RKP"
table(Idents(xenium_RKP.obj))
xenium_RKP.obj <- RenameIdents(xenium_RKP.obj, `1` = "Tumor_1", `2` = "Myeloid_cell_1", 
                               `3` = "Fibroblast-like",`4` = "Tumor_2", `5` = "Myeloid_cell_2", 
                               `6` = "Endothelial_cell", 
                               `7` = "Tumor_3/Tff2+", `8` = "exhausted_T_cell", 
                               `9` = "Others_1", `10` = "Others_2" )
xenium_RKP.obj$celltype <- Idents(xenium_RKP.obj)

saveRDS(xenium_RKP.obj, "RKP_clustered.rds")
#xenium_RKP.obj = readRDS("RKP_clustered.rds")



#### Tumor cell subset #####
EKP_tumor <- c("Tumor_1","Tumor_2","Tumor_3/Lgr5+")
RKP_tumor <- c("Tumor_1","Tumor_2","Tumor_3/Tff2+")

xenium_EKP_tumor <- subset(xenium_EKP.obj, idents == EKP_tumor)
xenium_RKP_tumor <- subset(xenium_RKP.obj, idents == RKP_tumor)

DimPlot(xenium_EKP.obj)
DimPlot(xenium_EKP_tumor)

DimPlot(xenium_RKP.obj)
DimPlot(xenium_RKP_tumor)

DotPlot(xenium_EKP_tumor, features = c("Axin2","Tbx3",'Lgr5','Tff2'), dot.scale = 8, cols = "RdYlBu") +RotatedAxis()+ labs(y="group", x = "Genes") + theme(axis.text=element_text(size=13)) +theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill=NA))
DotPlot(xenium_RKP_tumor, features = c("Axin2","Tbx3",'Lgr5','Tff2'), dot.scale = 8, cols = "RdYlBu") +RotatedAxis()+ labs(y="group", x = "Genes") + theme(axis.text=element_text(size=13)) +theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill=NA))



xenium_EKP_tumor <- RenameIdents(xenium_EKP_tumor, `Tumor_1` = "EKP_Tumor_1", `Tumor_2` = "EKP_Tumor_2", 
                               `Tumor_3/Lgr5+` = "EKP_Tumor_3")
xenium_EKP_tumor$celltype_EKP <- Idents(xenium_EKP_tumor)
xenium_RKP_tumor <- RenameIdents(xenium_RKP_tumor, `Tumor_1` = "RKP_Tumor_1", `Tumor_2` = "RKP_Tumor_2", 
                                 `Tumor_3/Tff2+` = "RKP_Tumor_3")
xenium_RKP_tumor$celltype_RKP <- Idents(xenium_RKP_tumor)

saveRDS(xenium_EKP_tumor, "EKP_tumor.rds")
saveRDS(xenium_RKP_tumor, "RKP_tumor.rds")
#xenium_EKP_tumor = readRDS("EKP_tumor.rds")
#xenium_RKP_tumor = readRDS("RKP_tumor.rds")


#### Tumor cells of EKP and RKP integration ####
###### Merge done in other workstation by Jinho ######

options(future.globals.maxSize = 2000000000000)

Tumor <- merge(xenium_EKP_tumor, y=xenium_RKP_tumor, add.cell.ids=c("EKP","RKP"))
TUmor$orig.ident
Tumor$celltype
table(Tumor$celltype)
head(colnames(Tumor))
tail(colnames(Tumor))

RKP_EKP <- merge(xenium_EKP.obj, y=xenium_RKP.obj, add.cell.ids=c("EKP","RKP"))
#########################################################################################

EKP_RKP = readRDS("EPK_RPK.rds")
table(EKP_RKP$celltype)
table(Idents(EKP_RKP))
head(colnames(EKP_RKP))      
tail(colnames(EKP_RKP))
unique(sapply(X = strsplit(colnames(EKP_RKP), split = "_"), FUN = "[", 1))

head(EKP_RKP@meta.data)
split_names = strsplit(colnames(EKP_RKP), split = "_")
genotype <- sapply(split_names, '[',1)
head(genotype)
EKP_RKP <- AddMetaData(EKP_RKP, metadata = genotype, col.name = "genotype")
head(EKP_RKP@meta.data)
tail(EKP_RKP@meta.data)

combined_metadata <- paste(EKP_RKP@meta.data$genotype, EKP_RKP@meta.data$celltype, sep = "_")
EKP_RKP <- AddMetaData(EKP_RKP, metadata = combined_metadata, col.name = "genotype_celltype")
head(EKP_RKP@meta.data)


###################### sctransform is left but no xenium in assay ###########################

EKP_RKP <- RunPCA(EKP_RKP, npcs = 30, features = rownames(EKP_RKP))
EKP_RKP <- RunUMAP(EKP_RKP, dims = 1:30)
EKP_RKP <- FindNeighbors(EKP_RKP, reduction = "pca", dims = 1:30)
EKP_RKP <- FindClusters(EKP_RKP, resolution = 0.3)


DimPlot(EKP_RKP, reduction = "umap",group.by = "genotype_celltype", label = F,label.size = 5)
DimPlot(EKP_RKP, reduction = "umap",group.by = "celltype",label = F,label.size = 5)
DimPlot(EKP_RKP, reduction = "umap", label = F,label.size = 5)

Idents(EKP_RKP) = "genotype_celltype"
EKP_RKP <- RenameIdents(EKP_RKP, 
                               `EKP_Myeloide_cell_2` = "Myeloid", `EKP_Myeloid_cell_1` = "Myeloid", 
                               `EKP_Fibroblast_like` = "Fibroblast",`EKP_Myeloid_cell_3` = "Myeloid",
                               `EKP_exhausted_T_cell` = "T cell", `EKP_Endothelial_cell` = "Endothelial", 
                               `EKP_Tumor_3/Lgr5+` = "EKP_Tumor_3", `EKP_Tumor_1` = "EKP_Tumor_1",
                               `EKP_Tumor_2` = "EKP_Tumor_2", `EKP_effector_T_cell` = "T cell",
                               
                               `RKP_Myeloid_cell_2` = "Myeloid", `RKP_Endothelial_cell` = "Endothelial", 
                               `RKP_Myeloid_cell_1` = "Myeloid",`RKP_Fibroblast-like` = "Fibroblast",
                               `RKP_Others_2` = "Others_2", `RKP_exhausted_T_cell` = "T cell", 
                               `RKP_Tumor_1` = "RKP_Tumor_1", `RKP_Tumor_2` = "RKP_Tumor_2",
                               `RKP_Others_1` = "Others_1", `RKP_Tumor_3/Tff2+` = "RKP_Tumor_3"
                               )

table(Idents(EKP_RKP))
EKP_RKP$genotype_celltype = Idents(EKP_RKP)

my_levels <- c("EKP_Tumor_1","EKP_Tumor_2","EKP_Tumor_3",
               "RKP_Tumor_1","RKP_Tumor_2","RKP_Tumor_3",
               "T cell","Myeloid","Fibroblast","Endothelial","Others_1","Others_2")
EKP_RKP$genotype_celltype <- factor(x = EKP_RKP$genotype_celltype, levels = my_levels)

color_list = c("#23BD3D","#5DDEB1","#096B08",
               "#2FBDEC","#346CF0","#4136EE",
               "#FF0408","#F2E827","#9F6668","#FF7400","#A9A7A9","#838685")
DimPlot(EKP_RKP, reduction = "umap",group.by = "genotype_celltype",label = F,label.size = 5, cols = color_list)

saveRDS(EKP_RKP, "EKP_RKP_merged.rds")
#EKP_RKP = readRDS("EKP_RKP_merged.rds")

################### Subsetting Tumor and Immune cells because they are not separated well ####################
Tumor <- subset(EKP_RKP, idents = c("RKP_Tumor_1","RKP_Tumor_2","RKP_Tumor_3","EKP_Tumor_1","EKP_Tumor_2","EKP_Tumor_3"))
Immn <- subset(EKP_RKP, idents = c("Myeloid","T cell"))

Tumor <- RunPCA(Tumor, npcs = 30, features = rownames(Tumor))
Tumor <- RunUMAP(Tumor, dims = 1:30)
Tumor <- FindNeighbors(Tumor, reduction = "pca", dims = 1:30)
Tumor <- FindClusters(Tumor, resolution = 0.3)
DimPlot(Tumor, reduction = "umap", label = F,label.size = 5)
DimPlot(Tumor, reduction = "umap",group.by = "genotype_celltype", label = F,label.size = 5)
DimPlot(Tumor, reduction = "umap",group.by = "genotype_celltype", label = T,label.size = 5)

table(Tumor$genotype_celltype)
my_levels <- c("EKP_Tumor_1","EKP_Tumor_2","EKP_Tumor_3",
               "RKP_Tumor_1","RKP_Tumor_2","RKP_Tumor_3",
               "T cell","Myeloid","Fibroblast","Endothelial","Others_1","Others_2")
Tumor$genotype_celltype <- factor(x = Tumor$genotype_celltype, levels = my_levels)

color_list = c("#23BD3D","#5DDEB1","#096B08",
               "#2FBDEC","#346CF0","#4136EE",
               "#FF0408","#F2E827","#9F6668","#FF7400","#A9A7A9","#838685")
DimPlot(Tumor, reduction = "umap",group.by = "genotype_celltype",label = F,label.size = 5, cols = color_list)



saveRDS(Tumor, "EKP_RKP_tumor.rds")
Tumor = readRDS("EKP_RKP_tumor.rds")

Immn <- RunPCA(Immn, npcs = 30, features = rownames(Tumor))
Immn <- RunUMAP(Immn, dims = 1:30)
Immn <- FindNeighbors(Immn, reduction = "pca", dims = 1:30)
Immn <- FindClusters(Immn, resolution = 0.3)
DimPlot(Immn, reduction = "umap", label = F,label.size = 5)
DimPlot(Immn, reduction = "umap",group.by = "genotype_celltype", label = F,label.size = 5)
DimPlot(Immn, reduction = "umap",group.by = "genotype_celltype", label = T,label.size = 5)

table(Idents(Immn))
table(Immn$genotype_celltype)
Idents(Immn)= Immn$genotype_celltype


color_list = c("#FFBD22","#FF39BC")
DimPlot(Immn, reduction = "umap",group.by = "genotype_celltype",label = F,label.size = 5, cols = color_list)





saveRDS(Immn, "EKP_RKP_immune.rds")
#Immn = readRDS("EKP_RKP_immune.rds")

########################################### Find markers ######################################################
########################################### Failed ######################################################



DotPlot(Tumor, features = c("Mki67",'Top2a','Stmn1','Aqp5','Lgr5',"Tff2","Muc1","Muc4"), dot.scale = 8, group.by = "genotype", cols = "RdYlBu")+coord_flip()+theme(axis.text.x = element_text(angle = 45, hjust=1))
DotPlot(Tumor, features = c("Mki67",'Top2a','Stmn1','Aqp5','Lgr5',"Tff2","Muc1","Muc4"), dot.scale = 8, group.by = "genotype_celltype", cols = "RdYlBu")+coord_flip()+theme(axis.text.x = element_text(angle = 45, hjust=1))

VlnPlot(Tumor, features = c("Mki67",'Top2a','Stmn1','Aqp5','Lgr5',"Tff2","Muc1","Muc4"), group.by = "celltype", split.by = "genotype", 
        cols = c("ivory", "hotpink"))


FeaturePlot(Tumor, features = c("Mki67",'Top2a'), split.by = "genotype", max.cutoff = 5, 
            min.cutoff = 1, cols = c("lightcyan2", "red1"))
FeaturePlot(Tumor, features = c("Aqp5",'Tff2'), split.by = "genotype", max.cutoff = 5, 
            min.cutoff = 0.5, cols = c("lightcyan2", "red1"))
FeaturePlot(Tumor, features = c("Muc1",'Muc4'), split.by = "genotype", max.cutoff = 5, 
            min.cutoff = 1, cols = c("lightcyan2", "red1"))

head(Tumor@assays$SCT@misc)
all_genes <- rownames(Tumor@assays$SCT@data)
DoHeatmap(Tumor, group.by = "genotype")


DefaultAssay(Tumor)
PrepSCTFindMarkers(Tumor, assay="SCT")
Tumor_DEG <- FindAllMarkers(Tumor)
head(Tumor_DEG)

######################################################################################
########################## AddModuleScore ###########################################

DefaultAssay(Tumor)

wnt = list(c("Axin2","Tbx3"))
tgfb = list(c("Serpine1","Smad7"))
hh = list(c("Gli","Ptch1"))
hippo = list(c("Ccn2","Ccn1"))
bmp = list(c("Id1","Smad6"))
notch = list(c("Hey1","Hey2","Hes1"))

APP_wnt <- AddModuleScore(object = Tumor, features = wnt, name = "wnt", assay = "SCT", nbin=4)
head(APP_wnt)
DotPlot(APP_wnt, features = "wnt1", cols = "RdYlBu",scale.min = 0,scale.max = 20, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")

APP_tgfb <- AddModuleScore(object = Tumor, features = tgfb, name = "tgfb", assay = "SCT", nbin=4)
head(APP_tgfb)
DotPlot(APP_tgfb, features = "tgfb1", cols = "RdYlBu",scale.min = 0,scale.max = 20, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")

APP_hh <- AddModuleScore(object = Tumor, features = hh, name = "hh", assay = "SCT", nbin=4)
head(APP_hh)
DotPlot(APP_hh, features = "hh1", cols = "RdYlBu",scale.min = 0,scale.max = 20, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")

APP_bmp <- AddModuleScore(object = Tumor, features = bmp, name = "bmp", assay = "SCT", nbin=4)
head(APP_bmp)
DotPlot(APP_bmp, features = "bmp1", cols = "RdYlBu",scale.min = 0,scale.max = 20, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")

APP_notch <- AddModuleScore(object = Tumor, features = notch, name = "notch", assay = "SCT", nbin=4)
head(APP_notch)
DotPlot(APP_notch, features = "notch1", cols = "RdYlBu",scale.min = 0,scale.max = 20, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")

####################################################################################################################################################################################


head(Tumor@meta.data)
Idents(Tumor) <- "genotype_celltype"
table(Idents(Tumor))


RKP_tumor <- subset(x= Tumor, idents = c("RKP_Tumor_1","RKP_Tumor_2","RKP_Tumor_3"))
EKP_tumor <- subset(x= Tumor, idents = c("EKP_Tumor_1","EKP_Tumor_2","EKP_Tumor_3"))


ImageFeaturePlot(RKP_tumor, fov="fov_RKP", features = c("Cd80", "Cd86", "Ctla4"), max.cutoff = c(0.1,0.1,0.1), size = 1.5, cols = c("white", "red"))
ImageFeaturePlot(EKP_tumor, fov="fov", features = c("Cd80", "Cd86", "Ctla4"), max.cutoff = c(0.5,0.5,0.5), size = 1.5, cols = c("white", "red"))

Idents(Immn) <- 'genotype'
table(Idents(Immn))


RKP_immn <- subset(x= Immn, idents = c("RKP"))
EKP_immn <- subset(x= Immn, idents = c("EKP"))

ImageFeaturePlot(RKP_immn, fov="fov_RKP", features = c("Ctla4","Pdcd1"), max.cutoff = c(0.5,0.5), size = 1.5, cols = c("white", "red"))
ImageFeaturePlot(EKP_immn, fov="fov", features = c("Ctla4","Pdcd1"), max.cutoff = c(0.5,0.5), size = 1.5, cols = c("white", "red"))


DotPlot(Tumor, features = c("Cd80",'Cd86','Ctla4'), dot.scale = 8, group.by = "genotype", split.by = "genotype_celltype", cols = "RdYlBu")+coord_flip()+theme(axis.text.x = element_text(angle = 45, hjust=1))
DotPlot(Immn, features = c("Cd80",'Cd86','Ctla4'), dot.scale = 8, group.by = "genotype", split.by = "genotype_celltype", cols = "RdYlBu")+coord_flip()+theme(axis.text.x = element_text(angle = 45, hjust=1))
DotPlot(Immn, features = c("Pdcd1"), dot.scale = 8, group.by = "genotype", split.by = "genotype_celltype", cols = "RdYlBu")+coord_flip()+theme(axis.text.x = element_text(angle = 45, hjust=1))


head(Immn@meta.data)
Idents(Immn) <- "genotype"
table(Idents(Immn))

RKP_immn <- subset(x= Immn, idents = c("RKP"))
EKP_immn <- subset(x= Immn, idents = c("EKP"))

ImageFeaturePlot(RKP_immn, fov="fov_RKP", features = c("Cd80", "Cd86", "Ctla4"), max.cutoff = c(3,3,0.1), size = 3, border.color="#F8F5F9",border.size = NA, cols = c("#D3D3D300", "#EF02FA"))
ImageFeaturePlot(EKP_immn, fov="fov", features = c("Cd80", "Cd86", "Ctla4"), max.cutoff = c(3,3,0.1), size = 3, border.color="#F8F5F9",border.size = NA, cols = c("#D3D3D300", "#EF02FA"))

DotPlot(Tumor, features = c("Aqp5",'Runx1','Wfdc2','Tff2'), dot.scale = 8, group.by = "genotype", split.by = "genotype_celltype", cols = "RdYlBu")+coord_flip()+theme(axis.text.x = element_text(angle = 45, hjust=1))
DotPlot(Tumor, features = c("Aqp5",'Runx1','Wfdc2','Tff2'), dot.scale = 8, group.by = "genotype", cols = "RdYlBu")+coord_flip()+theme(axis.text.x = element_text(angle = 45, hjust=1))
DotPlot(Tumor, features = c("Serpine1",'Smad7'), dot.scale = 8, group.by = "genotype", split.by = "genotype", cols = "RdYlBu")+coord_flip()+theme(axis.text.x = element_text(angle = 45, hjust=1))


ImageFeaturePlot(EKP_tumor, fov="fov", features = c("Muc1", "Muc4", "Aqp5"), max.cutoff = c(3,3, 2.5), size = 2, border.color="#F8F5F9",border.size = NA, cols = c("#D3D3D3", "#EF02FA"))                                                                                          
ImageFeaturePlot(RKP_tumor, fov="fov_RKP", features = c("Muc1", "Muc4", "Aqp5"), max.cutoff = c(3,3, 2.5), size = 2, border.color="#F8F5F9",border.size = NA, cols = c("#D3D3D3", "#EF02FA"))                                                                                          


VlnPlot(Immn, features = c("Cd80",'Cd86','Ctla4'),pt.size = 0, group.by = "genotype_celltype", split.by = "genotype", 
        cols = c("ivory", "hotpink"))
VlnPlot(Immn, features = c("Pdcd1"),pt.size = 0, group.by = "genotype_celltype", split.by = "genotype", 
        cols = c("ivory", "hotpink"))
VlnPlot(Immn, features = c("Prf1","Gzmb"),pt.size = 0, group.by = "genotype_celltype", split.by = "genotype", 
        cols = c("ivory", "hotpink"))

VlnPlot(Tumor, features = c("Serpine1",'Smad7','Smad6'),pt.size = 0, group.by = "genotype_celltype", split.by = "genotype", 
        cols = c("ivory", "hotpink"))
VlnPlot(Tumor, features = c("Serpine1",'Smad7'),pt.size = 0, group.by = "genotype", split.by = "genotype", 
        cols = c("green", "blue"))


table(Immn@meta.data$genotype_celltype)




############################### pValue (t-test) ####################################
#################### Tumor cells ############################
ekp_tumor.mat = EKP_tumor@assays$SCT@counts
ekp.cells = colnames(EKP_tumor@assays$SCT@counts)
ekp.genes = rownames(EKP_tumor@assays$SCT@counts)
rownames(ekp_tumor.mat) = ekp.genes
colnames(ekp_tumor.mat) = paste0("EKP_", ekp.cells)

rkp_tumor.mat = RKP_tumor@assays$SCT@counts
rkp.cells = colnames(RKP_tumor@assays$SCT@counts)
rkp.genes = rownames(RKP_tumor@assays$SCT@counts)
rownames(rkp_tumor.mat) = rkp.genes
colnames(rkp_tumor.mat) = paste0("RKP_", rkp.cells)

output.df = data.frame()
for(i in 1:length(rkp.genes)) {
  cat(i, "\t", length(rkp.genes), "\n")
  tmp.gene = rkp.genes[i]
  ekp.array = as.numeric(ekp_tumor.mat[tmp.gene,])
  rkp.array = as.numeric(rkp_tumor.mat[tmp.gene,])
  test.result = t.test(ekp.array, rkp.array)
  pval = test.result$p.value
  output.df[i,1] = tmp.gene
  output.df[i,2] = pval
  
  if(mean(ekp.array) > mean(rkp.array)) {
    output.df[i,3] = "EKP"
  }
  if(mean(ekp.array) < mean(rkp.array)) {
    output.df[i,3] = "RKP"
  }
}

colnames(output.df) = c("Gene", "pval", "Group")
output.df[which(output.df$Gene == "Serpine1"),]
output.df[which(output.df$Gene == "Smad7"),]


#################### all immune cells ############################
ekp_imn.mat = EKP_immn@assays$SCT@counts
ekp.cells = colnames(EKP_immn@assays$SCT@counts)
ekp.genes = rownames(EKP_immn@assays$SCT@counts)
rownames(ekp_imn.mat) = ekp.genes
colnames(ekp_imn.mat) = paste0("EKP_", ekp.cells)

rkp_imn.mat = RKP_immn@assays$SCT@counts
rkp.cells = colnames(RKP_immn@assays$SCT@counts)
rkp.genes = rownames(RKP_immn@assays$SCT@counts)
rownames(rkp_imn.mat) = rkp.genes
colnames(rkp_imn.mat) = paste0("RKP_", rkp.cells)

output.df = data.frame()
for(i in 1:length(rkp.genes)) {
  cat(i, "\t", length(rkp.genes), "\n")
  tmp.gene = rkp.genes[i]
  ekp.array = as.numeric(ekp_imn.mat[tmp.gene,])
  rkp.array = as.numeric(rkp_imn.mat[tmp.gene,])
  test.result = t.test(ekp.array, rkp.array)
  pval = test.result$p.value
  output.df[i,1] = tmp.gene
  output.df[i,2] = pval
  
  if(mean(ekp.array) > mean(rkp.array)) {
    output.df[i,3] = "EKP"
  }
  if(mean(ekp.array) < mean(rkp.array)) {
    output.df[i,3] = "RKP"
  }
}

colnames(output.df) = c("Gene", "pval", "Group")
output.df[which(output.df$Gene == "Cd80"),]
output.df[which(output.df$Gene == "Cd86"),]
output.df[which(output.df$Gene == "Ctla4"),]
output.df[which(output.df$Gene == "Pdcd1"),]


#################### Myeloid cells ############################
table(EKP_immn$genotype_celltype)
Idents(EKP_immn) = "genotype_celltype"
Idents(RKP_immn) = "genotype_celltype"
ekp_my = subset(EKP_immn, idents = c("Myeloid"))
rkp_my = subset(RKP_immn, idents = c("Myeloid"))


ekp_imn.mat = ekp_my@assays$SCT@counts
ekp.cells = colnames(ekp_my@assays$SCT@counts)
ekp.genes = rownames(ekp_my@assays$SCT@counts)
rownames(ekp_imn.mat) = ekp.genes
colnames(ekp_imn.mat) = paste0("EKP_", ekp.cells)

rkp_imn.mat = rkp_my@assays$SCT@counts
rkp.cells = colnames(rkp_my@assays$SCT@counts)
rkp.genes = rownames(rkp_my@assays$SCT@counts)
rownames(rkp_imn.mat) = rkp.genes
colnames(rkp_imn.mat) = paste0("RKP_", rkp.cells)

output.df = data.frame()
for(i in 1:length(rkp.genes)) {
  cat(i, "\t", length(rkp.genes), "\n")
  tmp.gene = rkp.genes[i]
  ekp.array = as.numeric(ekp_imn.mat[tmp.gene,])
  rkp.array = as.numeric(rkp_imn.mat[tmp.gene,])
  test.result = t.test(ekp.array, rkp.array)
  pval = test.result$p.value
  output.df[i,1] = tmp.gene
  output.df[i,2] = pval
  
  if(mean(ekp.array) > mean(rkp.array)) {
    output.df[i,3] = "EKP"
  }
  if(mean(ekp.array) < mean(rkp.array)) {
    output.df[i,3] = "RKP"
  }
}

colnames(output.df) = c("Gene", "pval", "Group")
output.df[which(output.df$Gene == "Cd80"),]
output.df[which(output.df$Gene == "Cd86"),]
output.df[which(output.df$Gene == "Ctla4"),]

#################### T cells ############################
table(EKP_immn$genotype_celltype)
Idents(EKP_immn) = "genotype_celltype"
Idents(RKP_immn) = "genotype_celltype"
ekp_tcell = subset(EKP_immn, idents = c("T cell"))
rkp_tcell = subset(RKP_immn, idents = c("T cell"))

ImageFeaturePlot(rkp_tcell, fov="fov_RKP", features = c("Ctla4","Pdcd1"), max.cutoff = c(0.5,0.5), size = 3, cols = c("white", "red"))
ImageFeaturePlot(ekp_tcell, fov="fov", features = c("Ctla4","Pdcd1"), max.cutoff = c(0.5,0.5), size = 3, cols = c("white", "red"))



ekp_imn.mat = ekp_tcell@assays$SCT@counts
ekp.cells = colnames(ekp_tcell@assays$SCT@counts)
ekp.genes = rownames(ekp_tcell@assays$SCT@counts)
rownames(ekp_imn.mat) = ekp.genes
colnames(ekp_imn.mat) = paste0("EKP_", ekp.cells)

rkp_imn.mat = rkp_tcell@assays$SCT@counts
rkp.cells = colnames(rkp_tcell@assays$SCT@counts)
rkp.genes = rownames(rkp_tcell@assays$SCT@counts)
rownames(rkp_imn.mat) = rkp.genes
colnames(rkp_imn.mat) = paste0("RKP_", rkp.cells)

output.df = data.frame()
for(i in 1:length(rkp.genes)) {
  cat(i, "\t", length(rkp.genes), "\n")
  tmp.gene = rkp.genes[i]
  ekp.array = as.numeric(ekp_imn.mat[tmp.gene,])
  rkp.array = as.numeric(rkp_imn.mat[tmp.gene,])
  test.result = t.test(ekp.array, rkp.array)
  pval = test.result$p.value
  output.df[i,1] = tmp.gene
  output.df[i,2] = pval
  
  if(mean(ekp.array) > mean(rkp.array)) {
    output.df[i,3] = "EKP"
  }
  if(mean(ekp.array) < mean(rkp.array)) {
    output.df[i,3] = "RKP"
  }
}

colnames(output.df) = c("Gene", "pval", "Group")
output.df[which(output.df$Gene == "Cd80"),]
output.df[which(output.df$Gene == "Cd86"),]
output.df[which(output.df$Gene == "Ctla4"),]

