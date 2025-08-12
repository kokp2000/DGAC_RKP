rm(list = ls())
##########################################################################################
####################merge ekp and rkp
##########################################################################################
library(reticulate)
reticulate::py_install(packages = 'umap-learn')
reticulate::py_install(packages= "numpy")

library(CellChat)
library(ggplot2)                  
library(patchwork)
library(igraph)

setwd('D:/KP/rkp_RKP')

ekp <- readRDS('ekp_after_cellchat.rds') 
rkp <- readRDS('rkp_after_cellchat.rds')



levels(ekp@idents)
levels(rkp@idents)


#Lift up CellChat object and merge together
#Since there are additional two populations (i.e., dermal DC and pericytes) specific to E14.5 
#compared to E13.5, we lift up cellchat.E13 by lifting up the cell groups to the same cell 
#labels as E14.5. liftCellChat will only update the slot related to the cell-cell 
#communication network, including slots object@net, object@netP and object@idents.

# Define the cell labels to lift up
group.new = levels(ekp@idents)
ekp <- liftCellChat(ekp, group.new)
rkp <- liftCellChat(rkp, group.new)


saveRDS(ekp, "ekp_for_merge.rds")
saveRDS(rkp, "rkp_for_merge.rds")



#> The CellChat object will be lifted up using the cell labels FIB-A, FIB-B, FIB-P, DC, Pericyte, MYL, Immune, ENDO, Muscle, MELA, Basal-P, Basal, Spinious
#> Update slots object@net, object@netP, object@idents in a single dataset...

object.list_12 <- list(ekp = ekp, rkp = rkp)
cellchat_12 <- mergeCellChat(object.list_12, add.names = names(object.list_12), cell.prefix = TRUE)
cellchat_12
saveRDS(cellchat_12, 'cellchat.ekp_rkp_merged.rds')


###############################################################
#Compare the total number of interactions and interaction strength
#To answer on question on whether the cell-cell communication is enhanced or not, 
#CellChat compares the the total number of interactions and interaction strength of 
#the inferred cell-cell communication networks from different biological conditions.
gg1 <- compareInteractions(cellchat_12, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_12, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#To identify the interaction between which cell populations showing significant changes, 
#CellChat compares the number of interactions and interaction strength among different cell populations.

#The differential number of interactions or interaction strength in the cell-cell communication network 
#between two datasets can be visualized using circle plot, where red (or blue) colored edges represent 
#increased (or decreased) signaling in the second dataset compared to the first one.

# red: increased in rkp, blue: decreased in rkp
par(mfrow = c(1,1), xpd=TRUE)

netVisual_diffInteraction <- function(cellchat_12, comparison = c(1,2), 
                                      measure = c("count", "weight", "count.merged", "weight.merged"), 
                                      color.use = NULL, 
                                      color.edge = c('#b2182b','#2166ac'), 
                                      title.name = NULL, 
                                      sources.use = NULL, 
                                      targets.use = NULL, 
                                      remove.isolate = FALSE, 
                                      top = 1,
                                      weight.scale = FALSE, 
                                      vertex.weight = 20, 
                                      vertex.weight.max = NULL, 
                                      vertex.size.max = 15, 
                                      vertex.label.cex=1,
                                      vertex.label.color= "black",
                                      edge.weight.max = NULL, 
                                      edge.width.max=8, 
                                      alpha.edge = 0.6, 
                                      label.edge = FALSE,
                                      edge.label.color='black',
                                      edge.label.cex=0.8,
                                      edge.curved=0.2,
                                      shape='circle',
                                      layout=in_circle(), 
                                      margin=0.2,
                                      arrow.width=1,
                                      arrow.size = 0.2
){
  options(warn = -1)
  measure <- match.arg(measure)
  
  obj1_raw <- cellchat_12@net[[comparison[1]]][[measure]]
  obj2_raw <- cellchat_12@net[[comparison[2]]][[measure]]
  
  shared_label = intersect(rownames(obj1_raw), rownames(obj2_raw))
  obj1 = obj1_raw[shared_label, shared_label]
  obj2 = obj2_raw[shared_label, shared_label]
  
  net.diff <- obj2 - obj1
  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  } else if (measure %in% c("weight", "weight.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
    net[is.na(net)] <- 0
  }
  
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  
  net[abs(net) < stats::quantile(abs(net), probs = 1-top, na.rm= T)] <- 0
  
  g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5
  
  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  
  # ERROR: Length of new attribute value must be 1 or 45
  # igraph::V(g)$size<-vertex.weight
  
  # Fix below
  igraph::V(g)$size<-vertex.weight
  
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  #igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1],color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, alpha.edge)
  
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }
  
  
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    # add this line
    igraph::E(g)$loop.angle = rep(0, nrow(edge.start))
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+2
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
  if (!is.null(title.name)) {
    text(0,1.5,title.name, cex = 1.1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}

netVisual_diffInteraction(cellchat_12, weight.scale = T)
netVisual_diffInteraction(cellchat_12, weight.scale = T, measure = "weight")







#We can also show differential number of interactions or interaction strength in a greater details using 
#a heatmap. The top colored bar plot represents the sum of column of values displayed in the heatmap 
#(incoming signaling). The right colored bar plot represents the sum of row of values (outgoing signaling). 
#In the colorbar, red (or blue) represents increased (or decreased) signaling in the second dataset compared 
#to the first one.
######### modified_source code

netVisual_heatmap <- function(cellchat_12, comparison = c(1,2), measure = c("count", "weight"), signaling = NULL, slot.name = c("netP", "net"), color.use = NULL, color.heatmap = c("#2166ac","#b2182b"),
                              title.name = NULL, width = NULL, height = NULL, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE,
                              sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, row.show = NULL, col.show = NULL){
  # obj1 <- cellchat_12.list[[comparison[1]]]
  # obj2 <- cellchat_12.list[[comparison[2]]]
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(cellchat_12@net[[1]])) {
    message("Do heatmap based on a merged cellchat_12 \n")
    obj1 <- cellchat_12@net[[comparison[1]]][[measure]]
    obj2 <- cellchat_12@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    
    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Differential number of interactions"
      }
    } else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Differential interaction strength"
      }
    }
    legend.name = "Relative values"
  } else {
    message("Do heatmap based on a single cellchat_12 \n")
    if (!is.null(signaling)) {
      net.diff <- slot(cellchat_12, slot.name)$prob[,,signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    } else if (!is.null(measure)) {
      net.diff <- cellchat_12@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      } else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }
  
  net <- net.diff
  
  
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    if (length(idx) > 0) {
      net <- net[-idx, ]
      net <- net[, -idx]
    }
  }
  
  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[ ,col.show]
    color.use <- color.use[col.show]
  }
  
  
  if (min(mat) < 0) {
    color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)), c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), 0, round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
    # color.heatmap.use = colorRamp3(c(seq(min(mat), -(max(mat)-min(max(mat)))/9, length.out = 4), 0, seq((max(mat)-min(max(mat)))/9, max(mat), length.out = 4)), RColorBrewer::brewer.pal(n = 9, name = color.heatmap))
  } else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)), color.heatmap)
    } else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min(mat), max(mat)), color.heatmap)
    } else if (length(color.heatmap) == 1) {
      color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
  }
  # col_fun(as.vector(mat))
  
  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  } else {
    mat[mat == 0] <- NA
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = legend.name,
                bottom_annotation = col_annotation, left_annotation =row_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows,cluster_columns = cluster.rows,
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                # width = unit(width, "cm"), height = unit(height, "cm"),
                column_title = title.name,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
                row_title = "Sources (Sender)",row_title_gp = gpar(fontsize = font.size.title),row_title_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                            border = NA, #at = colorbar.break,
                                            legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
  )
  #  draw(ht1)
  return(ht1)
}



######## The differential network analysis only works for pairwise datasets. 

# red: increased in rkp, blue: decreased in rkp
gg1 <- netVisual_heatmap(cellchat_12)
#> Do heatmap based on a merged cellchat_12
gg2 <- netVisual_heatmap(cellchat_12, measure = "weight")
#> Do heatmap based on a merged cellchat_12
gg1 + gg2


#To better control the node size and edge weights of the inferred networks across different datasets, we compute 
#the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) 
#across all datasets.
weight.max <- getMaxWeight(object.list_12, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list_12)) {
  netVisual_circle(object.list_12[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list_12)[i]))
}

#Identify signaling groups based on their functional similarit
library(umap)
cellchat_12 <- computeNetSimilarityPairwise(cellchat_12, type = "functional")
cellchat_12 <- netEmbedding(cellchat_12, type = "functional")
cellchat_12 <- netClustering(cellchat_12, type = "functional")
netVisual_embeddingPairwise(cellchat_12, type = "functional", label.size = 3.5)
#> Compute the distance of signaling networks between datasets 1 2
gg1 <- rankNet(cellchat_12, mode = "comparison", stacked = T, do.stat = FALSE)
gg2 <- rankNet(cellchat_12, mode = "comparison", stacked = F, do.stat = FALSE)
gg1 + gg2


table(cellchat_12@idents$ekp)
table(cellchat_12@idents$rkp)


#Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
group.new
netVisual_bubble(cellchat_12, sources.use = 11, targets.use = c(5:7,9),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat_12, sources.use = 11, targets.use = c(2,3,10,11,13), comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat_12, sources.use = c(5,6,7,9), targets.use = c(2,3,10,12), comparison = c(1, 2), angle.x = 45)

#> Comparing communications on a merged object

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "rkp"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat_12_1 <- identifyOverExpressedGenes(cellchat_12, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat_12_1, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat_12_1, net = net, datasets = "rkp",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat_12_1, net = net, datasets = "ekp",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat_12_1)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_12_1)

#We then visualize the upgulated and down-regulated signaling ligand-receptor pairs using bubble plot or chord diagram

#   [1] "B cell"           "CD4 T cell"       "CD8 T cell"       "Fibroblast"       "M1 Macrophage"    "M2 Macrophage"    "Macrophage"       "Mast cell"       
#   [9] "Monocyte"         "T cell"           "effector T cell"  "epithelial"       "exhausted T cell"

#c(12), targets.use = c(5:7,9)
# upregulated signaling
pairLR.use.up = net.up[, "interaction_name", drop = F]
# epi to Myeloid
gg1 <- netVisual_bubble(cellchat_12_1, pairLR.use = pairLR.use.up, sources.use = c(11), targets.use = c(5:7,9), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("upregulated in ", names(object.list_12)[1]))
gg1
# T to Myeloid  
gg1 <- netVisual_bubble(cellchat_12_1, pairLR.use = pairLR.use.up, sources.use = c(2,3,10,12), targets.use = c(5:7,9,11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("upregulated in ", names(object.list_12)[1]))
gg1
# Myeloid to T
gg1 <- netVisual_bubble(cellchat_12_1, pairLR.use = pairLR.use.up, sources.use = c(5:7,9), targets.use = c(2,3,10,11,12), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("upregulated in ", names(object.list_12)[1]))
gg1






























