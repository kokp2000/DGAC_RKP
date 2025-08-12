rm(list = ls())

phate.file = "D:/KP/EKP_RKP/phate/PHATE.Default.output.meta.txt"
ccat.file = "D:/KP/EKP_RKP/phate/CCAT.txt"
dist.file = "D:/KP/EKP_RKP/phate/graphtools.txt"
#aci.scvelo.file = "/mnt/e/JH/Research/COLON-CRACD_2024-04-01/workspace/scanpy/ACI_scvelo.txt"
#wi.scvelo.file = "/mnt/e/JH/Research/COLON-CRACD_2024-04-01/workspace/scanpy/WI_scvelo.txt"
meta.df = read.table("D:/KP/EKP_RKP/files_for_phate/EKP_RKP_metaData.csv", header = T, as.is = T, sep = ",")
rownames(meta.df) = meta.df[,1]

phate.df = read.delim(phate.file, header = T, as.is = T)
ccat.df = read.delim(ccat.file, header = T, as.is = T)
rownames(ccat.df) = ccat.df[,1]
dist.df = read.csv(dist.file, header=  T, as.is = T, sep = ",")
rownames(dist.df) = dist.df[,1]
#aci.scvelo.df = read.table(aci.scvelo.file, header = T, as.is = T, sep = ",", row.names = 1)
#rownames(aci.scvelo.df) = paste0(rownames(aci.scvelo.df), "-0")
#wi.scvelo.df = read.table(wi.scvelo.file, header = T, as.is = T, sep = ",", row.names = 1)
#rownames(wi.scvelo.df) = paste0(rownames(wi.scvelo.df), "-1")
#scvelo.df = rbind(aci.scvelo.df, wi.scvelo.df)

inter.cells = intersect(rownames(phate.df), rownames(ccat.df))
phate.df[1:5,]



inter.cells = intersect(inter.cells, rownames(meta.df))
inter.cells[1:10]

output.df = data.frame(PHATE1 = phate.df[inter.cells,"PHATE1"], PHATE2 = phate.df[inter.cells,"PHATE2"], Dist = dist.df[inter.cells,'dist_scaled'], CCAT = ccat.df[inter.cells,"CCAT"], VelocityLength = meta.df[inter.cells,"velocity_length"], celltype = meta.df[inter.cells,'celltype'], type = meta.df[inter.cells,'type'])
rownames(output.df) = inter.cells
output.df$Scaled_VelocityLength = (output.df$VelocityLength)^(-1)/max((output.df$VelocityLength)^(-1))
cluster.array = unique(output.df$celltype)
for(i in 1:length(cluster.array)) {
  tmp.cluster = cluster.array[i]
  tmp.cells = rownames(output.df)[which(output.df$celltype == tmp.cluster)]
  output.df[tmp.cells,"Med_CCAT"] = median(output.df[tmp.cells,"CCAT"])
  output.df[tmp.cells,"Med_Scaled_VelocityLength"] = median(output.df[tmp.cells,"Scaled_VelocityLength"])
}
output.df$Valley = output.df$Med_CCAT
output.df$Ridge = output.df$Med_Scaled_VelocityLength * output.df$Dist
output.df$VR_score = 0.9*output.df$Valley + 0.1*output.df$Ridge


write.table(output.df, "D:/KP/EKP_RKP/phate/Waddington_Input_PHATE_by_celltype.txt", row.names = T, col.names = T, sep = "\t", quote = F)
write.csv(output.df, "D:/KP/EKP_RKP/phate/Waddington_Input_PHATE_by_celltype.csv")
