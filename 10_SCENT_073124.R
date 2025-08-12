rm(list = ls())
library(SCENT)
library(data.table)

count.df = as.data.frame(fread("D:/KP/EKP_RKP/EKP_RKP_raw_counts.csv"))
count.df_1 = count.df
rownames(count.df_1) = count.df_1[,1]
count.df_1 = count.df_1[, -1]
count.df = count.df_1

count.df_2 = as.data.frame(t(count.df)) # column: cell ID, rownames: gene name
count.df_2[1:5,1:5]

gene.df = read.table("D:/KP/EKP_RKP/EKP_RKP_genes.csv", header = T, as.is = T, sep = ",")
meta.df = read.table("D:/KP/EKP_RKP/EKP_RKP_metadata.csv", header = T, as.is = T, sep = ",")
rownames(meta.df) = meta.df[,1]
rownames(gene.df) = gene.df[,2]

output.dir = "D:/KP/EKP_RKP/phate"


# Normlize to Count Per Million
cpm.df = sapply(as.data.frame(count.df_2), function(x) {1000000 * (x/sum(x))})
rownames(cpm.df) = rownames(count.df_2)

cpm.df[1:5,1:5]

# Gene symbol to entrez id
library(org.Mm.eg.db)
Mm <- org.Mm.eg.db
my.symbols <- rownames(cpm.df)
my.df = select(Mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
b.index = which(is.na(my.df[,2]))
if(length(b.index) > 0) {
  id.cpm.df = cpm.df[-b.index,]
  rownames(id.cpm.df) = my.df[-b.index,2]
}

id.cpm.df[1:5, 1:5]

# Convert to HUMAN
or.file = "D:/KP/gene_orthologs" # Obtain from NCBI "https://ftp.ncbi.nlm.nih.gov/gene/DATA/" , download 'gene_orthologs.gz'
or.df = as.data.frame(fread(or.file))
colnames(or.df)[1] = 'tax_id'
hm.or.df = subset(or.df, Other_tax_id == "10090" & tax_id == "9606")

id.array = rownames(id.cpm.df)
hm.array = sapply(id.array, function(x){return(subset(hm.or.df, Other_GeneID == x)[1,"GeneID"])})

b.index = which(is.na(hm.array))
if(length(b.index) > 0) {
  f.hm.array = hm.array[-b.index]
}
f.id.cpm.df = id.cpm.df[-b.index,]
rownames(f.id.cpm.df) = f.hm.array

# Run SCENT
data(net13Jun12)
log.f.id.cpm.df <- log2(f.id.cpm.df+1);
ccat.v <- CompCCAT(exp = log.f.id.cpm.df, ppiA = net13Jun12.m)
output.df = data.frame(Cell = colnames(log.f.id.cpm.df), CCAT = ccat.v)
write.table(output.df, paste0(output.dir, "/CCAT.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
