rm(list=ls())

install.packages('pak')
pak::pkg_install("r-lib/rlang")
devtools::install_github("igordot/msigdbr")
BiocManager::install("fgsea")

library(data.table)
library(fgsea)
library(ggplot2)
library(BiocParallel)
library(msigdbr)
library(xlsx)
library(dplyr)


setwd("D:/KP/EKP_RKP/fgsea")
######################################################################################################
####################### fgsea using Wilcoxon gene rank list from Scanpy ##############################
######################################################################################################

###################### RKP to EKP ################################################################

### KEGG ###

DEG <- read.csv("D:/KP/EKP_RKP/CD8_Tcell_rank_by_type_fgsea.csv", sep = ',')
head(DEG)
DEG_1 <- DEG %>% select(RKP_names, RKP_score, RKP_pvals)
DEG_1 %>% filter(RKP_pvals < 0.05)
DEG_1$gene <- DEG_1$RKP_names
row.names(DEG_1) <- DEG_1$gene

##### rank with log2FC may not correct ####
#### I used new rank (avg_log2FC * -log10(p.adj)) ####calculate and make new column in Excel ####
ranks_1 <- DEG_1$RKP_score
#ranks_1 <- DEG_1$rank
names(ranks_1) <- row.names(DEG_1)
head(ranks_1)
fgsea.p1<-barplot(sort(ranks_1, decreasing = TRUE))

msigdbr_show_species()
genesets_1 = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
msigdbr_list = split(x = genesets_1$gene_symbol, f = genesets_1$gs_name)
fgseaRes_1 <- fgsea(msigdbr_list, ranks_1, minSize=15, maxSize = 500)
head(fgseaRes_1[order(padj, -abs(NES)), ], n=10)

topUp <- fgseaRes_1 %>% filter(ES > 0) %>% top_n(5, wt=-padj)
topDown <- fgseaRes_1 %>% filter(ES < 0) %>% top_n(5, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% arrange(-ES)
fgsea.p2<-plotGseaTable(msigdbr_list[topPathways$pathway], ranks_1, fgseaRes_1, gseaParam = 0.5)

fgseaResTidy <- fgseaRes_1 %>%  as_tibble() %>% arrange(desc(NES))

fwrite(fgseaResTidy, file="RKP_CD8_Tcell_GSEA_results_KEGG.csv", sep=",", sep2=c("", " ", ""))
fgseaResTidy = read.csv('d:/KP/EKP_RKP/fgsea/Tcell/RKP_CD8_Tcell_GSEA_results_KEGG.csv', sep=",")


ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

plotEnrichment(msigdbr_list[["KEGG_ERBB_SIGNALING_PATHWAY"]],
               ranks_1) + labs(title="KEGG_ERBB_SIGNALING_PATHWAY")





#### REACTOME  ####

genesets_2 = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
msigdbr_list_2 = split(x = genesets_2$gene_symbol, f = genesets_2$gs_name)
msigdbr_list_2
fgseaRes_2 <- fgsea(pathways=msigdbr_list_2, ranks_1, minSize=5, maxSize = 500)
head(fgseaRes_2[order(padj, -abs(NES)), ], n=10)


topUp_2 <- fgseaRes_2 %>% filter(ES > 0) %>% top_n(5, wt=-padj)
topDown_2 <- fgseaRes_2 %>% filter(ES < 0) %>% top_n(5, wt=-padj)
topPathways_2 <- bind_rows(topUp_2, topDown_2) %>% arrange(-ES)

## To remove previoud plots, use this code ## 
##dev.off() 
fgsea.p3<-plotGseaTable(msigdbr_list_2[topPathways_2$pathway], ranks_1, fgseaRes_2, gseaParam = 0.5)



fgseaResTidy_2 <- fgseaRes_2 %>%  as_tibble() %>% arrange(desc(NES))
fgseaResTidy_2 %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

fwrite(fgseaResTidy_2, file="RKP_CD8_Tcell_GSEA_results_REACTOME.csv", sep=",", sep2=c("", " ", ""))
fgseaResTidy_2 = read.csv('d:/KP/EKP_RKP/fgsea/Tcell/RKP_CD8_Tcell_GSEA_results_REACTOME.csv', sep=",")


ggplot(fgseaResTidy_2, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="REACTOME Hallmark pathways") + 
  theme_minimal()

plotEnrichment(msigdbr_list_2[["REACTOME_RUNX3_REGULATES_YAP1_MEDIATED_TRANSCRIPTION"]],
               ranks_1) + labs(title="REACTOME_RUNX3_REGULATES_YAP1_MEDIATED_TRANSCRIPTION")
plotEnrichment(msigdbr_list_2[["REACTOME_RUNX3_REGULATES_IMMUNE_RESPONSE_AND_CELL_MIGRATION"]],
               ranks_1) + labs(title="REACTOME_RUNX3_REGULATES_IMMUNE_RESPONSE_AND_CELL_MIGRATION")
plotEnrichment(msigdbr_list_2[["REACTOME_REGULATION_OF_RUNX3_EXPRESSION_AND_ACTIVITY"]],
               ranks_1) + labs(title="REACTOME_REGULATION_OF_RUNX3_EXPRESSION_AND_ACTIVITY")
plotEnrichment(msigdbr_list_2[["REACTOME_RUNX1_REGULATES_TRANSCRIPTION_OF_GENES_INVOLVED_IN_DIFFERENTIATION_OF_HSCS"]],
               ranks_1) + labs(title="REACTOME_RUNX1_REGULATES_TRANSCRIPTION_OF_GENES_INVOLVED_IN_DIFFERENTIATION_OF_HSCS")
plotEnrichment(msigdbr_list_2[["REACTOME_REGULATION_OF_RUNX2_EXPRESSION_AND_ACTIVITY"]],
               ranks_1) + labs(title="REACTOME_REGULATION_OF_RUNX2_EXPRESSION_AND_ACTIVITY")

plotEnrichment(msigdbr_list_2[["REACTOME_TRANSCRIPTIONAL_REGULATION_BY_RUNX3"]],
               ranks_1) + labs(title="REACTOME_TRANSCRIPTIONAL_REGULATION_BY_RUNX3")
plotEnrichment(msigdbr_list_2[["REACTOME_TRANSCRIPTIONAL_REGULATION_BY_RUNX1"]],
               ranks_1) + labs(title="REACTOME_TRANSCRIPTIONAL_REGULATION_BY_RUNX1")
plotEnrichment(msigdbr_list_2[["REACTOME_TRANSCRIPTIONAL_REGULATION_BY_RUNX2"]],
               ranks_1) + labs(title="REACTOME_TRANSCRIPTIONAL_REGULATION_BY_RUNX2")
plotEnrichment(msigdbr_list_2[["REACTOME_RUNX1_INTERACTS_WITH_CO_FACTORS_WHOSE_PRECISE_EFFECT_ON_RUNX1_TARGETS_IS_NOT_KNOWN"]],
               ranks_1) + labs(title="REACTOME_RUNX1_INTERACTS_WITH_CO_FACTORS_WHOSE_PRECISE_EFFECT_ON_RUNX1_TARGETS_IS_NOT_KNOWN")

plotEnrichment(msigdbr_list_2[["REACTOME_SIGNALING_BY_ERBB2"]],
               ranks_1) + labs(title="REACTOME_SIGNALING_BY_ERBB2")
plotEnrichment(msigdbr_list_2[["REACTOME_MAPK_FAMILY_SIGNALING_CASCADES"]],
               ranks_1) + labs(title="REACTOME_MAPK_FAMILY_SIGNALING_CASCADES")


####### GOBP
genesets_5 = msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
msigdbr_list_5 = split(x = genesets_5$gene_symbol, f = genesets_5$gs_name)
msigdbr_list_5
fgseaRes_5 <- fgsea(pathways=msigdbr_list_5, ranks_1, minSize=5, maxSize = 500)
head(fgseaRes_5[order(padj, -abs(NES)), ], n=10)


topUp_5 <- fgseaRes_5 %>% filter(ES > 0) %>% top_n(5, wt=-padj)
topDown_5 <- fgseaRes_5 %>% filter(ES < 0) %>% top_n(5, wt=-padj)
topPathways_5 <- bind_rows(topUp_5, topDown_5) %>% arrange(-ES)

## To remove previoud plots, use this code ## 
dev.off() 
fgsea.p5<-plotGseaTable(msigdbr_list_5[topPathways_5$pathway], ranks_1, fgseaRes_5, gseaParam = 0.5)



fgseaResTidy_5 <- fgseaRes_5 %>%  as_tibble() %>% arrange(desc(NES))
fgseaResTidy_5 %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

fwrite(fgseaResTidy_5, file="RKP_CD8_Tcell_GSEA_results_GOBP.csv", sep=",", sep2=c("", " ", ""))
fgseaTidy_5 <- read.csv('d:/KP/EKP_RKP/fgsea/Tcell/RKP_CD8_Tcell_GSEA_results_GOBP.csv',sep=',')

plotEnrichment(msigdbr_list_5[["GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_6_PRODUCTION"]],
               ranks_1) + labs(title="GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_6_PRODUCTION")
plotEnrichment(msigdbr_list_5[["GOBP_SMAD_PROTEIN_SIGNAL_TRANSDUCTION"]],
               ranks_1) + labs(title="GOBP_SMAD_PROTEIN_SIGNAL_TRANSDUCTION")



ggplot(fgseaResTidy_5, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="RKPvs.others GOBP pathways") + 
  theme_minimal()
