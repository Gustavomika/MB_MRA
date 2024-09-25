library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(limma)
library(RTN)
library(ggpubr)
library(ComplexHeatmap)


# set to your working directory
setwd("C:/Users/guga_/Documents/tayrone")

# GO enrichment analisis of regulon identified in the graphs from "V1_rede_plots" script


load("./Dados/Processados/V1/V1_rede_sub_7/rtni_final_sub_7.rdata")
load("./Dados/Processados/V1/V1_MRA_SHH/MRA_SHH_7.rdata")
SHH = master
load("./Dados/Processados/V1/V1_MRA_G3/MRA_G3_7.rdata")
G3 = master
load("./Dados/Processados/V1/V1_MRA_G4/MRA_G4_7.rdata")
G4 = master

j = c(SHH$Regulon, G3$Regulon, G4$Regulon)
j = unique(j)




## Over representation analysis of regulons identified visualizing "V1_rede_plots" graphs

# Bottom (G4 regulons present in the branch)
phyper(28,104,1477,77,lower.tail= FALSE)
#phyper(37,104,1477,131,lower.tail= FALSE)

# Upper (G3 regulons present in the branch)
phyper(19,87,1494,82,lower.tail= FALSE)
#phyper(23,87,1494,131,lower.tail= FALSE)

# Gold (regulators shared for subgroups SHH, Group 3, and Group 4 in the branch)
phyper(6,21,1560,20,lower.tail= FALSE)


# G4 regulons present in the branch
esquerda_menor = c("CXXC4", "ZNF462", "RARB","POU3F2","SIX6","PBX1","FEZF1","ZSCAN20",
                   "DLX2","POU6F2","VEZF1","LIN28B","KAT7","GLI3","ZNF311","SP110","MSC","ZNF516",
                   "PPARA","TBX18","ZNF334","NEUROG1","OTX2","FOXG1","AKAP8L",
                   "SATB2","RREB1","EOMES","PBX3")

# G3 regulons present in the branch
direita_menor = c("IRX6","PBX4","ZNF385B","PURG","DPF3","ZMAT4","TCF7","INSM2","ASCL1","RFX4","NPAS3","DBX2",
                  "HES6","NACC2","MKX","ZFAT","HIVEP3","ARNT2","PRDM8","TSHZ3","CASZ1")

# regulators shared for subgroups SHH, Group 3, and Group 4 in the branch
dourado_maior = c("ASCL1","RFX4","NPAS3","DBX2","HES6","TBX19","TCF20","MKX","NACC2","MSX2")

####### G4 regulons GOs ######




regulons <- tni.get(rtni_final_sub, what = "regulons.and.mode", idkey = "ID")
regulons = regulons[which(names(regulons) %in% esquerda_menor)] 


i = 1
genes = NULL

while (i <= length(regulons)) {
  
  genes = c(genes, names(regulons[[i]]))
  i = i + 1
  
}

genes = genes[!duplicated(genes)]



ensen_degs = mapIds(org.Hs.eg.db, alias2Symbol(genes),"ENSEMBL", "SYMBOL")

go_enrich <- enrichGO(
  gene = ensen_degs,
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  keyType = 'ENSEMBL'
)



df_go <- go_enrich@result
df_go_esquerda_menor <- df_go[df_go$p.adjust <= 0.01,]
df_go_esquerda_menor = cbind(rep("A", length(df_go_esquerda_menor[,1])), df_go_esquerda_menor)
colnames(df_go_esquerda_menor)[1] = "Cluster"

entrez_degs = mapIds(org.Hs.eg.db, alias2Symbol(genes),"ENTREZID", "SYMBOL")


kegg_enrich <- enrichKEGG(
  gene = entrez_degs,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
)


df_kegg <- kegg_enrich@result
df_kegg_esquerda_menor <- df_kegg[df_kegg$p.adjust < 0.05,]
df_kegg_esquerda_menor = cbind(rep("A", length(df_kegg_esquerda_menor[,1])), df_kegg_esquerda_menor)
colnames(df_kegg_esquerda_menor)[1] = "Cluster"

##### G3 regulons GOs #####


regulons <- tni.get(rtni_final_sub, what = "regulons.and.mode", idkey = "ID")
regulons = regulons[which(names(regulons) %in% direita_menor)] 


i = 1
genes = NULL

while (i <= length(regulons)) {
  
  genes = c(genes, names(regulons[[i]]))
  i = i + 1
  
}

genes = genes[!duplicated(genes)]



ensen_degs = mapIds(org.Hs.eg.db, alias2Symbol(genes),"ENSEMBL", "SYMBOL")

go_enrich <- enrichGO(
  gene = ensen_degs,
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  keyType = 'ENSEMBL'
)



df_go <- go_enrich@result
df_go_direita_menor <- df_go[df_go$p.adjust <= 0.01,]
df_go_direita_menor = cbind(rep("B", length(df_go_direita_menor[,1])), df_go_direita_menor)
colnames(df_go_direita_menor)[1] = "Cluster"

entrez_degs = mapIds(org.Hs.eg.db, alias2Symbol(genes),"ENTREZID", "SYMBOL")


kegg_enrich <- enrichKEGG(
  gene = entrez_degs,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
)


df_kegg <- kegg_enrich@result
df_kegg_direita_menor <- df_kegg[df_kegg$p.adjust < 0.05,]
df_kegg_direita_menor = cbind(rep("B", length(df_kegg_direita_menor[,1])), df_kegg_direita_menor)
colnames(df_kegg_direita_menor)[1] = "Cluster"


########### Shared regulons GOs # #####


regulons <- tni.get(rtni_final_sub, what = "regulons.and.mode", idkey = "ID")
regulons = regulons[which(names(regulons) %in% dourado_maior)] 


i = 1
genes = NULL

while (i <= length(regulons)) {
  
  genes = c(genes, names(regulons[[i]]))
  i = i + 1
  
}

genes = genes[!duplicated(genes)]



ensen_degs = mapIds(org.Hs.eg.db, alias2Symbol(genes),"ENSEMBL", "SYMBOL")

go_enrich <- enrichGO(
  gene = ensen_degs,
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  keyType = 'ENSEMBL'
)



df_go <- go_enrich@result
df_go_dourado_maior <- df_go[df_go$p.adjust <= 0.01,]
df_go_dourado_maior = cbind(rep("C", length(df_go_dourado_maior[,1])), df_go_dourado_maior)
colnames(df_go_dourado_maior)[1] = "Cluster"

entrez_degs = mapIds(org.Hs.eg.db, alias2Symbol(genes),"ENTREZID", "SYMBOL")


kegg_enrich <- enrichKEGG(
  gene = entrez_degs,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
)


df_kegg <- kegg_enrich@result
df_kegg_dourado_maior <- df_kegg[df_kegg$p.adjust < 0.05,]
df_kegg_dourado_maior = cbind(rep("C", length(df_kegg_dourado_maior[,1])), df_kegg_dourado_maior)
colnames(df_kegg_dourado_maior)[1] = "Cluster"


########## plots ############

library(stringr)
toPlot =  rbind(df_go_esquerda_menor, df_go_direita_menor, df_go_dourado_maior)
toPlot$Description <- str_wrap(toPlot$Description, width = 30)
toPlot$Cluster = factor(toPlot$Cluster, levels = c("C", "B", "A"))


## run: 
data(gcSample)
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)


xx@compareClusterResult = toPlot

p = dotplot(xx, showCategory=7, size = "count")+
  scale_y_discrete(guide = guide_axis(n.dodge = 2))+
  scale_y_discrete(guide = guide_axis(angle = 40))+
  coord_flip()+ theme(legend.position="left")


ggexport(filename = "./Dados/Processados/V1/V1_GO_masters/menor.png",
         plot = p , device = "png",
         width = 4500,
         height = 1500, res = 330)

toPlot2 =  rbind(df_go_esquerda, df_go_direita, df_go_dourado_maior)



## run: 
data(gcSample)
xx1 <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)


xx1@compareClusterResult = toPlot2

p1 = dotplot(xx1, showCategory=7, size = "count")


ggexport(filename = "./Dados/Processados/V1/V1_GO_masters/maior.png",
         plot = p1 , device = "png",
         width = 1700,
         height = 2700, res = 330)


toPlot3 =  rbind(df_kegg_esquerda_menor, df_kegg_direita_menor, df_kegg_dourado_maior)


## run: 
data(gcSample)
xx2 <- compareCluster(gcSample, fun="enrichKEGG",
                      organism="hsa", pvalueCutoff=0.05)


xx2@compareClusterResult = toPlot3

p2 = dotplot(xx2, showCategory=7, size = "count")


ggexport(filename = "./Dados/Processados/V1/V1_GO_masters/kegg_menor.png",
         plot = p2 , device = "png",
         width = 1700,
         height = 2700, res = 330)

toPlot4 =  rbind(df_kegg_esquerda, df_kegg_direita_menor, df_kegg_dourado_maior)


## run: 
data(gcSample)
xx3 <- compareCluster(gcSample, fun="enrichKEGG",
                      organism="hsa", pvalueCutoff=0.05)


xx3@compareClusterResult = toPlot4

p3 = dotplot(xx3, showCategory=7, size = "count")


ggexport(filename = "./Dados/Processados/V1/V1_GO_masters/kegg_maior.png",
         plot = p3 , device = "png",
         width = 1700,
         height = 2700, res = 330)



# Load network, masters


load("./Dados/Processados/V1/V1_MRA_SHH/MRA_SHH_7.rdata")
SHH = master
load("./Dados/Processados/V1/V1_MRA_G3/MRA_G3_7.rdata")
G3 = master
load("./Dados/Processados/V1/V1_MRA_G4/MRA_G4_7.rdata")
G4 = master

ALL = rbind(SHH, G3, G4)
master = ALL[!duplicated(ALL$Regulon),]

# Create graph

g <- tni.graph(rtni_final_sub, regulatoryElements = master$Regulon, tnet="ref", gtype="amapDend")

lt = list(SHH = SHH$Regulon, G3 = G3$Regulon, G4 = G4$Regulon)
mt = list_to_matrix(lt)

# Setting colors

unicos = mt[which(rowSums(mt) == 1),]
dual = mt[which(rowSums(mt) == 2),]
tres = mt[which(rowSums(mt) == 3),]


regs_SHH = names(which(unicos[,1] == 1))
regs_G3 = names(which(unicos[,2] == 1))
regs_G4 = names(which(unicos[,3] == 1))
regs_SHH_G3 = names(which(dual[,3] == 0))
regs_SHH_G4 = names(which(dual[,2] == 0))
regs_G3_G4 = names(which(dual[,1] == 0))
regs_SHH_G3_G4 = rownames(tres)


# Genes SHH

regulons <- tni.get(rtni_final_sub, what = "regulons.and.mode", idkey = "ID")
regulons = regulons[which(names(regulons) %in% regs_SHH)] 


i = 1
genes = NULL

while (i <= length(regulons)) {
  
  genes = c(genes, names(regulons[[i]]))
  i = i + 1
  
}

genes = genes[!duplicated(genes)]

SHH = mapIds(org.Hs.eg.db, alias2Symbol(genes),"ENTREZID", "SYMBOL")

## Genes G3

regulons <- tni.get(rtni_final_sub, what = "regulons.and.mode", idkey = "ID")
regulons = regulons[which(names(regulons) %in% regs_G3)] 


i = 1
genes = NULL

while (i <= length(regulons)) {
  
  genes = c(genes, names(regulons[[i]]))
  i = i + 1
  
}

genes = genes[!duplicated(genes)]

G3 = mapIds(org.Hs.eg.db, alias2Symbol(genes),"ENTREZID", "SYMBOL")

## Genes G4

regulons <- tni.get(rtni_final_sub, what = "regulons.and.mode", idkey = "ID")
regulons = regulons[which(names(regulons) %in% regs_G4)] 


i = 1
genes = NULL

while (i <= length(regulons)) {
  
  genes = c(genes, names(regulons[[i]]))
  i = i + 1
  
}

genes = genes[!duplicated(genes)]

G4 = mapIds(org.Hs.eg.db, alias2Symbol(genes),"ENTREZID", "SYMBOL")


## Genes SHH/G3/G4

regulons <- tni.get(rtni_final_sub, what = "regulons.and.mode", idkey = "ID")
regulons = regulons[which(names(regulons) %in% regs_SHH_G3_G4)] 


i = 1
genes = NULL

while (i <= length(regulons)) {
  
  genes = c(genes, names(regulons[[i]]))
  i = i + 1
  
}

genes = genes[!duplicated(genes)]

ALL_3 = mapIds(org.Hs.eg.db, alias2Symbol(genes),"ENTREZID", "SYMBOL")



clusters1 <- list(SHH, G3,G4, ALL_3)

names(clusters1) <- c('SHH', 'G3',"G4", "SHH/G3/G4")


CompareGO_BP = compareCluster(clusters1, fun="enrichGO", pvalueCutoff=0.01,
                              pAdjustMethod="BH", OrgDb=org.Hs.eg.db, ont="BP",
                              readable=T)

dotplot(CompareGO_BP, showCategory=3, title="GO - BP 7 MRs Network 1")



## Not run: 
data(gcSample)
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)


xx@compareClusterResult = CompareGO_BP@compareClusterResult
dotplot(xx, showCategory=3, title="GO - BP 7 MRs Network 1")







