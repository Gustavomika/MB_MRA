library(RTN)

# Master regulation analysis of SHH

# set to your working directory
setwd("C:/Users/guga_/Documents/tayrone")

# Load network, DEGs and logFC

load("./Dados/Processados/V1/V1_rede_sub_7/rtni_final_sub_7.rdata")
load("./Dados/Processados/V1/V1_GSE167447_GSE85217_SHH/DEGs.rdata")
load("./Dados/Processados/V1/V1_GSE167447_GSE85217_SHH/ALL.rdata")

# Retrieve logFCs

logFC = ALL$logFC
names(logFC) = rownames(ALL)

# Master regulation analysis from RTN

rtna <- tni2tna.preprocess(object = rtni_final_sub, 
                           phenotype = logFC, 
                           hits = rownames(DEGs))

rtna <- tna.mra(rtna, pAdjustMethod = "bonferroni")

mra <- tna.get(rtna, what="mra", ntop = -1)
head(mra)

master = mra[which(mra$Adjusted.Pvalue <= 0.01),]

# Save master from SHH

save(mra, master, file = "./Dados/Processados/V1/V1_MRA_SHH/MRA_SHH_7.rdata")

# Just a plot

g <- tni.graph(rtni_final_sub, regulatoryElements = master$Regulon, tnet="ref", gtype="amapDend")

library(RedeR)
rdp <- RedPort()
calld(rdp)
addGraph(rdp, g$g, layout=NULL)
addLegend.color(rdp, g, type="edge")
addLegend.shape(rdp, g)
relax(rdp, ps = TRUE)
