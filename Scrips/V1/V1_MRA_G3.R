library(RTN)

# Master regulation analysis of G3

# set to your working directory
setwd("C:/Users/guga_/Documents/tayrone")

# Load network, DEGs and logFC

load("./Dados/Processados/V1/V1_rede_sub_7/rtni_final_sub_7.rdata")
load("./Dados/Processados/V1/V1_GSE167447_GSE85217_G3/DEGs.rdata")
load("./Dados/Processados/V1/V1_GSE167447_GSE85217_G3/ALL.rdata")

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

# Save master from G3
 
save(mra, master, file = "./Dados/Processados/V1/V1_MRA_G3/MRA_G3_7.rdata")

