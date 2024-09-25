
library(RTN)
library(limma)
library(snow)

# Script for MB network construction

# set to your working directory
setwd("C:/Users/guga_/Documents/tayrone")

# Load normalized case samples and human TFs list

load("./Dados/Processados/V1/V1_express_norm_rede/V1_signature_gexp.rdata")
TFs = read.table("./Dados/Brutos/TFs/TFs.txt")

# Filter TFs and gexp to make sure all genes in one
# object are in the other
# alias2SymbolTable function upddates gene symbols from both data

filtro = which(is.na(alias2SymbolTable(rownames(signature_gexp), species = "Hs")))
filtro2 = which(duplicated(alias2SymbolTable(rownames(signature_gexp), species = "Hs"), incomparables = NA))
signature_gexp = signature_gexp[-c(filtro, filtro2),]

rownames(signature_gexp) = alias2SymbolTable(rownames(signature_gexp), species = "Hs")


filtroTFs = which(is.na(alias2SymbolTable(TFs$V1, species = "Hs")))
TFs = as.data.frame(TFs[-filtroTFs,])


######

# network construction from RTN

rtni_sub <- tni.constructor(expData = as.matrix(signature_gexp), 
                        regulatoryElements = TFs$`TFs[-filtroTFs, ]`)

options(cluster=snow::makeCluster(3, "SOCK"))

rtni1_sub <- tni.permutation(rtni_sub, nPermutations = 1000, pValueCutoff = 1e-7)
rtni2_sub <- tni.bootstrap(rtni1_sub)

stopCluster(getOption("cluster"))


rtni_final_sub <- tni.dpi.filter(rtni2_sub)

# save

save(rtni_sub, file = "./Dados/Processados/V1/V1_rede_sub_7/rtni_sub_7.rdata")
save(rtni1_sub, rtni2_sub, rtni_final_sub, file = "./Dados/Processados/V1/V1_rede_sub_7/rtni_sub_tudo_7.rdata")
save(rtni_final_sub, file = "./Dados/Processados/V1/V1_rede_sub_7/rtni_final_sub_7.rdata")





