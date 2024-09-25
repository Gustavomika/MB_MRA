library(RTN)
library(RTNsurvival)
library(pheatmap)
library(openxlsx)
library(snow)
library(survival)
library(maxstat)
library(survminer)
library(RedeR)
library(broom)

# This script identify and plot the risk master regulators using cox regression and kaplan-meier curves

# set to your working directory
setwd("C:/Users/guga_/Documents/tayrone")


load("./Dados/Processados/V1/V1_rede_sub_7/rtni_final_sub_7.rdata")
load("./Dados/Processados/V1/V1_MRA_SHH/MRA_SHH_7.rdata")
SHH = master
load("./Dados/Processados/V1/V1_MRA_G3/MRA_G3_7.rdata")
G3 = master
load("./Dados/Processados/V1/V1_MRA_G4/MRA_G4_7.rdata")
G4 = master

j = c(SHH$Regulon, G3$Regulon, G4$Regulon)
j = unique(j)


load("./Dados/Processados/V1/V1_rede_sub_7/rtni_final_sub_7.rdata")
design = read.xlsx("./Dados/Brutos/MetaDados/1-s2.0-S1535610817302015-mmc2.xlsx", startRow = 2)
tabelao = read.csv("./Dados/Brutos/MetaDados/tabelao.txt", sep = "\t")

colnames(design)[5] = "Histology"

back = rtni_final_sub@colAnnotation[["ID"]]

rtni_final_sub@colAnnotation[["ID"]] = substr(rtni_final_sub@colAnnotation[["ID"]], 1, 10)

design[,2] = tabelao$Accession
design = design[which(design$Age %in% rtni_final_sub@colAnnotation[["ID"]]),]


rtni_final_sub@colAnnotation[["ID"]] = back
rownames(design) = back

design$Histology = as.factor(design$Histology)
design$Subtype = as.factor(design$Subtype)
colnames(design)[6] = "Met"

design = design[which(design$Subgroup != "WNT"),]

rtns <- tni2tnsPreprocess(rtni_final_sub, survivalData = design, time = 8, event = 7, 
                          regulatoryElements = j,
                          samples = rownames(design))

options(cluster=snow::makeCluster(7, "SOCK"))
rtns <- tnsGSEA2(rtns)
stopCluster(getOption("cluster"))

rtns1 = rtns

rtns <- tnsCox(rtns)




ann_colors = list(
  Met = c("blue3", "red3"),
  Histology = c(Classic = "#1B9E37", Desmoplastic = "#D95F02",LCA = "cyan",MBEN = "#E7298A"),
  Subtype = c(Group3_alpha = "#7570B3",Group3_beta =  "#E7298A", Group3_gamma = "#66A61E",
              Group4_alpha = "darkgreen",Group4_beta = "coral1",Group4_gamma = "blue",
              SHH_alpha = "blueviolet",SHH_beta = "brown1",SHH_gamma = "darkslategrey",SHH_delta = "cyan",
              WNT_alpha = "darkgreen",WNT_alpha = "darkgrey"),
  Subgroup = c(Group3 = "yellow2",Group4 =  "green3", SHH = "red3")
)

survival.data <- tnsGet(rtns1, "survivalData")
annotation_col <- survival.data[,c("Subgroup", "Subtype", "Histology", "Met")]
enrichmentScores <- tnsGet(rtns1, "regulonActivity")
p = pheatmap(t(enrichmentScores$dif),
             show_colnames = FALSE, 
             annotation_colors = ann_colors,
             annotation_col  = annotation_col,
             annotation_legend = FALSE,
             main = "Medulloblastoma Master Regulons")

save_pheatmap <- function(x, filename, width=480, height=960) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename,width = width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


save_pheatmap(p, filename = "./Dados/Processados/V1/V1_survival_ALL/ALL_masters.png",
              width = 1500)




#############

kk =rtns@results[["Cox"]][["Table"]]

rtns_teste = tni2tnsPreprocess(rtni_final_sub, survivalData = design, time = 8, event = 7, samples = rownames(design),
                               keycovar = c("Subtype","Histology", "Met"),
                               regulatoryElements = kk$Regulons[kk$Adjusted.Pvalue < 0.01])


rtns_teste <- tnsGSEA2(rtns_teste)
rtns_teste <- tnsCox(rtns_teste)


enrichmentScores <- tnsGet(rtns_teste, "regulonActivity")
survival.data <- tnsGet(rtns_teste, "survivalData")
annotation_col <- survival.data[,c("Subgroup", "Histology")]
pp = pheatmap(t(enrichmentScores$dif), 
              show_colnames = FALSE,
              annotation_col  = annotation_col,
              annotation_colors = ann_colors,
              main = "Medulloblastoma Risk Master Regulons",
              annotation_legend = T, 
              clustering_method = "ward.D",
              fontsize = 12)

save_pheatmap(pp, filename = "./Dados/Processados/V1/V1_survival_ALL/ALL_masters_risk3.png",
              width = 1000,height = 333)




g <- tni.graph(rtni_final_sub, regulatoryElements = kk$Regulons[kk$Adjusted.Pvalue < 0.01],
               tnet="ref", gtype="amap", amapFilter = "phyper")


rdp <- RedPort()
calld(rdp)
addGraph(rdp, g, layout=NULL)



dado = t(enrichmentScores$dif)


dado = dado[pp[["tree_row"]][["order"]],pp[["tree_col"]][["order"]]]
survival.data = survival.data[pp[["tree_col"]][["order"]],]


dado = t(dado)
survival.data = survival.data[,1:2]
dado = cbind(survival.data,dado)


dado = dado[-which(is.na(dado$time)),]
dado = dado[-which(is.na(dado$event)),]


# Kaplan Meier

i = 3

while (i <= length(dado[1,])) {
  

res.cat = dado[,c(1,2,i)]
colnames(res.cat)[3] = "gene"

res.cat[which(res.cat$gene > 0),3] = "high"
res.cat[which(res.cat$gene < 0),3] = "low"

km_gene_opt = survfit(Surv(time, event) ~ gene , data = res.cat)


# Obtenção do raw p 

raw_p = survdiff(Surv(time, event) ~ gene, data = res.cat)
raw_p = raw_p$pvalue

# Plot
k2 =  ggsurvplot(km_gene_opt, conf.int=TRUE, pval.method = T,pval=paste("
                                                                        
                                                                        
                                                                        
p  = ",format(round(raw_p, 5), nsmall = 5)), risk.table=TRUE,
legend.labs=c("High", "Low"), legend.title="Activity",  
palette=c("orchid2","dodgerblue2"), 
title=paste(colnames(dado)[i], "Regulon Activity"), 
risk.table.height=.19)


ggexport(filename = paste("./Dados/Processados/V1/V1_survival_ALL/",colnames(dado)[i], "_kaplan.png", sep = ""),
         plot = k2 , device = "png",
         width = 1700*0.81,
         height = 2700*0.81, res = 330)


i = i + 1

}




res.cat$gene[1:122] = "high"
res.cat$gene[123:528] = "low"


km_gene_opt = survfit(Surv(time, event) ~ gene , data = res.cat)


# Obtenção do raw p 

raw_p = survdiff(Surv(time, event) ~ gene, data = res.cat)
raw_p = raw_p$pvalue

# Plot
k2 =  ggsurvplot(km_gene_opt, conf.int=TRUE, pval.method = T,pval=paste("

                                                                        
                                                                        
                                                                        
                                                                        
p-value (Log-Rank)  = ",format(round(raw_p, 5), nsmall = 5)), risk.table=TRUE,
legend.labs=c("High", "Low"), legend.title="Expression",  
palette=c("orchid2","dodgerblue2"), 
title=paste("Medulloblastoma Risk Regulons Activity"), 
risk.table.height=.15) 


ggexport(filename = paste("./Dados/Processados/V1/V1_survival_ALL/ALL_kaplan.png", sep = ""),
         plot = k2 , device = "png",
         width = 1700,
         height = 2700, res = 330)




