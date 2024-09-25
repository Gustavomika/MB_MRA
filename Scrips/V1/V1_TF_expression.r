library(oligo)
library(data.table)
library(arrayQualityMetrics)
library(dplyr)
library(limma)
library(ggplot2)
library(ggstatsplot)
library(ggpubr)


# Different expression analyses of control (fetal and adult cerebellum) vs MB samples

# set to your working directory
setwd("C:/Users/guga_/Documents/tayrone")



###################


# Loads cell files, outlier samples and signature samples

celFiles<- list.files("./Dados/Brutos/Amostras/GSE167447", full.names = T)
celFiles = celFiles[-1]


# Get only fetal samples


celFiles1 <- list.files("./Dados/Brutos/Amostras/GSE85217", full.names = T)
celFiles1 = celFiles1[-1]

load("./Dados/Processados/V1/V1_express_norm/V1_outlier.rdata")


tabelao = read.csv("./Dados/Brutos/MetaDados/tabelao.txt", sep = "\t")
tabelao = tabelao[-outlier,]

# Remove outliers and signature samples

celFiles1 = celFiles1[-outlier]


celFiles1 = celFiles1[-which(tabelao$Subgroup == "WNT")]
tabelao = tabelao[-which(tabelao$Subgroup == "WNT"),]



celFiles2 = c(celFiles, celFiles1)

# Checks the array quality

gexp <- read.celfiles(celFiles2)
gexp <- rma(gexp)



gexp <- exprs(gexp)
gexp <- as.data.frame(gexp)


# Reads the GPL dictionary 

platform_data <- as.data.frame(data.table::fread(file = "./Dados/Brutos/AnotacoesPlataforma/GPL11532-32230.txt",
                                                 header = T, skip = 12))

# Gets probe IDs for construction of dictionary

dictionary <- as.data.frame(platform_data$ID, stringsAsFactors = FALSE) 

# Gene_asigment column will be used to retrieve gene symbols

gene_symbol <- platform_data$gene_assignment 

# Splits string in order to obtain only gene symbol

gene_symbol <- strsplit(gene_symbol, " // ") 

# Collects second value of each element, which is HGNC gene symbol

gene_symbol <- as.data.frame(sapply(gene_symbol, function(x) x[2])) 

# Creates dictionary with probe ids and gene symbols

dictionary <- cbind(dictionary, gene_symbol)
colnames(dictionary) <- c("probe_id", "gene_symbol")

# Removes rows that gene symbol is not available, probably because the probe does not belong to a gene

dictionary <- as.data.frame(dictionary[!is.na(dictionary$gene_symbol), ])

rownames(dictionary) <- dictionary$probe_id

# Removes temporary objects

rm(gene_symbol, platform_data) 


dictionary = dictionary[-which(is.na(alias2SymbolTable(dictionary$gene_symbol)) == T),]
dictionary$gene_symbol = alias2SymbolTable(dictionary$gene_symbol)



# Subset lines from expression according to dictionary


gexp = gexp[rownames(gexp) %in% dictionary[,1], ]
dictionary = dictionary[dictionary[,1] %in% rownames(gexp),]

# Orders the data

gexp = gexp[order(rownames(gexp)),]
dictionary <- dictionary[order(dictionary[,1]), ]
rownames(dictionary) <- dictionary$probe_id


# Must be true

all(dictionary[,1] == rownames(gexp))

# There`re genes with more than one probe
# The goal of this section is to select the probes with the highest variance of each gene
# And remove of others from both dictionary and express data


# Gets the duplicated genes

uniqueSymbols  = as.character(unique(dictionary[which(duplicated(dictionary[,2])),2]))

# Returns a vector of probe IDs with greatest variance 

get_son = function() {
  
  gg =  cbind(gexp, dictionary[,2], rownames(gexp))
  k = 1
  sondas_corretas = vector(mode = "character", length = length(uniqueSymbols))
  
  
  while (k <= length(uniqueSymbols)) {
    
    
    buffer = gg[gg$`dictionary[, 2]` == uniqueSymbols[k],]
    
    l = length(buffer[,1])
    
    maxvar = 0
    i = 1
    
    while (i <= l)
    {
      variancia = var(as.numeric(buffer[i,1:length(gexp)]))
      if (variancia > maxvar)
      {
        maxvar = variancia
        sondas_corretas[k] = rownames(buffer[i,])
      }
      i = i + 1
    }
    
    k = k +1
  }
  
  return(sondas_corretas)
}

sondas = get_son()


# Creates signature_dictionary as a subset of dictionary

repeated <- which(dictionary$gene %in% uniqueSymbols)
signature_dictionary <- dictionary[-repeated,]
signature_dictionary = rbind(signature_dictionary, dictionary[which(dictionary$probe_id %in% sondas),])

# Creates signature_gexp according to signature_dictionary

signature_gexp = gexp[rownames(gexp) %in% signature_dictionary$probe_id,]


# Order

signature_gexp = signature_gexp[order(rownames(signature_gexp)),]
signature_dictionary = signature_dictionary[order(rownames(signature_dictionary)),]
rownames(signature_gexp) = signature_dictionary$gene_symbol


# Main master regulators identified (risk MRs and shared)


genes = c(which(rownames(signature_gexp) == "BHLHE41"),
          which(rownames(signature_gexp) == "RFX4"),
          which(rownames(signature_gexp) == "NPAS3"),
          which(rownames(signature_gexp) == "REL"),
          which(rownames(signature_gexp) == "ZFAT"),
          which(rownames(signature_gexp) == "MYC"),
          which(rownames(signature_gexp) == "ZSCAN5A"),
          which(rownames(signature_gexp) == "PAX6"),
          which(rownames(signature_gexp) == "ZNF157"),
          which(rownames(signature_gexp) == "ARNT2"),
          which(rownames(signature_gexp) == "HIVEP3"))


genes = signature_gexp[genes,]
genes = t(genes)


genes = as.data.frame(genes)



genes[,12] = c(rep("Fetal Control",8), rep("Adult Control", 5), tabelao$Subgroup)

genes$V12[which(genes$V12 == "Group3")] = "Group 3"
genes$V12[which(genes$V12 == "Group4")] = "Group 4"

genes$V12 = as.factor(genes$V12)


genes$V12 <- ordered(genes$V12, levels = c("Fetal Control","Adult Control","SHH","Group 3","Group 4"))

 
i = 1

k2 = ggbetweenstats(
    data = genes,
    x = V12,
    y = !!colnames(genes)[i],
    type = "nonparametric",
    p.adjust.method = "bonferroni",
    title = paste(colnames(genes)[i], "Groups Comparision"),
    xlab = "",
    ylab = paste(colnames(genes)[i], "Log Expression"),
    pairwise.comparisons = T,
    results.subtitle = F,
    bf.message = F,
    pairwise.display = "none"
    
    )+ 
    theme_bw()+
    theme(legend.position="none")+
  geom_signif(annotations = '***', y_position = 12.4, xmin = 3, xmax = 5, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '**', y_position = 12.2, xmin = 3, xmax = 4, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '***', y_position = 12, xmin = 1, xmax = 3, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '**', y_position = 11.8, xmin = 1, xmax = 5, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '**', y_position = 11.6, xmin = 1, xmax = 4, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '***', y_position = 11.4, xmin = 2, xmax = 3, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '**', y_position = 11.2, xmin = 2, xmax = 5, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '**', y_position = 11, xmin = 2, xmax = 4, textsize = 4.5, vjust = 0.9)+ 
  scale_color_manual(values=c("#998789", "#998729", "red3","yellow3","green4"))

  


ggexport(filename = paste("./Dados/Processados/V1/V1_TF_expression/",colnames(genes)[i], "_adult.png", sep = ""),
         plot = k2 ,device = "png",
         width = 1700*1.0,
         height = 2100*1.0, res = 390)



i = 2

k2 = ggbetweenstats(
  data = genes,
  x = V12,
  y = !!colnames(genes)[i],
  type = "nonparametric",
  p.adjust.method = "bonferroni",
  title = paste(colnames(genes)[i], "Groups Comparision"),
  xlab = "",
  ylab = paste(colnames(genes)[i], "Log Expression"),
  pairwise.comparisons = T,
  results.subtitle = F,
  bf.message = F,
  pairwise.display = "none"
  
)+ 
  theme_bw()+
  theme(legend.position="none")+
  geom_signif(annotations = '***', y_position = 13.1, xmin = 3, xmax = 5, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '***', y_position = 12.8, xmin = 1, xmax = 3, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '***', y_position = 12.5, xmin = 1, xmax = 5, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '***', y_position = 12.2, xmin = 1, xmax = 4, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '**', y_position = 11.9, xmin = 2, xmax = 3, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '*', y_position = 11.6, xmin = 2, xmax = 5, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '**', y_position = 11.3, xmin = 2, xmax = 4, textsize = 4.5, vjust = 0.9)+ 
  scale_color_manual(values=c("#998789", "#998729", "red3","yellow3","green4"))
  
  
  
  
  ggexport(filename = paste("./Dados/Processados/V1/V1_TF_expression/",colnames(genes)[i], "_adult.png", sep = ""),
           plot = k2 ,device = "png",
           width = 1700*1.0,
           height = 2100*1.0, res = 390)




i = 3

k2 = ggbetweenstats(
  data = genes,
  x = V12,
  y = !!colnames(genes)[i],
  type = "nonparametric",
  p.adjust.method = "bonferroni",
  title = paste(colnames(genes)[i], "Groups Comparision"),
  xlab = "",
  ylab = paste(colnames(genes)[i], "Log Expression"),
  pairwise.comparisons = T,
  results.subtitle = F,
  bf.message = F,
  pairwise.display = "none"
  
)+ 
  theme_bw()+
  theme(legend.position="none")+
geom_signif(annotations = '***', y_position = 10.4, xmin = 1, xmax = 3, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '***', y_position = 10.2, xmin = 1, xmax = 5, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '***', y_position = 10, xmin = 1, xmax = 4, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '**', y_position = 9.8, xmin = 2, xmax = 3, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '**', y_position = 9.6, xmin = 2, xmax = 5, textsize = 4.5, vjust = 0.9)+
  geom_signif(annotations = '**', y_position = 9.4, xmin = 2, xmax = 4, textsize = 4.5, vjust = 0.9)+ 
  scale_color_manual(values=c("#998789", "#998729", "red3","yellow3","green4"))
  
  
  
  ggexport(filename = paste("./Dados/Processados/V1/V1_TF_expression/",colnames(genes)[i], "_adult.png", sep = ""),
           plot = k2 ,device = "png",
           width = 1700*1.0,
           height = 2100*1.0, res = 390)




save(signature_gexp, file = "./Dados/Processados/V1/V1_TF_expression/signature_gexp.rdata")

