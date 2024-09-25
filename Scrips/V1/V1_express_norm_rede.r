library(oligo)
library(data.table)
library(arrayQualityMetrics)
library(dplyr)
library(limma)
library(ggplot2)

# Renormalization of case samples without outliers and signature samples

# set to your working directory
setwd("C:/Users/guga_/Documents/tayrone")

# Loads cell files, outlier samples and signature samples

celFiles<- list.files("./Dados/Brutos/Amostras/GSE85217", full.names = T)
load("./Dados/Processados/V1/V1_express_norm/V1_outlier.rdata")
load("./Dados/Processados/V1/V1_signature_samples/tabelao_filt.rdata")

celFiles = celFiles[-outlier]
celFiles = celFiles[-as.numeric(rownames(tabelao_filt))]


gexp <- read.celfiles(celFiles)

# Robust Multiarray Average for data normalization

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

# Save

save(signature_gexp, file = "./Dados/Processados/V1/V1_express_norm_rede/V1_signature_gexp.rdata")


