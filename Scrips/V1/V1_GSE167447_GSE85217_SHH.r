library(oligo)
library(data.table)
library(arrayQualityMetrics)
library(dplyr)
library(limma)
library(ggplot2)


# Different expression analyses of control vs SHH samples

# set to your working directory
setwd("C:/Users/guga_/Documents/tayrone")


# Loads cell files (control)

celFiles<- list.files("./Dados/Brutos/Amostras/GSE167447", full.names = T)

# Get only fetal samples

celFiles = celFiles[1:8]


###################


# Loads cell files, outlier samples and signature samples

celFiles1 <- list.files("./Dados/Brutos/Amostras/GSE85217", full.names = T)
load("./Dados/Processados/V1/V1_express_norm/V1_outlier.rdata")
load("./Dados/Processados/V1/V1_signature_samples/tabelao_filt.rdata")

# Remove outliers and signature samples

celFiles1 = celFiles1[-outlier]
celFiles1 = celFiles1[as.numeric(rownames(tabelao_filt[which(tabelao_filt$Subgroup == "SHH"),]))]



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



#################################



## Different expression analyses from limma


labels = c(rep("Control", 8), rep("Case", 8))


design <- model.matrix(~0+labels)
colnames(design) <- c("Case", "Control")

# A linear model (MArrayLM) is fitted to the expression data for each probe

fit <- lmFit(signature_gexp, design)

# Defines case-control as contrasts variable, which will be used on next steps

contrasts <- makeContrasts(Case-Control, levels=design) 

# Applies function to rank genes in order of evidence for differential expression

ct.fit = contrasts.fit(fit, contrasts)
ct.fit <- eBayes(ct.fit)


# Checks whether each statistic should be considered significantly different from zero


DEGs <- topTable(ct.fit, adjust.method = "BH", p.value=  0.05, number = 9999999, lfc = 1.5)
ALL <- topTable(ct.fit, adjust.method = "BH", p.value=  1, number = 9999999, lfc = 0)

# Save

save(DEGs, file = "./Dados/Processados/V1/V1_GSE167447_GSE85217_SHH/DEGs.rdata")
save(ALL, file = "./Dados/Processados/V1/V1_GSE167447_GSE85217_SHH/ALL.rdata")


