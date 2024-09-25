
# Script to select 8 samples of each MB subgroup, that will be the signature
# samples considering its subtypes incidence  

# set to your working directory
setwd("C:/Users/guga_/Documents/tayrone")

# Load meta data and outlier samples 

tabelao = read.csv("./Dados/Brutos/MetaDados/tabelao.txt", sep = "\t")
load("./Dados/Processados/V1/V1_express_norm/V1_outlier.rdata")

# Remove outlier samples 

tabelao = tabelao[-outlier,]
rownames(tabelao) = 1:756


# Calculate the number of each subtype in each subgroup

SHH_alpha = round(length(which(tabelao$Subtype == "SHH_alpha"))*8/length(which(tabelao$Subgroup == "SHH")))
SHH_beta = round(length(which(tabelao$Subtype == "SHH_beta"))*8/length(which(tabelao$Subgroup == "SHH")))
SHH_gamma = round(length(which(tabelao$Subtype == "SHH_gamma"))*8/length(which(tabelao$Subgroup == "SHH")))
SHH_delta = round(length(which(tabelao$Subtype == "SHH_delta"))*8/length(which(tabelao$Subgroup == "SHH")))

Group3_alpha = round(length(which(tabelao$Subtype == "Group3_alpha"))*8/length(which(tabelao$Subgroup == "Group3")))
Group3_beta = round(length(which(tabelao$Subtype == "Group3_beta"))*8/length(which(tabelao$Subgroup == "Group3")))
Group3_gamma = round(length(which(tabelao$Subtype == "Group3_gamma"))*8/length(which(tabelao$Subgroup == "Group3")))

Group4_alpha = round(length(which(tabelao$Subtype == "Group4_alpha"))*8/length(which(tabelao$Subgroup == "Group4")))
Group4_beta = round(length(which(tabelao$Subtype == "Group4_beta"))*8/length(which(tabelao$Subgroup == "Group4")))
Group4_gamma = round(length(which(tabelao$Subtype == "Group4_gamma"))*8/length(which(tabelao$Subgroup == "Group4")))


# Get samples index

set.seed(1)
SHH_alpha_samples = which(tabelao$Subtype == "SHH_alpha")[sample(length(which(tabelao$Subtype == "SHH_alpha")),SHH_alpha) ]

set.seed(1)
SHH_beta_samples = which(tabelao$Subtype == "SHH_beta")[sample(length(which(tabelao$Subtype == "SHH_beta")),SHH_beta) ]

set.seed(1)
SHH_gamma_samples = which(tabelao$Subtype == "SHH_gamma")[sample(length(which(tabelao$Subtype == "SHH_gamma")),SHH_gamma) ]

set.seed(1)
SHH_delta_samples = which(tabelao$Subtype == "SHH_delta")[sample(length(which(tabelao$Subtype == "SHH_delta")),SHH_delta) ]



set.seed(1)
Group3_alpha_samples = which(tabelao$Subtype == "Group3_alpha")[sample(length(which(tabelao$Subtype == "Group3_alpha")),Group3_alpha) ]

set.seed(1)
Group3_beta_samples = which(tabelao$Subtype == "Group3_beta")[sample(length(which(tabelao$Subtype == "Group3_beta")),Group3_beta) ]

set.seed(1)
Group3_gamma_samples = which(tabelao$Subtype == "Group3_gamma")[sample(length(which(tabelao$Subtype == "Group3_gamma")),Group3_gamma) ]



set.seed(1)
Group4_alpha_samples = which(tabelao$Subtype == "Group4_alpha")[sample(length(which(tabelao$Subtype == "Group4_alpha")),Group4_alpha) ]

set.seed(1)
Group4_beta_samples = which(tabelao$Subtype == "Group4_beta")[sample(length(which(tabelao$Subtype == "Group4_beta")),Group4_beta) ]

set.seed(1)
Group4_gamma_samples = which(tabelao$Subtype == "Group4_gamma")[sample(length(which(tabelao$Subtype == "Group4_gamma")),Group4_gamma) ]

# Create data frame of signature samples meta data

tabelao_filt = tabelao[c(SHH_alpha_samples,
                    SHH_beta_samples,
                    SHH_gamma_samples,
                    SHH_delta_samples,
                    Group3_alpha_samples,
                    Group3_beta_samples,
                    Group3_gamma_samples,
                    Group4_alpha_samples,
                    Group4_beta_samples,
                    Group4_gamma_samples),]


# save

save(tabelao_filt, file = "./Dados/Processados/V1/V1_signature_samples/tabelao_filt.rdata")





