library(RTN)
library(RColorBrewer)
library(ComplexHeatmap)
library(igraph)
library(RedeR)


# This script create and plot graphs through RedeR package. RedeR uses a Java aplication for network visualization,
# so it have to be done manually in the aplication

# set to your working directory
setwd("C:/Users/guga_/Documents/tayrone")

# Load network, masters

load("./Dados/Processados/V1/V1_rede_sub_7/rtni_final_sub_7.rdata")
load("./Dados/Processados/V1/V1_MRA_SHH/MRA_SHH_7.rdata")

# Create graph

g <- tni.graph(rtni_final_sub, regulatoryElements = master$Regulon, tnet="ref", gtype="amapDend")

# Set node colors

x = -log(master$Adjusted.Pvalue[match(igraph::V(g$g)$name, master$Regulon)])

rbPal <- colorRampPalette(c('#F2E5E5','firebrick4'))
colors <- rbPal(10)[as.numeric(cut(x,breaks = 10))]

# Set nnode size

y = master$Regulon.Size[match(igraph::V(g$g)$name, master$Regulon)]
y[which(is.na(y))] = 10


igraph::V(g$g)$nodeFontSize[which(igraph::V(g$g)$nodeFontSize == 20)] = 30
igraph::V(g$g)$nodeColor = colors
igraph::V(g$g)$nodeSize = y

# Call net

rdp <- RedPort()
calld(rdp)
addGraph(rdp, g$g, layout=NULL)

scl <- rbPal(10)
leg <- c("5","15","25","35","45","55","65","75","85","95")
addLegend.color(rdp, colvec=scl, labvec=leg, title="-Log(Adjusted P-value)", size = 40,
                ftsize = 15)

addLegend.size(obj = rdp, g$g, title = "Regulon Size",
               position = "bottomright", sizevec = as.numeric(c(20,50,100,200,300)),
               intersp = 7)


#### Second Plot

# Load network, masters

load("./Dados/Processados/V1/V1_rede_sub_7/rtni_final_sub_7.rdata")

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



igraph::V(g$g)$nodeColor[!is.na(match(igraph::V(g$g)$name, regs_SHH))] = "red"
igraph::V(g$g)$nodeColor[!is.na(match(igraph::V(g$g)$name, regs_G3))] = "yellow"
igraph::V(g$g)$nodeColor[!is.na(match(igraph::V(g$g)$name, regs_G4))] = "green"
igraph::V(g$g)$nodeColor[!is.na(match(igraph::V(g$g)$name, regs_SHH_G3))] = "orange3"
igraph::V(g$g)$nodeColor[!is.na(match(igraph::V(g$g)$name, regs_SHH_G4))] = "pink3"
igraph::V(g$g)$nodeColor[!is.na(match(igraph::V(g$g)$name, regs_G3_G4))] = "yellowgreen"
igraph::V(g$g)$nodeColor[!is.na(match(igraph::V(g$g)$name, regs_SHH_G3_G4))] = "gold4"

# Set size

y = master$Regulon.Size[match(igraph::V(g$g)$name, master$Regulon)]
y[which(is.na(y))] = 10


igraph::V(g$g)$nodeFontSize[which(igraph::V(g$g)$nodeFontSize == 20)] = 38
igraph::V(g$g)$nodeSize = y

# Call net


rdp <- RedPort()
calld(rdp)
addGraph(rdp, g$g, layout=NULL)

scl <- c("red","yellow","green", "gold4","pink3","orange3","yellowgreen")
leg <- c("SHH","G3","G4","SHH/G3/G4","SHH/G4","SHH/G3","G3/G4")
addLegend.color(rdp, colvec=scl, labvec=leg, title="Groups", size = 60,
                ftsize = 15, position = "topleft")

addLegend.size(obj = rdp, g$g, title = "Regulon Size",
               position = "bottomright", sizevec = as.numeric(c(20,50,100,200,300)),
               intersp = 7)


g2 <- tni.graph(rtni_final_sub, tnet="ref", gtype="amapDend")


graf = g2$g

save(graf, file = "./Dados/Processados/V1/V1_rede_plots/graf.rdata")
load("./Dados/Processados/V1/V1_rede_plots/graf.rdata")

coords <- layout_as_tree(graf, circular = T, flip.y = F)

igraph::V(graf)$nodeSize = igraph::V(graf)$nodeSize/10
E(graf)$edgeWidth <- 1

# Regulon Font
igraph::V(graf)$nodeFontSize[which(igraph::V(graf)$nodeFontSize == 20)] = 0
# Fonte nos internos
igraph::V(graf)$nodeFontSize[which(igraph::V(graf)$nodeFontSize == 1)] = 0

# Fonte mestres
igraph::V(graf)$nodeFontSize[!is.na(match(igraph::V(graf)$name, master$Regulon))] = 5

# Cor fonte dos mestres

igraph::V(graf)$nodeFontColor = "black"

igraph::V(graf)$nodeFontColor[!is.na(match(igraph::V(graf)$name, regs_SHH))] = "red4"
igraph::V(graf)$nodeFontColor[!is.na(match(igraph::V(graf)$name, regs_G3))] = "yellow4"
igraph::V(graf)$nodeFontColor[!is.na(match(igraph::V(graf)$name, regs_G4))] = "green4"
igraph::V(graf)$nodeFontColor[!is.na(match(igraph::V(graf)$name, regs_SHH_G3))] = "orange3"
igraph::V(graf)$nodeFontColor[!is.na(match(igraph::V(graf)$name, regs_SHH_G4))] = "pink3"
igraph::V(graf)$nodeFontColor[!is.na(match(igraph::V(graf)$name, regs_G3_G4))] = "yellowgreen"
igraph::V(graf)$nodeFontColor[!is.na(match(igraph::V(graf)$name, regs_SHH_G3_G4))] = "gold4"

# Cor nos
igraph::V(graf)$nodeColor[!is.na(match(igraph::V(graf)$name, regs_SHH))] = "red"
igraph::V(graf)$nodeColor[!is.na(match(igraph::V(graf)$name, regs_G3))] = "yellow"
igraph::V(graf)$nodeColor[!is.na(match(igraph::V(graf)$name, regs_G4))] = "green"
igraph::V(graf)$nodeColor[!is.na(match(igraph::V(graf)$name, regs_SHH_G3))] = "orange"
igraph::V(graf)$nodeColor[!is.na(match(igraph::V(graf)$name, regs_SHH_G4))] = "pink"
igraph::V(graf)$nodeColor[!is.na(match(igraph::V(graf)$name, regs_G3_G4))] = "yellowgreen"
igraph::V(graf)$nodeColor[!is.na(match(igraph::V(graf)$name, regs_SHH_G3_G4))] = "gold4"



rdp <- RedPort()
calld(rdp)
addGraph(rdp, graf, layout=coords)
relax(rdp, p1=2, p2=100, p3=1, p4=10, p5=100, p6=4, p7=0)


g1 <-  att.addv(g1, to = "nodeFontSize", value = 20)



####################################




