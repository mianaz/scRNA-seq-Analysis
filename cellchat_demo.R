# installation
# dependencies
#install.packages('NMF')
#devtools::install_github("jokergoo/circlize")
#devtools::install_github("jokergoo/ComplexHeatmap")
#reticulate::py_install("umap-learn")
#devtools::install_github("jinworks/CellChat")

#load libraries
library(Seurat)
library(CellChat)
library(dplyr)
library(ggplot2)
library(patchwork)
options(stringsAsFactors = FALSE)

#load data
chen_crpc <- readRDS("~/Desktop/chen_crpc.rds")
#we have samples of 2 pathologies: crpc and prad
#table(chen_cohort@meta.data$pathology)
#let's first just use the crpc subset (1 patient)
#chen_prad <- subset(chen_cohort, pathology=="PRAD")
#chen_crpc <- subset(chen_cohort, pathology=="CRPC")
#make sure we have cell meta data
chen_crpc$celltype <- Idents(chen_crpc)
table(chen_crpc$celltype)
#saveRDS(chen_crpc, "~/Desktop/chen_crpc.rds")

# create cellchat object
#cellchat <- createCellChat(chen_crpc, group.by = "celltype", assay = "SCT")
#generate input (normalized data)
data.input=chen_crpc[["SCT"]]@data # normalized data
meta=data.frame(labels=chen_crpc$celltype, row.names=names(chen_crpc$celltype))
meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.human # change to mouse for mouse data
showDatabaseCategory(CellChatDB)
#exclude non-protein signaling 
CellChatDB.use <- subsetDB(CellChatDB)
# if using all signaling pathways,use
# CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use

# subset to signaling genes only
cellchat <- subsetData(cellchat)
#future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#project data to protein-protein network (optional)
cellchat <- projectData(cellchat, PPI.human)
#compute communication probability (this is the longest step)
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = FALSE) # change raw.use to FALSE when using PPI
#filter out cell-cell communication in types with <10 cells
cellchat <- filterCommunication(cellchat, min.cells = 5)
# infer pathway
cellchat <- computeCommunProbPathway(cellchat)
# create net
cellchat <- aggregateNet(cellchat)
# visualize network
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
