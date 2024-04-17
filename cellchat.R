library(CellChat)
library(patchwork)
library(Seurat)
library(ggpubr)
library(cowplot)
options(stringsAsFactors = FALSE)


ov1 <- readRDS("EPI_230316_NMF_topic9_recluster.rds")
ov1

meta1 <- read.csv("10_OC_initial_analysis/EPI_20230923_FT_anno_result.csv",row.names = 1)
ov1@meta.data <- meta1

Idents(ov1) <- "Tissue1"
ov1 <- subset(ov1,idents=c("Fallopian tube"))


endo <- readRDS("OV-Endo-1127.rds")
cell <- merge(ov1,endo)
cell$major <- c(ov1$FT_anno2,endo$celltype2)


cell$major <- as.character(cell$major)
data.input  <- cell
meta <- cell@meta.data
unique(meta$major)

cellchat <- createCellChat(object = data.input, meta = meta,group.by = "major")
cellchat
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use


cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 10) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
table(cellchat@idents)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


pathways.show <- c("TGFb") 
pdf("cellchat_circle_1.pdf",width = 10,height = 10)
vertex.receiver = seq(1,10) # a numeric vector. 
netVisual_aggregate(cellchat,signaling = pathways.show,  vertex.receiver = vertex.receiver)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat,signaling = pathways.show, layout = "circle")
dev.off()

pdf("cellchat_vlnplot_1.pdf",width = 8,height = 6)
plotGeneExpression(cellchat, signaling = "TGFb")
dev.off()