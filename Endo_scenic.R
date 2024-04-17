
rm(list=ls())
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(SummarizedExperiment)


library(SCENIC)
loom <- open_loom('out_SCENIC.loom') 

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')

close_loom(loom)

rownames(regulonAUC)
names(regulons)


meta1 <- read.csv("Endo_celltype2_meta.csv",row.names = 1)


sub_regulonAUC <- regulonAUC[,match(rownames(meta1),colnames(regulonAUC))]
dim(sub_regulonAUC)
meta1 
identical(colnames(sub_regulonAUC), rownames(meta1))

cellClusters <- data.frame(row.names = rownames(meta1), 
                           seurat_clusters = as.character(meta1$celltype2))
cellTypes <- data.frame(row.names = rownames(meta1), 
                        celltype = meta1$celltype2)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 

head(cellTypes) 
sub_regulonAUC[1:4,1:2] 
dim(sub_regulonAUC)



selectedResolution <- "celltype" 
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution]) 


sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)

regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))


regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 

dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)



pheatmap(regulonActivity_byGroup_Scaled)
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellTypes[colnames(sub_regulonAUC), 
                                        selectedResolution]) 
rss=na.omit(rss) 
rssPlot <- plotRSS(rss)


library(dplyr) 
rss=regulonActivity_byGroup_Scaled
head(rss)
library(dplyr) 
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss),
                 cluster =   colnames(rss)[i],
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster) 
n = rss[top5$path,] 


pdf("Endo_1205-TF-2.pdf",width = 6,height = 10)
pheatmap(n,
         annotation_row = rowcn,
         show_rownames = T,border_color = NA)
dev.off()

