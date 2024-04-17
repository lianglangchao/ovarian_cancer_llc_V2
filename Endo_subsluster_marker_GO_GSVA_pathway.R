
library(Seurat)
library(ggplot2)

endo <- readRDS("06_endothelial/OV-Endo.rds")

meta1 <- read.csv("06_endothelial/Endo_celltype1_meta.csv",row.names = 1)
endo@meta.data <- meta1

DimPlot(endo,group.by = "celltype1",label = T)

endo1 <- endo[,which(endo@assays$RNA@counts["PECAM1",] > 0)]
endo2 <- endo1[,which(endo1@assays$RNA@counts["PTPRC",] == 0)]

DimPlot(endo2,group.by = "celltype1",label = T)


endo2 <- NormalizeData(endo2)
endo2 <- FindVariableFeatures(endo2, selection.method = "vst", nfeatures = 2000)
endo2 <- ScaleData(endo2, verbose = T)
endo2 <- RunPCA(endo2,features = VariableFeatures(object = endo2))


library(harmony)
endo2 <- RunHarmony(endo2, group.by.vars = "Patient_source",max.iter.harmony=15)
endo2 <- RunUMAP(endo2, reduction = "harmony", dims = 1:30)
endo2 <- FindNeighbors(endo2, reduction = "harmony", dims = 1:30) %>% FindClusters()
endo2 <- FindNeighbors(endo2, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 0.5)
endo2 <- FindNeighbors(endo2, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 0.3)

markers1 <- FindAllMarkers(endo2,logfc.threshold = 0.25,min.pct = 0.25,only.pos = T)
write.csv(markers1,"Endo_cluster_markers_1127.csv",quote = F)

gene3 <- c("ACKR1","GJA5","ESM1","PROX1","CA4","MKI67")
pdf("Endo_major_ct_featureplot1.pdf",width = 40,height = 7)
FeaturePlot(endo2,features = gene3,ncol = 6,cols = c("#f0f0f0","#000000"))
dev.off()

pdf("Endo_major_ct_featureplot2.pdf",width = 40,height = 7)
FeaturePlot(endo2,features = gene3,ncol = 6,cols = c("#f0f0f0","#000000"),raster = T)
dev.off()


endo2$celltype2 <- "endo2"
endo2$celltype2[which(endo2$RNA_snn_res.0.8 %in% c(2,3,5,6,7,13,17,19))] = "E1_Endo_ACKR1"
endo2$celltype2[which(endo2$RNA_snn_res.0.8 %in% c(4))] = "E2_Endo_GJA5"
endo2$celltype2[which(endo2$RNA_snn_res.0.8 %in% c(1,11))] = "E3_Endo_ESM1"
endo2$celltype2[which(endo2$RNA_snn_res.0.8 %in% c(0,16))] = "E4_Endo_PROX1"
endo2$celltype2[which(endo2$RNA_snn_res.0.8 %in% c(15))] = "E5_Endo_CA4"
endo2$celltype2[which(endo2$RNA_snn_res.0.8 %in% c(10))] = "E6_Endo_MKI67"
endo2$celltype2[which(endo2$RNA_snn_res.0.8 %in% c(8))] = "E7_Endo_RGS5"
endo2$celltype2[which(endo2$RNA_snn_res.0.8 %in% c(9))] = "E8_Endo_COL1A1"
endo2$celltype2[which(endo2$RNA_snn_res.0.8 %in% c(12))] = "E9_Endo_WFDC2"
endo2$celltype2[which(endo2$RNA_snn_res.0.8 %in% c(14))] = "E10_Endo_IGKC"
endo2$celltype2[which(endo2$RNA_snn_res.0.8 %in% c(18))] = "E11_Endo_CAPS"

endo2$celltype2 <- factor(endo2$celltype2,levels = c("E1_Endo_ACKR1","E2_Endo_GJA5","E3_Endo_ESM1","E4_Endo_PROX1",
                                                     "E5_Endo_CA4","E6_Endo_MKI67","E7_Endo_RGS5","E8_Endo_COL1A1",
                                                     "E9_Endo_WFDC2","E10_Endo_IGKC","E11_Endo_CAPS"))

library(RColorBrewer)
display.brewer.all()

cols1<-brewer.pal(12, "Set3")
mycolors1 <- colorRampPalette(cols1)(12)

pdf("Endo_celltype2_dimplot1.pdf",width = 8,height = 7)
DimPlot(endo2,group.by = "celltype2",label = T,cols = mycolors1,repel = T)
dev.off()

write.csv(endo2@meta.data,"Endo_celltype2_meta.csv",quote = F)
saveRDS(endo2,"OV-Endo-1127.rds")





##############################celltype1_marker_GO##############################
library(org.Hs.eg.db)
library(clusterProfiler)

Idents(endo2) <- "celltype2"
markers1 <- FindAllMarkers(endo2,logfc.threshold = 0.25,min.pct = 0.25,only.pos = T)
write.csv(markers1,"Endo_celltype2_markers.csv",quote = F)


marker1 <-read.csv("Endo_celltype2_markers.csv",row.names=1)
dir.create("Endo_celltype2_GO")
for (i in unique(marker1$cluster)) {
  col<-c("SYMBOL","ENTREZID")
  marker2 <- marker1[marker1$cluster==i,]
  marker3 <- marker2[marker2$p_val_adj < 0.05,]
  da1 <- AnnotationDbi::select(org.Hs.eg.db,columns=col,keytype="SYMBOL",keys=as.character(marker3$gene))
  ego<-enrichGO(OrgDb = "org.Hs.eg.db",gene=da1$ENTREZID,ont="BP",pvalueCutoff=0.05,readable=TRUE)
  write.csv(ego,paste("./Endo_celltype2_GO/",i,"_GO.csv",sep=""),row.names = F)
  
}

GO2 <- NULL
for (i in unique(marker1$cluster)) {
  GO1 <- read.csv(paste("./Endo_celltype2_GO/",i,"_GO.csv",sep=""))
  GO1$celltype <- i
  GO2 <- rbind(GO1[1:5,],GO2)
}

GO2$Description <- factor(GO2$Description,levels = unique(GO2$Description))
GO2$celltype <- factor(GO2$celltype,levels = c("E1_Endo_ACKR1","E2_Endo_GJA5","E3_Endo_ESM1","E4_Endo_PROX1",
                                               "E5_Endo_CA4","E6_Endo_MKI67","E7_Endo_RGS5","E8_Endo_COL1A1",
                                               "E9_Endo_WFDC2","E10_Endo_IGKC","E11_Endo_CAPS"))

p1 <- ggplot(GO2, aes(x = celltype, y = Description,color=-log10(p.adjust),size=Count)) +
  geom_point()+
  scale_color_gradient(low = "#fcbba1",high = "#99000d")+
  labs(x="",y="")+
  theme_bw()+
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45))

pdf("Endo_celltype2_marker_GO1.pdf",width = 10,height = 10)
print(p1)
dev.off()




GO1 <- read.csv("./Endo_celltype2_GO/E3_Endo_ESM1_GO.csv")
GO1$celltype <- "E3_Endo_ESM1"
GO2 <- read.csv("./Endo_celltype2_GO/E11_Endo_CAPS_GO.csv")
GO2$celltype <- "E11_Endo_CAPS"
GO3 <- rbind(GO1[1:10,],GO2[1:10,])
GO3$Description <- factor(GO3$Description,levels = unique(GO3$Description))
GO3$celltype <- factor(GO3$celltype,levels = c("E3_Endo_ESM1","E11_Endo_CAPS"))

p1 <- ggplot(GO3, aes(x = celltype, y = Description,color=-log10(p.adjust),size=Count)) +
  geom_point()+
  scale_color_gradient(low = "#fcbba1",high = "#99000d")+
  scale_size_continuous(range=c(5,10))+
  labs(x="",y="")+
  theme_bw()+
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45))

pdf("Endo_celltype2_marker_GO2.pdf",width = 5.5,height = 7)
print(p1)
dev.off()






######################celltype1_GSVA_Hallmark-----------------------

library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)


genesets <- msigdbr(species = "Homo sapiens", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)

Idents(endo2) <- "celltype2" 
expr <- AverageExpression(endo2, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  
expr <- as.matrix(expr)


gsva.res <- gsva(expr, genesets, method="ssgsea") 
saveRDS(gsva.res, "endo_celltype2_gsva.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "endo_celltype2_gsva.csv", row.names = F)

pdf("endo_celltype2_gsva_heatmap1.pdf",width = 8,height = 12)
pheatmap::pheatmap(gsva.res, show_colnames = T, scale = "row",cluster_cols = F)
dev.off()

######################celltype1_dotplot-----------------------

gene3 <- c("ACKR1","GJA5","ESM1","PROX1",
           "CA4","MKI67","RGS5",
           "COL1A1","WFDC2","IGKC",
           "CAPS","PECAM1","PTPRC")

endo2$celltype2 <- factor(endo2$celltype2,levels = levels(endo2$celltype2)[11:1])

pdf("Endo_celltype2_dotplot1.pdf",width = 9,height = 4)
DotPlot(endo2,group.by = "celltype2",features = gene3)+
  theme(axis.text.x = element_text(face = "italic",vjust = 1, hjust = 1, angle = 45))+
  scale_color_gradient2(low = "#2166AC",mid = "#F7F7F7",high = "#B2182B")
dev.off()

#############################################################




#######################barplot-----------------------------------
endo2$celltype2 <- factor(endo2$celltype2,levels = levels(endo2$celltype2)[11:1])

da1 <- endo2@meta.data[,c("celltype2","Tissue1")]
da2 <- as.data.frame(table(da1$celltype2,da1$Tissue1))
colnames(da2)[1:2] <- c("celltype2","Tissue1")

p1 = ggplot( da2, aes( x = Tissue1, y = Freq, fill = celltype2))+
  geom_bar( position = "fill",stat = "identity",width = 0.7)+theme_bw()+
  scale_fill_manual(values=mycolors1) + 
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf("endo_celltype2_tissue1_barplot1.pdf",width = 8,height = 6)
p1
dev.off()

da3 <- da2[da2$Tissue1 %in% c("Fallopian tube","Normal Ovarium","Ovarium"),]
p1 = ggplot( da3, aes( x = Tissue1, y = Freq, fill = celltype2))+
  geom_bar( position = "fill",stat = "identity",width = 0.6)+theme_bw()+
  scale_fill_manual(values=mycolors1) + 
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf("endo_celltype2_tissue1_barplot2.pdf",width = 5,height = 6)
p1
dev.off()






###################celltype marker heatmap################################################
library(dplyr)

marker1 <- read.csv("Endo_celltype2_markers.csv",row.names = 1)

marker1 %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

da1 <- AverageExpression(endo2,features = top5$gene,group.by = "celltype2")[[1]]

library(pheatmap)
pdf("Endo_celltype2_top5_heatmap1.pdf",width = 11,height = 4)
pheatmap(t(da1),border_color = NA,cluster_cols = F,cluster_rows = F,scale = "column",color = colorRampPalette(colors = c("#01665e","#f5f5f5","#8c510a"))(100))
dev.off()



######################celltype1_pathway_500-----------------------
endo <- endo2

library(progeny)
library(tidyr)
library(tibble)
library(pheatmap)

Idents(endo) <- "celltype2"
CellsClusters <- data.frame(Cell = names(Idents(endo)), 
                            CellType = as.character(Idents(endo)),
                            stringsAsFactors = FALSE)

## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
endo <- progeny(endo, scale=FALSE, organism="Human", top=500, perm=1, 
                return_assay = TRUE)
endo@assays$progeny
## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
endo <- Seurat::ScaleData(endo, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(endo, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 
dim(progeny_scores_df)

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
dim(summarized_progeny_scores)
## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("#2166ac", "#f7f7f7","#b2182b"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df),fontsize=12, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (500)", angle_col = 45,
                        treeheight_col = 0,  border_color = NA)

pdf("Endo_celltype2_pathway_heatmap_500.pdf",width = 6.5,height = 6)
progeny_hmap
dev.off()

