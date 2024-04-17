
library(Seurat)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
library(CIBERSORT)
source("sourcecibersort.R")
library(org.Hs.eg.db)
col<-c("ENTREZID","SYMBOL")

data1 <- read.table("./GSE223426_bulk_RNAseq_rawcounts.txt",sep = "\t",header = T,fill = T,check.names = F)
data1 <- data1[,-28]
colnames(data1)[1] <- "geneID"

da1 <- AnnotationDbi::select(org.Hs.eg.db,columns=col,keytype="ENTREZID",keys=as.character(data1$geneID))
da2 <- da1[!duplicated(da1$SYMBOL),]
da3 <- na.omit(da2)
data2 <- data1[data1$geneID %in% da3$ENTREZID,]
rownames(data2) <- data2$geneID
data2 <- data2[da3$ENTREZID,]
table(data2$geneID == da3$ENTREZID)
data2$geneID <- da3$SYMBOL
write.table(data2,"GSE223426_bulk_RNAseq_rawcounts_process.txt",sep = "\t",row.names = F,quote = F)



endo <- readRDS("OV-Endo-1127.rds")
Idents(endo) <- "celltype2"
vargene1 <- VariableFeatures(endo)
vargene2 <- read.csv("d:/project/OV2/09_endothelial_231127/Endo_celltype2_markers.csv")
vargene3 <- vargene2[vargene2$p_val_adj < 0.01 & vargene2$avg_log2FC > 1,]
vargene4 <- unique(vargene3$gene)
da1 <-AverageExpression(endo,features = vargene4,assays = "RNA",slot = "count")[[1]]
write.table(da1,"Endo-marker-0412-1.txt",row.names = T,sep = "\t",quote = F)



TME.results = CIBERSORT("Endo-marker-0412-1.txt", 
                        "GSE223426_bulk_RNAseq_rawcounts_process.txt" , 
                        perm = 1000, 
                        QN = T)
TME.results
da1 <- as.data.frame(TME.results[,1:14])
da1$group <- c(rep("Ovarian tumor",10),rep("Fallopian tube",12),rep("Normal Ovarium",4))
da1$group <- factor(da1$group,levels = c("Fallopian tube","Normal Ovarium","Ovarian tumor"))

library(reshape2)
da2 <- melt(da1,id.vars = "group")
library(ggpubr)
library(ggplot2)
compaired <-list(c("Normal Ovarium","Fallopian tube"),c("Normal Ovarium","Ovarian tumor"),c("Fallopian tube","Ovarian tumor"))
pall <- list()
for (i in unique(da2$variable)) {
  #i= "E3_Endo_ESM1"
  
  da3 <- da2[da2$variable== i,]
  
  p1 <- ggplot(da3,mapping = aes(x = group,y = value,color=group))+
    geom_boxplot(position=position_dodge(1))+ggtitle(i)+
    scale_color_manual(values = c("#3288bd","#4daf4a","#d53e4f"))+
    #描点
    geom_dotplot(aes(fill=group),binaxis = "y", stackdir = "center",
                 position = position_dodge(1))+
    scale_fill_manual(values = c("#3288bd","#4daf4a","#d53e4f"))+
    stat_compare_means(comparisons=compaired,method = "t.test",label = "p.signif")+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45,vjust = 1,hjust = 1))
  
  pall[[i]] <- p1 
}

library(cowplot)
pdf("boxplot1_0412_1000_1.pdf",width = 20,height = 12)
plot_grid(plotlist = pall,ncol = 5)
dev.off()



