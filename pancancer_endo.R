library(Seurat)

obj1 <- readRDS("GSE210347_counts.Rds")
obj2 <- CreateSeuratObject(counts = obj1, min.cells = 0,min.features = 0)
meta1 <- read.table("GSE210347_meta.txt",header = T,sep = "\t")
obj2@meta.data <- meta1
Idents(obj2) <- "celltype"
obj3 <- subset(obj2,idents="Endothelium")

obj3 <- NormalizeData(obj3)
obj3 <- FindVariableFeatures(obj3, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(obj3)
obj3 <- ScaleData(obj3, features = all.genes)
obj3 <- RunPCA(obj3, features = VariableFeatures(object = obj3))
obj3 <- FindNeighbors(obj3, dims = 1:30)
obj3 <- FindClusters(obj3, resolution = 0.5)
obj3 <- RunUMAP(obj3, dims = 1:30)
saveRDS(obj3, file = "pancancer_EC.rds")



setwd("pan_cancer_data")
library(Seurat)

obj1 <- readRDS("pancancer_EC.rds")
rownames(obj1@meta.data) <- colnames(obj1)
obj1 <- FindNeighbors(obj1, dims = 1:30)
obj1 <- FindClusters(obj1, resolution = 0.5)
table(obj1$seurat_clusters)

pdf("EC_dimplot1.pdf",width = 7,height = 6)
DimPlot(obj1,reduction = "umap",group.by = "seurat_clusters",raster = T,label = T)
dev.off()

tip_gene <- c("ESM1","COL4A1","SPARC","COL4A2","IGKC",
              "KDR","FLT1","KIT","CXCR4","ANGPT2")

pdf("EC_vlnplot_1218.pdf",width = 10*4,height = 4*3)
VlnPlot(obj1,group.by = "seurat_clusters",features = tip_gene,pt.size = 0,ncol = 3)
dev.off()

meta1 <- read.table("GSE210347_meta.txt",header = T,sep = "\t")
meta2 <- meta1[meta1$celltype == "Endothelium",]

Idents(obj1) <- "seurat_clusters"
tip14 <- subset(obj1,idents=14)
Idents(tip14) <- "tissue"
tip14 <- subset(tip14,idents=c("bladder","gastric","ICC","lung","PDAC"))


gene1 <- c("EFNA1","HBEGF","HSPG2","SEMA4C","TNFSF10",
           "LAMA4","LAMA5","LAMB1","LAMB2","LAMC1","TGFB1","TGFB2")

pdf("EC_vlnplot2_1218.pdf",width = 4*6,height = 4*2)
VlnPlot(tip14,group.by = "tissue",features = gene1,pt.size = 0,ncol = 6)
dev.off()