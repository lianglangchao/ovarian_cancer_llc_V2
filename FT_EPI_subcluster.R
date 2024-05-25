library(SeuratDisk)
library(Seurat)

FT <- readRDS("EPI_0916_allFT.rds")

FT1 <- CreateSeuratObject(counts=FT@assays$RNA@counts)
FT1[["percent.mt"]] <- PercentageFeatureSet(FT1, pattern = "^MT-")
#FT1 <- NormalizeData(FT1, normalization.method = "LogNormalize", scale.factor = 10000)


SaveH5Seurat(FT1, filename = "EPI_0916_allFT.h5Seurat",overwrite=T)
Convert("EPI_0916_allFT.h5Seurat", dest = "h5ad",assay="RNA",overwrite=T)

adata = sc.read_h5ad('EPI_0916_allFT.h5ad')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
notMT= adata.var['mt'][adata.var['mt']==False]
adata = adata[:, notMT.index]
var_names=adata.var.index.tolist()
n_vars=len(var_names)
n_vars
var=pd.DataFrame(data=var_names,index=var_names)
var=var.rename(columns={0: 'index'})
adata1=ad.AnnData(adata.X.astype(np.float32),obs=adata.obs,var=var,dtype='float32')
count_adat_fn='./FT_counts_without_MT.h5ad'
sc.write(count_adat_fn, adata1)
numiter=200
numhvgenes=3000
output_directory = './'
run_name = 'EPI_cNMF'
K = ' '.join([str(i) for i in range(3,26)])
seed = 141
countfn='./FT_counts_without_MT.h5ad'
cnmf_obj.prepare(counts_fn=countfn, components=np.arange(3,26), n_iter=200, seed=141, num_highvar_genes=3000)


#finally select 8#
k=8
W <- read.table("EPI_cNMF/EPI_cNMF.usages.k_8.dt_0_05.consensus.txt",header = T)
colnames(W) <- paste("module",1:k,sep="_")
W[1:5,1:5]
H <- read.table("EPI_cNMF/EPI_cNMF.gene_spectra_score.k_8.dt_0_05.txt",row.names = 1,header = T,sep = "\t")
rownames(H) <- paste("module",1:k,sep="_")
H <- data.frame(t(H))
H <- H[-which(is.na(apply(H,MARGIN = 1,sum))),]
H[1:5,1:5]

hvb_gene <- read.table("EPI_cNMF/EPI_cNMF.overdispersed_genes.txt")
VariableFeatures(FT) <- hvb_gene$V1
FT <- ScaleData(FT)

FT@reductions$nmf <- FT@reductions$pca
FT@reductions$nmf@cell.embeddings <- as.matrix(W) 
FT@reductions$nmf@feature.loadings <- t(H)  

set.seed(219)
FT <- RunUMAP(FT, reduction = 'nmf', dims = 1:k) %>% 
  FindNeighbors(reduction = 'nmf', dims = 1:k) %>% FindClusters(resolution = 0.1)
ncluster <- length(unique(FT$'RNA_snn_res.0.1'))
FT <- RunTSNE(FT, reduction = 'nmf', dims = 1:k)

gc()
color <- c(ggsci::pal_d3("category20")(20),brewer.pal(7, "Dark2")[-5])
pdf("./FT_NMF_recluster_dimplot1.pdf",width = 7,height = 6)
set.seed(199)
print(DimPlot(FT,group.by = "RNA_snn_res.0.1",label = T,cols=sample(color,size = 10,replace = F)))
set.seed(1111)
print(DimPlot(FT,group.by = "Stage1",label = T,cols=sample(color,size = 12,replace = F)))
set.seed(199)
print(DimPlot(FT,group.by = "RNA_snn_res.0.15",label = T,cols=sample(color,size = 10,replace = F)))
set.seed(199)
print(DimPlot(FT,reduction = "tsne",group.by = "RNA_snn_res.0.1",label = T,cols=sample(color,size = length(unique(FT$seurat_clusters)),replace = F)))
set.seed(1111)
print(DimPlot(FT,reduction = "tsne",group.by = "Stage1",label = T,cols=sample(color,size = 12,replace = F)))
set.seed(199)
print(DimPlot(FT,reduction = "tsne",group.by = "RNA_snn_res.0.15",label = T,cols=sample(color,size = 10,replace = F)))
dev.off()


FT <- AddMetaData(FT,W)
wd <- ceiling(sqrt(k))
ht <- ceiling(k/wd)
png("./FT_NMF_umap_topic_score_featureplot.png",width=500*wd,height=500*ht)
print(FeaturePlot(FT,features = paste("module",1:k,sep="_"),ncol = wd,label=T))
dev.off()
try(dev.off())

png("./FT_NMF_tsne_topic_score_featureplot.png",width=500*wd,height=500*ht)
print(FeaturePlot(FT,reduction = "tsne",features = paste("module",1:k,sep="_"),ncol = wd,label=T))
dev.off()
try(dev.off())

png("./FT_NMF_recluster_topic_dotplot.png",width=(ncluster*50+100),height=(50*k+100))
print(DotPlot(FT,features = paste("module",1:k,sep="_"))+
        scale_color_gradient2(low = "#084594",mid = "#ffffb3",high = "#e41a1c",midpoint = 0)+
        coord_flip())
dev.off()
try(dev.off())

H1 <- as.data.frame(t(H))
H1$module <- row.names(H1)
mydata <- melt(H1)
colnames(mydata) <- c("module","gene","value")
mydata %>%
  group_by(module) %>%
  top_n(n = 100, wt = value) -> myselect
myselect1 <- as.data.frame(myselect)
myselect1$gene <- as.character(myselect1$gene)

gid <- bitr(unique(myselect1$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')
markers <- full_join(myselect1, gid, by=c('gene' = 'SYMBOL'))
x = compareCluster(ENTREZID ~ module, data = markers, fun='enrichGO', OrgDb='org.Hs.eg.db',ont="BP")
p1 <- dotplot(x,showCategory=6, label_format=40) + theme(axis.text.x = element_text(angle=45, hjust=1)) 
png("./FT_topic_function_anno.png",width=(100*k+500),height=(150*k+100))
print(p1)
dev.off()
try(dev.off())
write.csv(x,"FT_all_230918_NMF_topic8_enrich.csv")

png("./FT_NMF_umap_cell_class_featureplot.png",width=1000,height=1500)
p2 <- FeaturePlot(FT,features = c("PAX8","KRT7","FOXJ1","CCDC17","CCDC78","CAPS"),ncol = 2,label=T)
print(p2)
dev.off()

saveRDS(FT,"./FT_all_230922_NMF_topic8_res0.1_recluster.rds")
markers <- FindAllMarkers(FT, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers,"FT_all_each_cluster_allDEGs.csv",quote = F)


FT <- readRDS("FT_all_230922_NMF_topic8_res0.1_recluster.rds")
library(dplyr)
DEG <- read.csv("FT_all_each_cluster_allDEGs.csv",row.names=1)
DEG[which(DEG$p_val_adj <= 0.01 & DEG$pct.1 >= 0.5 & DEG$avg_log2FC >= 1),] %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> myselect

myselect1 <- as.data.frame(myselect)
myselect1$gene <- as.character(myselect1$gene)

library(clusterProfiler)
gid <- bitr(unique(myselect1$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')
markers <- full_join(myselect1, gid, by=c('gene' = 'SYMBOL'))
x = compareCluster(ENTREZID ~ cluster, data = markers, fun='enrichGO', OrgDb='org.Hs.eg.db',ont="BP")
p1 <- dotplot(x,showCategory=6, label_format=40) + theme(axis.text.x = element_text(angle=45, hjust=1))
k=8
png("./FT_each_cluster_function_anno.png",width=(100*k+500),height=(150*k+100))
print(p1)
dev.off()



#FT process script#

EPI1 <- readRDS("EPI_230316_NMF_topic9_recluster.rds")
FT <- readRDS("10_OC_initial_analysis/FT_all_230918_NMF_topic8_recluster.rds")
head(FT@meta.data)
EPI1@meta.data$FT_anno <-  EPI1@meta.data$RNA_snn_res.0.15
EPI1@meta.data$FT_anno <- as.character(EPI1@meta.data$FT_anno)
EPI1@meta.data[row.names(FT@meta.data),"FT_anno"] <- paste("FTE",as.character(FT$RNA_snn_res.0.151),sep="")
head(EPI1@meta.data[row.names(FT@meta.data),"FT_anno"])
EPI1@meta.data$FT_anno <- factor(EPI1@meta.data$FT_anno, levels=c("0","1","2","3","4","5","6","7","8","9","10","11","FTE0","FTE1","FTE2","FTE3","FTE4","FTE5","FTE6","FTE7","FTE8","FTE9"))
EPI1 <- AddMetaData(EPI1,FT@meta.data[,paste("module",1:8,sep="_")])


cluster_Palette <- c('#33a02c','#ff69b4','#1f78b4','#ff4500','#b2df8a','#a6cee3','#a020f0','#ff7f00','#66cdaa','#db7093','#1e90ff','#fdbf6f','#6a3d9a','#ffd700','#b15928','#8dd3c7','#ffffb3','#bebada','#ffc1c1','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999','#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd','#fddaec','#f2f2f2','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#ffff33','#a65628','#f781bf','#999999','#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#fddaec','#f2f2f2','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#c51b7d')



pdf("./10_OC_initial_analysis/EPI_FT_reNMF_recluster_dimplot1_split.pdf",width = 14,height = 12)
print(DimPlot(EPI1,group.by = "FT_anno",label = T,cols=cluster_Palette,split.by="FT_anno",ncol=5))
dev.off()

pdf("./10_OC_initial_analysis/EPI_FT_reNMF_recluster_dimplot1.pdf",width = 7,height = 6)
print(DimPlot(EPI1,group.by = "FT_anno",label = T,cols=cluster_Palette))
dev.off()

pdf("./10_OC_initial_analysis/EPI_FT_reNMF_recluster_module_color.pdf",width = 20,height = 20)
FeaturePlot(EPI1,features=c(paste("module",1:8,sep="_")),label=T)
dev.off()

library(progeny)
library(tidyr)
library(tibble)
library(pheatmap)

Idents(EPI1) <- EPI1$FT_anno
CellsClusters <- data.frame(Cell = names(Idents(EPI1)), 
                            CellType = as.character(Idents(EPI1)),
                            stringsAsFactors = FALSE)

## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
EPI1 <- progeny(EPI1, scale=FALSE, organism="Human", top=500, perm=1, 
                return_assay = TRUE)
EPI1@assays$progeny
EPI1 <- Seurat::ScaleData(EPI1, assay = "progeny") 

progeny_scores_df <- 
  as.data.frame(t(GetAssayData(EPI1, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 
dim(progeny_scores_df)


pdf("10_OC_initial_analysis/all_EPI_FT_anno_celltype1_pathway_heatmap_500.pdf",width = 7,height = 7)
progeny_hmap
dev.off()


saveRDS(EPI1@assays$progeny,"all_EPI_FT_anno_celltype1_pathway_progeny_data.rds")


EPI1 <- readRDS("../EPI_230316_NMF_topic9_recluster.rds")
meta <- read.csv("EPI_20230923_FT_anno_result.csv",row.names=1,header=T)
progeny <- readRDS("../all_EPI_FT_anno_celltype1_pathway_progeny_data.rds")

EPI1 <- AddMetaData(EPI1,meta[,c("FT_anno","FT_anno1","FT_anno2")])
EPI1@assays$progeny <- progeny

progeny_scores_df <- 
  as.data.frame(t(GetAssayData(EPI1, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 
dim(progeny_scores_df)

CellsClusters <- data.frame(Cell = names(Idents(EPI1)), 
                            CellCluster = as.character(EPI1@meta.data$FT_anno),
                            CellType1 = as.character(EPI1@meta.data$FT_anno1),
                            CellType2 = as.character(EPI1@meta.data$FT_anno2),
                            stringsAsFactors = FALSE)


## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellCluster) %>%
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


pdf("EPI_FT_anno_cellcluster_pathway_heatmap_500.pdf",width = 7,height = 7)
progeny_hmap
dev.off()





## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType2) %>%
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


pdf("EPI_FT_anno_celltype2_pathway_heatmap_500.pdf",width = 7,height = 7)
progeny_hmap
dev.off()











