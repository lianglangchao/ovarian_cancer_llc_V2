library(infercnv)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(monocle)


EPI1 <- readRDS("EPI_0209_NMF_cluster.rds")
EPI1@meta.data[which(is.na(EPI1$percent.mt)),"percent.mt"] <- 0
EPI1@meta.data[which(is.na(EPI1$Stage2)),"Stage2"] <- "Unknown"
set.seed(12345)
select_CNV_cell <- sample(colnames(EPI1),30000,replace = F)

EPI1 <- EPI1[,select_CNV_cell]
color <- colorRampPalette(c(brewer.pal(8, "Accent"),brewer.pal(8, "Set3"),brewer.pal(8, "Dark2")))(64)

pdf("./05_EPI_monocle2_test/EPI_3w_select_dimplot1.pdf",width = 7,height = 6)
set.seed(123)
DimPlot(EPI1,group.by = "RNA_snn_res.0.7",label = T,cols=sample(color,size = 40,replace = F))
dev.off()

data <- as(as.matrix(EPI1@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = EPI1@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.1,
                              expressionFamily = negbinomial.size())
HSMM <- monocle_cds

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

dim(fData(HSMM))

#Filtering low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 3 )
print(head(fData(HSMM)))
print(head(pData(HSMM)))

expressed_genes <- rownames(EPI1)
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~percent.mt")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))

HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)


#Trajectory step 2: reduce data dimensionality  
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')



HSMM <- orderCells(HSMM)

a <- pData(HSMM)
write.csv(a,"./05_EPI_monocle2_test/pdata.csv",quote=F)

pdf("./05_EPI_monocle2_test/EPI_select_Pseudotime_0213.pdf",width = 6,height = 4)
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
dev.off()

pdf("./05_EPI_monocle2_test/EPI_select_Stage_0213.pdf",width = 20,height = 4)
plot_cell_trajectory(HSMM, color_by = "Stage2")+ facet_wrap("~Stage2", nrow = 2)
dev.off()

pdf("./05_EPI_monocle2_test/EPI_select_state_0213.pdf",width = 6,height = 4)
plot_cell_trajectory(HSMM, color_by = "State")
dev.off()

saveRDS(HSMM,"./05_EPI_monocle2_test/EPI_select_monocle2_0213.rds")

 EPI <- readRDS("EPI_0227_NMF_cluster.rds")
HSMM <- readRDS("./05_EPI_monocle2_test/EPI_select_monocle2_0213.rds")

to_be_tested <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% EPI@assays$RNA@var.features))
cds_subset <- HSMM[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")

head(diff_test_res[,c("gene_short_name", "pval", "qval")])


diff_test_res <- differentialGeneTest(HSMM[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))

saveRDS(diff_test_res,"./05_EPI_monocle2_test/Pseudotime_diff_gene_result.rds")
saveRDS(cds_subset,"./05_EPI_monocle2_test/Pseudotime_diff_gene_sub_object.rds")
saveRDS(HSMM,"./05_EPI_monocle2_test/EPI_select_monocle2_0323.rds")


HSMM@phenoData@data$Cluster <- EPI1@meta.data[row.names(HSMM@phenoData@data),"RNA_snn_res.0.15"]

to_be_tested <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% EPI1@assays$RNA@var.features))
cds_subset@phenoData@data$Cluster <- EPI1@meta.data[row.names(cds_subset@phenoData@data),"RNA_snn_res.0.15"]

cds_subset <- cds_subset[to_be_tested,pData(cds_subset)$Cluster %in% c(1,2,3)]
diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01))

sig_gene_names <- row.names(top_n(diff_test_res,n=100,wt=LogQ))

pdf("Fig1F_EPI_monocle_Heatmap_1.pdf",width=8,height = 4)
plot_pseudotime_heatmap(cds_subset[sig_gene_names,],
                        num_clusters = 5,
                        cores = 1,
                        show_rownames = F)
dev.off()


GO_results <- read.csv("GO_AllLists.csv")
GO_results %>%
  group_by(GeneList) %>%
  top_n(n = 5, wt = -LogP) -> myselect
GO_results %>%
  group_by(GeneList) %>%
  top_n(n = 10, wt = -LogP) -> myselect1

write.csv(GO_results[which(GO_results$Description %in% myselect1$Description),],"/Users/chaichaochao/Documents/BGI/1_研究生/2_华大/1_项目/18_妇科肿瘤数据挖掘/11_OV2_figure/Figure1_related/table/TableS1.2_monocle_function_enrich.csv",row.names = F)

mydata <- reshape2::dcast(GO_results[which(GO_results$Description %in% myselect$Description),],formula = Description~GeneList,value.var = "LogP")
mydata[is.na(mydata)] <- 0
row.names(mydata) <- mydata[,"Description"]
mydata <- mydata[,-1]
mydata<- -mydata
pheatmap::pheatmap(mydata)

pdf("Fig1G_EPI_monocle_GO_enrich_Heatmap.pdf",width=6,height = 5)
pheatmap::pheatmap(mydata)
dev.off()



