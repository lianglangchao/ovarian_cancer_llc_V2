
library(survival)
library(ggplot2)
library(survminer)
library(dplyr)


rm(list=ls())
choose_genes <- c("ESM1","COL4A1","SPARC","COL4A2","INSR",
                  "SPP1","ANGPT2","RGCC","IGFBP3","HSPG2")


da1 <- read.table("ov_tcga_pub_clinical_data.tsv",header = T,sep = "\t")
da2 <- read.table("ov_tcga_pub.tar/ov_tcga_pub/data_mrna_agilent_microarray.txt",header = T,sep = "\t",check.names = F)
da3 <- da2[da2$Hugo_Symbol %in% choose_genes,]
rownames(da3) <- da3$Hugo_Symbol
da4 <- da3[,-c(1:2)]
da4[11,] <- apply(da4,2,mean)
da5 <- as.data.frame(t(da4))
colnames(da5)[11] <- "ESM1_score"
identical(rownames(da5),da1$Sample.ID)

da1 <- cbind(da1,da5)
colnames(da1)
colnames(da1)[c(15,16)] <- c("OS_MONTHS","OS_STATUS")

da1$OS_STATUS <- substr(da1$OS_STATUS,start = 1,stop = 1)
da1$OS_STATUS <- as.numeric(da1$OS_STATUS)


# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(da1, 
                         time = "OS_MONTHS", 
                         event = "OS_STATUS", 
                         variables = "ESM1_score", 
                         minprop = 0.05)
summary(res.cut)
res.cat <- surv_categorize(res.cut)
head(res.cat)


fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ ESM1_score, data = res.cat)
summary(fit)

p1 <- ggsurvplot(fit,
                 pval = TRUE, conf.int = TRUE,
                 risk.table = TRUE, 
                 risk.table.col = "strata", 
                 linetype = "strata", 
                 surv.median.line = "hv", 
                 ggtheme = theme_bw(),
                 palette = c("#E7B800", "#2E9FDF","#e31a1c","#984ea3"),
)

p1



da1$ESM1_score <- res.cat$ESM1_score
res.cox <- coxph(Surv(OS_MONTHS,OS_STATUS) ~ ESM1_score,data = da1)
res.cox
summary(res.cox)


covariates <- c("Fraction.Genome.Altered", "Neoplasm.Histologic.Grade",  "Mutation.Count", "Platinum.Status", "TMB..nonsynonymous.","Tumor.Stage.2009","ESM1_score")
class(da1$Fraction.Genome.Altered)
class(da1$Neoplasm.Histologic.Grade) #character
class(da1$Mutation.Count)
class(da1$Platinum.Status)           #character
class(da1$TMB..nonsynonymous.)
class(da1$Tumor.Stage.2009)          #character
class(da1$ESM1_score)                #character

table(da1$Neoplasm.Histologic.Grade)
table(da1$Platinum.Status)
table(da1$Tumor.Stage.2009)
table(da1$ESM1_score)

da1$ESM1_score[which(da1$ESM1_score == "low")] <- 1
da1$ESM1_score[which(da1$ESM1_score == "high")] <- 2

da1$Neoplasm.Histologic.Grade[which(da1$Neoplasm.Histologic.Grade == "G2")] <- 1
da1$Neoplasm.Histologic.Grade[which(da1$Neoplasm.Histologic.Grade == "G3")] <- 2

da1$Platinum.Status[which(da1$Platinum.Status == "Resistant")] <- 1
da1$Platinum.Status[which(da1$Platinum.Status == "Sensitive")] <- 2
da1$Platinum.Status[which(da1$Platinum.Status == "Tooearly")] <- 3

da1$Tumor.Stage.2009[which(da1$Tumor.Stage.2009 == "IIA")] <- 1
da1$Tumor.Stage.2009[which(da1$Tumor.Stage.2009 == "IIB")] <- 2
da1$Tumor.Stage.2009[which(da1$Tumor.Stage.2009 == "IIC")] <- 3
da1$Tumor.Stage.2009[which(da1$Tumor.Stage.2009 == "IIIA")] <- 4
da1$Tumor.Stage.2009[which(da1$Tumor.Stage.2009 == "IIIB")] <- 5
da1$Tumor.Stage.2009[which(da1$Tumor.Stage.2009 == "IIIC")] <- 6
da1$Tumor.Stage.2009[which(da1$Tumor.Stage.2009 == "IV")] <- 7

da1$Neoplasm.Histologic.Grade <- as.numeric(da1$Neoplasm.Histologic.Grade)
da1$Platinum.Status <- as.numeric(da1$Platinum.Status)
da1$Tumor.Stage.2009 <- as.numeric(da1$Tumor.Stage.2009)
da1$ESM1_score <- as.numeric(da1$ESM1_score)




univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS_MONTHS,OS_STATUS)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = da1)})
univ_models
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)",           "wald.test", "p.value")
                         return(res)
                       })
class(univ_results)
## [1] "list"
str(univ_results)
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)
write.table(file="univariate_cox_result_1.txt",as.data.frame(res),quote=F,sep="\t")



res.cox <- coxph(Surv(OS_MONTHS,OS_STATUS) ~ ESM1_score + Mutation.Count+Platinum.Status + TMB..nonsynonymous.+Tumor.Stage.2009, data =  da1)
res.cox
summary(res.cox)
fit<- survfit(Surv(OS_MONTHS,OS_STATUS) ~ ESM1_score, data = da1)
ggsurvplot(fit, data = da1)
p1 <- ggsurvplot(fit,
                 pval = TRUE, conf.int = TRUE,
                 risk.table = TRUE, 
                 risk.table.col = "strata", 
                 linetype = "strata",
                 surv.median.line = "hv",
                 ggtheme = theme_bw(), 
                 palette = c("#E7B800", "#2E9FDF","#e31a1c","#984ea3"),
)

p1
pdf("E3_Endo_ESM1_multivariate_cox_1.pdf",width = 7,height = 5)
p1
dev.off()

x <- summary(res.cox)
pvalue=signif(as.matrix(x$coefficients)[,5],2)
HR=signif(as.matrix(x$coefficients)[,2],2)
low=signif(x$conf.int[,3],2)
high=signif(x$conf.int[,4],2)
multi_res=data.frame(p.value=pvalue,
                     HR=paste(HR," (",low,"-",high,")",sep=""),
                     stringsAsFactors = F
)
multi_res
write.table(file="multivariate_cox_result_1.txt",multi_res,quote=F,sep="\t")

pdf("E3_Endo_ESM1_multivariate_cox_forest.pdf",width = 8,height = 8)
ggforest(res.cox,main="hazard ratio",
         cpositions=c(0.02,0.22,0.4),
         fontsize=0.8,noDigits=2)
dev.off()