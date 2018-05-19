#######################################################
# Deseq Analysis for Naive_Rats (Taconic vs Envigo)   #
#######################################################

#loading libraries (computerome)
library("ggplot2")
library("reshape2")
library("DESeq2")
library("plyr")
library("pheatmap")
library("RColorBrewer")
library("WriteXLS")
library("regionReport")
library("ggfortify")

getwd()

# load sample table computerome.
load("/home/people/prasing/Master_Thesis/scripts/sample_table_naive.R")


# Directory with HTSeq count data for Naive rats
# HTSeq count directory for computerome
HTSeq_direc <- "/home/people/prasing/Master_Thesis/Datasets/study1_naive_diff"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# divide into TG and DRG
sample_table_TG <- sample_table[grep("TG",sample_table$sampleName),]

DEseq_data_TG <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table_TG, directory = HTSeq_direc, design = ~ breeder)
DEseq_data_TG
dim(DEseq_data_TG)
# Note:
# Here each row is gene and each column is a sample name. same for DRG to
#[1] 32662    12


# remove genes that have 0 to 1 count for all samples. 
filt_DEseq_data_TG <- DEseq_data_TG[ rowSums(counts(DEseq_data_TG))>1,]
dim(filt_DEseq_data_TG)
#[1] 22102    12

# Extracting the count dataset.
cou_data <- as.data.frame(counts(filt_DEseq_data_TG, normalized = FALSE))
dim(cou_data)


#Caculate the size_factor, normalization
# TG
Norm_filt_DEseq_data_TG <- estimateSizeFactors(filt_DEseq_data_TG)
dim(Norm_filt_DEseq_data_TG)
#[1] 22102    12

# Extracting the count dataset.
# filt_norm_TG_cou <- as.data.frame(counts(Norm_filt_DEseq_data_TG, normalized = TRUE))
# dim(filt_norm_TG_cou)

# Extracting counts per sample. 
filt_TG_cou_Taco <- cou_data[,1:6]
filt_TG_cou_Envi <- cou_data[,7:12]
 
 
#"""""""""""" Calcualting coffeficient of variation for (TG) tissue. before DE testing """"""""""""""
# Taco
cv_TG_taco_cou <- apply(filt_TG_cou_Taco ,1, function(x) ((sd(x)/mean(x))*100))
cv_melted_data_TG_taco <- melt(cv_TG_taco_cou)

# plotting cv without squaered
# TG_taco
p_cv_TG_taco <-  ggplot(cv_melted_data_TG_taco ,aes(x = value)) + geom_density()+
  xlab("Coefficient of Variation") +
  ylab("Density")+ ggtitle("Density-plot of Cofficient of Variation\n Naive Rats Taconic (TG)")+ theme_bw()

p_cv_TG_taco

# Envi
cv_TG_envi_cou <- apply(filt_TG_cou_Envi ,1, function(x) ((sd(x)/mean(x))*100))
cv_melted_data_TG_envi <- melt(cv_TG_envi_cou)
 
#TG_envi
p_cv_TG_envi <- ggplot(cv_melted_data_TG_envi ,aes(x = value)) + geom_density()+
  xlab("Coefficient of Variation") +
  ylab("Density")+ ggtitle("Density-plot of Cofficient of Variation\n Naive Rats Envigo (TG)")+ theme_bw()

p_cv_TG_envi
dev.off()

# calculating means (TG)
# Taco
calc_means_TG_taco <- rowMeans(filt_TG_cou_Taco, na.rm = TRUE, dims = 1)
calc_log2_means_TG_taco <- log2(calc_means_TG_taco)

# plotting calc means
# taco
qplot(calc_log2_means_TG_taco , cv_TG_taco_cou, xlab = "log2(means)", ylab = "CV", main = "CV versus log2(mean) plot of Taconic (TG)") + geom_smooth()
dev.off()

# envigo
calc_means_TG_envi <- rowMeans(filt_TG_cou_Envi, na.rm = TRUE, dims = 1)
calc_log2_means_TG_envi <- log2(calc_means_TG_envi)
 

# For envigo
qplot(calc_log2_means_TG_envi  , cv_TG_envi_cou , xlab = "log2(means)", ylab = "CV", main = "CV versus log2(mean) plot of Envigo (TG)") + geom_smooth()
dev.off()

# Compute the density (TG) .
# taco
dens_taco_TG <- density(calc_log2_means_TG_taco)
 
# plotting the density of means
# For Taconic
plot(dens_taco_TG, type="l", col="green",
     xlab="x", main="Density plot of the log2(mean) of Taconic (TG)",lwd=2)

# envigo
dens_envi_TG <- density(calc_log2_means_TG_envi)

# Combining both suppliers. 
plot(dens_envi_TG, type="l", col="red",
    xlab="x", main="Density plot of the log2(mean) of Envigo (TG)",lwd=2)
 
plot(dens_envi_TG, col= "red", lwd = 3, main="Density plot of the log2(mean) of Envigo & Taconic (TG)")
lines(dens_taco_TG, col = "blue", lwd = 3)
legend('topleft', c("taconic", "envigo"), lwd=c(2.5,2.5),col=c("blue", "red"), cex = 0.75)


# computing sample distances to see how similar are the samples to each other. 
rld <- rlog(cou_data, blind=FALSE)
head(assay(rld), 3)
 
# Making the PCA plots
data <- plotPCA(rld, intgroup=c("breeder"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color = breeder, shape = breeder, label = name)) +
   geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text()
dev.off()
 

#Estimating dispersions. 
disper_DEseq_data_TG <- estimateDispersions(Norm_filt_DEseq_data_TG)
dim(disper_DEseq_data_TG)
#[1] 22102    12

# plot dispersion ('This is the dispersion plot for TG both breeds')
# TG
plotDispEsts(disper_DEseq_data_TG)


# doing the negative binomial wald test
# TG
nbinom_DEseq_data_TG <- nbinomWaldTest(disper_DEseq_data_TG)
dim(nbinom_DEseq_data_TG)
#[1] 22102    12

# Adjusting the p-value to 0.05
# TG
res_TG_new_Pval <- results(nbinom_DEseq_data_TG, alpha=0.05)
res_TG_new_Pval
summary(res_TG_new_Pval)
# log2 fold change (MAP): breeder Taconic vs Envigo 
# Wald test p-value: breeder Taconic vs Envigo 
# DataFrame with 22102 rows and 6 columns
dim(res_TG_new_Pval)
#[1] 22102     6


# How many genes adjusted p-values were less than 0.05?
# TG
sum(res_TG_new_Pval$padj < 0.05, na.rm = TRUE)
#[1] 7188 differentially expressed. 7186 (DF genes now)

sig_gene_TG <- subset(res_TG_new_Pval, res_TG_new_Pval$padj < 0.05)
sig_gene_TG$geneID <- rownames(sig_gene_TG)
dim(sig_gene_TG)
#[1] 7188    7
summary(sig_gene_TG)
# out of 7188 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 3471, 48% 
# LFC < 0 (down)   : 3717, 52% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 0, 0% 
# (mean count < 1)

# this is to get the information about  result section of deseq. 
mcols(sig_gene_TG, use.names = TRUE)

# Helps to find out which gene is higher expressed in the breeders. 
topGene <- rownames(sig_gene_TG)[which.min(sig_gene_TG$padj)]
plotCounts(nbinom_DEseq_data_TG, gene = topGene, intgroup="breeder", returnData = FALSE)

# Helps plotting multiple gene counts according to their groups. 
par(mfrow=c(1,5))

for (i in 1:5){
  gene <- rownames(sig_gene_DRG)[i]
  main = gene
  DESeq2::plotCounts(nbinom_DEseq_data_DRG, gene=gene, intgroup="breeder", main = main)
}

#################### Merging DE genes with original DB 22102 genes. ################

# extracting the count genes
filt_norm_TG_cou <- as.data.frame(counts(Norm_filt_DEseq_data_TG, normalized = TRUE))
filt_norm_TG_cou$geneID <- rownames(filt_norm_TG_cou)


# Extracting counts per sample.
# Taco
all_TG_genes_Taco <- filt_norm_TG_cou[,1:6]
all_TG_genes_Taco$geneID <- filt_norm_TG_cou$geneID
all_TG_genes_Taco$log2foldchange <- res_TG_new_Pval$log2FoldChange
all_TG_genes_Taco$padj <- res_TG_new_Pval$padj

# Calculating CV after de testing.
# Taco
DE_TG_taco_cou <- apply(all_TG_genes_Taco[,1:6] ,1, function(x) ((sd(x)/mean(x))*100))
DE_TG_taco_cou_melted <- as.data.frame(DE_TG_taco_cou)
all_TG_genes_Taco$CV <- DE_TG_taco_cou_melted$DE_TG_taco_cou

# Calculating Means after de testing.
# Taco
DE_means_TG_taco <- rowMeans(all_TG_genes_Taco[,1:6], na.rm = TRUE, dims = 1)
DE_log2_means_TG_taco <- log2(DE_means_TG_taco)
DE_log2_means_TG_taco <- as.data.frame(DE_log2_means_TG_taco)
all_TG_genes_Taco$means <- DE_log2_means_TG_taco$DE_log2_means_TG_taco

# Extracting counts per sample.
# Envi
all_TG_genes_Envigo <- filt_norm_TG_cou[,7:12]
all_TG_genes_Envigo$geneID <- filt_norm_TG_cou$geneID
all_TG_genes_Envigo$log2foldchange <- res_TG_new_Pval$log2FoldChange
all_TG_genes_Envigo$padj <- res_TG_new_Pval$padj

# Calculating CV after de testing.
# Envi
DE_TG_envi_cou <- apply(all_TG_genes_Envigo[,1:6] ,1, function(x) ((sd(x)/mean(x))*100))
DE_melted_data_TG_envi <- melt(DE_TG_envi_cou)
all_TG_genes_Envigo$CV <- DE_melted_data_TG_envi$value

# Calculating Means after de testing.
# envigo
DE_means_TG_envi <- rowMeans(all_TG_genes_Envigo[,1:6], na.rm = TRUE, dims = 1)
DE_log2_means_TG_envi <- log2(DE_means_TG_envi)
DE_log2_means_TG_envi <- as.data.frame(DE_log2_means_TG_envi)
all_TG_genes_Envigo$means <- DE_log2_means_TG_envi$DE_log2_means_TG_envi


res_TG_taco_envi <- as.data.frame(res_TG_new_Pval)
res_TG_taco_envi$gene_id <- rownames(res_TG_taco_envi)
res_TG_taco_envi$cv_taco <- all_TG_genes_Taco$CV
res_TG_taco_envi$mean_taco <- all_TG_genes_Taco$means
res_TG_taco_envi$cv_envi <- all_TG_genes_Envigo$CV   
res_TG_taco_envi$mean_envi <- all_TG_genes_Envigo$means
res_TG_taco_envi <- subset(res_TG_taco_envi, res_TG_taco_envi$padj < 0.05)
res_TG_taco_envi <- res_TG_taco_envi[(res_TG_taco_envi$cv_taco <75 & res_TG_taco_envi$cv_envi <75), ]
res_TG_taco_envi <- res_TG_taco_envi[which(res_TG_taco_envi$mean_taco >5 | res_TG_taco_envi$mean_envi >5), ]
res_TG_taco_envi_up <- subset(res_TG_taco_envi, res_TG_taco_envi$log2FoldChange > 1)
res_TG_taco_envi_down <- subset(res_TG_taco_envi, res_TG_taco_envi$log2FoldChange < -1)

dim(res_TG_taco_envi)
dim(res_TG_taco_envi_up)
dim(res_TG_taco_envi_down)
View(res_TG_taco_envi)

#Identfying most significant up-reg genes:
sig_reg_tg_taco_envi_up <- as.data.frame(res_TG_taco_envi_up$gene_id)
View(sig_reg_tg_taco_envi_up)

write.table(sig_reg_tg_taco_envi_up, file = "sig_reg_tg_up_taco-envi.txt", sep = "\t", col.names = FALSE , row.names = FALSE, quote = FALSE)

sig_reg_tg_names_up <- read.table("sig_tg_taco_envi_up_info.txt", header = TRUE, sep = ",")
colnames(sig_reg_tg_names_up)[1] <- "gene_id"
View(sig_reg_tg_names_up)

sig_reg_tg_up <- merge.data.frame(res_TG_taco_envi_up, sig_reg_tg_names_up, by = "gene_id", all = TRUE)
dim.data.frame(sig_reg_tg_up)
View(sig_reg_tg_up)

sorted_sig_reg_tg_up <- sig_reg_tg_up[order(-sig_reg_tg_up$log2FoldChange),c(1,3,7,9,11,12,13,14)]
View(sorted_sig_reg_tg_up)
write.table(sorted_sig_reg_tg_up, file = "sig_tg_up_info.txt", sep = "\t", col.names = TRUE , row.names = FALSE, quote = FALSE)

#Identfying most significant down-reg genes:
sig_reg_tg_taco_envi_down <- as.data.frame(res_TG_taco_envi_down$gene_id)
View(sig_reg_tg_taco_envi_down)

write.table(sig_reg_tg_taco_envi_down, file = "sig_reg_tg_down_taco-envi.txt", sep = "\t", col.names = FALSE , row.names = FALSE, quote = FALSE)

sig_reg_tg_names_down <- read.table("sig_tg_taco_envi_down.txt", header = TRUE, sep = ",")
colnames(sig_reg_tg_names_down)[1] <- "gene_id"

sig_reg_tg_down <- merge.data.frame(res_TG_taco_envi_down, sig_reg_tg_names_down, by = "gene_id", all = TRUE)
dim.data.frame(sig_reg_tg_down)
View(sig_reg_tg_down)

sorted_sig_reg_tg_down <- sig_reg_tg_down[order(sig_reg_tg_down$log2FoldChange),c(1,3,7,9,11,12,13,14)]
View(sorted_sig_reg_tg_down)
write.table(sorted_sig_reg_tg_down, file = "sig_tg_down_info.txt", sep = "\t", col.names = TRUE , row.names = FALSE, quote = FALSE)


# volcano plot
#sig=ifelse(res_TG_taco_envi$padj<0.05, "FDR<0.05","no_sig")
qplot(res_TG_taco_envi$log2FoldChange  , res_TG_taco_envi$cv_envi , xlab = "log2_fold_change", ylab = "CV", main = "CV versus log2(foldchange) plot of Taconic\n genes (TG)")+ 
  annotate("text", x = -2, y = 200, label = paste0("",dim(res_TG_taco_envi)[1]," DE genes"), color="red")
dev.off()

# Mean and cv distribution plot. 
qplot(res_TG_taco_envi$mean_envi , res_TG_taco_envi$cv_envi , xlab = "log2(means)", ylab = "CV", main = "CV versus log2(mean) plot of Envigo\n genes (TG)") + geom_smooth() + 
  annotate("text", x = 16, y = 200, label = paste0("",dim(res_TG_taco_envi)[1]," genes"), color="red")


# Making the Venn Diagram:
#library("gplots", lib.loc="C:/Program Files/R/R-3.2.3/library")
library("gplots")


# up-regulated
venn_db_up <- list(res_TG_taco_envi_up$gene_id, res_DRG_taco_envi_up$gene_id)

venn_up <- venn(venn_db_up, show.plot = FALSE)
venn_up

plot.venn(venn_up)
# Note: Here A = up_sig_gene_TG & B = up_sig_gene_DRG
# Note: Out of 104 up_sig_gene_TG and 164 up_sig_gene_DRG, there are 6 common genes. 
#Note: Common genes = "ENSRNOG00000002730" "ENSRNOG00000005811" "ENSRNOG00000009589" "ENSRNOG00000038406" "ENSRNOG00000058039" "ENSRNOG00000061179"



# down-regulated
venn_db_down <- list(res_TG_taco_envi_down$gene_id, res_DRG_taco_envi_down$gene_id)

venn_down <- venn(venn_db_down, show.plot = FALSE)
venn_down

plot.venn(venn_down)
# Note: Here A = down_sig_gene_TG & B = down_sig_gene_DRG
# Note: Out of 51 down_sig_gene_TG and 474 down_sig_gene_DRG, there are 37 common genes. 
# Note: Common genes
# # "ENSRNOG00000002475" "ENSRNOG00000009614" "ENSRNOG00000010144" "ENSRNOG00000011955" "ENSRNOG00000013141" "ENSRNOG00000022572" "ENSRNOG00000025179"
# [8] "ENSRNOG00000026917" "ENSRNOG00000027574" "ENSRNOG00000027791" "ENSRNOG00000029443" "ENSRNOG00000030442" "ENSRNOG00000031045" "ENSRNOG00000031439"
# [15] "ENSRNOG00000031517" "ENSRNOG00000031579" "ENSRNOG00000032771" "ENSRNOG00000032777" "ENSRNOG00000032947" "ENSRNOG00000033100" "ENSRNOG00000033146"
# [22] "ENSRNOG00000034597" "ENSRNOG00000038106" "ENSRNOG00000046271" "ENSRNOG00000047276" "ENSRNOG00000048365" "ENSRNOG00000049229" "ENSRNOG00000050669"
# [29] "ENSRNOG00000051249" "ENSRNOG00000052174" "ENSRNOG00000055067" "ENSRNOG00000055727" "ENSRNOG00000058555" "ENSRNOG00000058589" "ENSRNOG00000058710"
# [36] "ENSRNOG00000061099" "ENSRNOG00000061960"


# making the MA-plot.
# All up_down_reg_genes (TG):
plotMA(sig_gene_TG, main = "MA plot of\n All up_down_reg genes (TG)")
dev.off()

topGene_TG <- rownames(sig_gene_TG)[which.min(sig_gene_TG$padj)]
topGene_TG

# All genes ma_plot genes TG: 
plotMA(sig_gene_TG,ylim=c(-5,5), main = "MA plot of\n All up_down_reg genes\n with min_padj (TG)")
with(sig_gene_TG[topGene_TG,], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene_TG, pos=2, col="dodgerblue")
})
dev.off()

# Getting the full information about the gene with max foldchange and min padj value. 
up_down_TG[which.max(up_down_TG$log2FoldChange),]
# Gene wit the max log2foldchange.
# baseMean log2FoldChange     lfcSE     stat       pvalue         padj             geneID
# ENSRNOG00000008137 207.9229          1.589 0.1893237 8.393031 4.737669e-17 1.693009e-15 ENSRNOG00000008137

up_down_TG[which.min(up_down_TG$padj),]
# Gene with the min p.adjust value. 
# baseMean log2FoldChange      lfcSE     stat        pvalue         padj             geneID
# ENSRNOG00000005811 441.8145       1.548018 0.07201954 21.49441 1.756064e-102 3.150204e-98 ENSRNOG00000005811


# Plotting raw p-values.
# All de genes.
h_sig_Tg <- hist(sig_gene_TG$pvalue, col="lightblue", main="Histogram of raw P-values\n of all DE genes TG", breaks=20, xlab="P-value")
dev.off()

# upregulated
h_up<- hist(up_sig_gene_TG$pvalue, col="lightblue", main="Histogram of raw P-values\n up-regulated genes", breaks=20, xlab="P-value")
dev.off()

# downregulated
h_down <- hist(down_sig_gene_TG$pvalue, col="lightblue", main="Histogram of raw P-values\n down-regulated genes", breaks=20, xlab="P-value")
dev.off()

# Plotting the volcano plots
# all de genes
volc_plot_all_de_TG <- plot(x=sig_gene_TG$log2FoldChange, y=-log10(sig_gene_TG$padj),
                     xlab="log2(Fold-Change)", ylab="-log10(adjusted P-value",
                     col=ifelse(sig_gene_TG$log2FoldChange < -1, "red", "black"), main="Volcano plot\n all de genes TG\n less than -1")

dev.off()


################################# DRG ##############################################

sample_table_DRG <- sample_table[grep("DRG",sample_table$sampleName),]

DEseq_data_DRG <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table_DRG, directory = HTSeq_direc, design = ~ breeder)
dim(DEseq_data_DRG)
# [1] 32662    11

# Note:
# Here there are 11 samples as one sample from Envigo "115" have
# very low RIN value. 
#[1] 32662    11

# remove genes that have 0 to 1 count for all samples.
filt_DEseq_data_DRG <- DEseq_data_DRG[ rowSums(counts(DEseq_data_DRG))>1,]
dim(filt_DEseq_data_DRG)
#[1] 21699    11

# Extracting the count dataset. 
DRG_cou <- as.data.frame(counts(Norm_filt_DEseq_data_DRG, normalized = TRUE))
class(filt_norm_DRG_cou) # [1] "data.frame"
dim(filt_norm_DRG_cou)
#[1] 21699    11

#Caculate the size_factor, normalization
Norm_filt_DEseq_data_DRG <- estimateSizeFactors(filt_DEseq_data_DRG)
dim(Norm_filt_DEseq_data_DRG)
#[1] 21699    11

# Extracting counts per sample. 
DRG_cou_Taco <- DRG_cou[,1:6]
DRG_cou_Envi <- DRG_cou[,7:11]


# with squared CV (DRG)
#Taco
cv_DRG_taco_cou_sq_cv <- apply(DRG_cou_Taco ,1, function(x) ((sd(x)/mean(x))^2))
cv_melted_data_DRG_taco_sq <- melt(cv_DRG_taco_cou_sq_cv)

#¤¤ plotting squared cv
# DRG_taco
p_cv_DRG_taco_sq <-  ggplot(cv_melted_data_DRG_taco_sq ,aes(x = value))+ geom_vline(xintercept = 1, color = "blue") + geom_density()+
  annotate("text", x = 7, y = 6, label = ("cutoff"), color="blue")+ xlab("Coefficient of Variation") +
  ylab("Density")+ ggtitle("Density-plot of Cofficient of Variation\n Naive Rats Taconic (DRG)")+ theme_bw() 

p_cv_DRG_taco_sq
dev.off()

#Envi
cv_DRG_envi_cou_sq_cv <- apply(DRG_cou_Envi ,1, function(x) ((sd(x)/mean(x))^2))
cv_melted_data_DRG_envi_sq <- melt(cv_DRG_envi_cou_sq_cv)

# TG_envi
p_cv_DRG_envi_sq <- ggplot(cv_melted_data_DRG_envi_sq ,aes(x = value))+ geom_vline(xintercept = 1, color = "blue") + geom_density()+
  annotate("text", x = 7, y = 6, label = ("cutoff"), color="blue")+ xlab("Coefficient of Variation") +
  ylab("Density")+ ggtitle("Density-plot of Cofficient of Variation\n Naive Rats Envigo (DRG)")+ theme_bw() 

p_cv_DRG_envi_sq
dev.off()

#"""""""""""" Calcualting coffeficient of variation for (DRG) tissue. before DE testing """"""""""""""
# Taco
cv_DRG_taco_cou <- apply(DRG_cou_Taco ,1, function(x) ((sd(x)/mean(x))*100))
cv_melted_data_DRG_taco <- melt(cv_DRG_taco_cou)

# plotting cv without squared
# TG_taco
p_cv_DRG_taco <-  ggplot(cv_melted_data_DRG_taco ,aes(x = value)) + geom_density()+
  xlab("Coefficient of Variation") +
  ylab("Density")+ ggtitle("Density-plot of Cofficient of Variation\n Naive Rats Taconic (DRG)")+ theme_bw() 

p_cv_DRG_taco
dev.off()
# Envi
cv_DRG_envi_cou <- apply(DRG_cou_Envi ,1, function(x) ((sd(x)/mean(x))*100))
cv_melted_data_DRG_envi <- melt(cv_DRG_envi_cou)

# TG_envi
p_cv_DRG_envi <- ggplot(cv_melted_data_DRG_envi ,aes(x = value)) + geom_density()+
  xlab("Coefficient of Variation") +
  ylab("Density")+ ggtitle("Density-plot of Cofficient of Variation\n Naive Rats Envigo (DRG)")+ theme_bw() 

p_cv_DRG_envi
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@ calc means (TG) , before de testing @@@@@@@@@@@@@@@@@@@
# Taco
calc_means_DRG_taco <- rowMeans(DRG_cou_Taco, na.rm = TRUE, dims = 1)
calc_log2_means_DRG_taco <- log2(calc_means_DRG_taco)

# plotting calc means
# For taco
qplot(calc_log2_means_DRG_taco , cv_DRG_taco_cou, xlab = "log2(means)", ylab = "CV", main = "CV versus log2(mean) plot of Taconic (DRG)") + geom_smooth() 
dev.off()

# envigo
calc_means_DRG_envi <- rowMeans(DRG_cou_Envi, na.rm = TRUE, dims = 1)
calc_log2_means_DRG_envi <- log2(calc_means_DRG_envi)


# For envigo
qplot(calc_log2_means_DRG_envi  , cv_DRG_envi_cou , xlab = "log2(means)", ylab = "CV", main = "CV versus log2(mean) plot of Envigo (DRG)") + geom_smooth()
dev.off()

# Compute the density (TG) .
# taco
dens_taco_DRG <- density(calc_log2_means_DRG_taco)

# plotting the density of means
# For Taconic
plot(dens_taco_DRG, type="l", col="green",
     xlab="x", main="Density plot of the log2(mean) of Taconic (DRG)",lwd=2)

# envigo
dens_envi_DRG <- density(calc_log2_means_DRG_envi)

# For Envigo
plot(dens_envi_DRG, type="l", col="red",
     xlab="x", main="Density plot of the log2(mean) of Envigo (DRG)",lwd=2)

# Combining both suppliers.	
plot(dens_envi_DRG, col= "red", lwd = 3, main="Density plot of the log2(mean) of Envigo & Taconic (DRG)") 
lines(dens_taco_DRG, col = "blue", lwd = 3)
legend('topleft', c("taconic", "envigo"), lwd=c(2.5,2.5),col=c("blue", "red"), cex = 0.75)

# computing sample distances to see how similar are the samples to each other. 
rld <- rlog(Norm_filt_DEseq_data_DRG, blind=FALSE)
head(assay(rld), 3)

# Making the PCA plots
data <- plotPCA(rld, intgroup= "RIN" , returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color = RIN ,label = name)) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
dev.off()

#Estimating dispersions. 
disper_DEseq_data_DRG <- estimateDispersions(Norm_filt_DEseq_data_DRG)
dim(disper_DEseq_data_DRG)
#[1] 21699    11

# plot dispersion ('This is the dispersion plot for DRG both breeds')
plotDispEsts(disper_DEseq_data_DRG)
dev.off()

# doing the negative binomial wald test
#DRG
nbinom_DEseq_data_DRG <- nbinomWaldTest(disper_DEseq_data_DRG)
dim(nbinom_DEseq_data_DRG)
#[1] 21699    11

# Adjusting the p-value to 0.05
#DRG
res_DRG_new_Pval <- results(nbinom_DEseq_data_DRG, alpha = 0.05)
res_DRG_new_Pval
# log2 fold change (MAP): breeder Taconic vs Envigo 
# Wald test p-value: breeder Taconic vs Envigo 
# DataFrame with 21699 rows and 6 columns

dim(res_DRG_new_Pval)
#[1] 21699    6
summary(res_DRG_new_Pval)

# checking the number of DE-genes
sum(res_DRG_new_Pval$padj < 0.05, na.rm=TRUE)
# Total of 6568 De genes. 

sig_gene_DRG <- subset(res_DRG_new_Pval, res_DRG_new_Pval$padj < 0.05)
sig_gene_DRG$geneID <- rownames(sig_gene_DRG)

dim(sig_gene_DRG)
#[1] 6568    7
summary(sig_gene_DRG)

# out of 6568 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 3165, 48% 
# LFC < 0 (down)   : 3403, 52% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 0, 0% 
# (mean count < 1)

# this is to get the information about  result section of deseq. 
mcols(sig_gene_DRG, use.names = TRUE)

# Helps to find out which gene is higher expressed in the breeders. 
topGene <- rownames(sig_gene_DRG)[which.min(sig_gene_DRG$padj)]
plotCounts(nbinom_DEseq_data_DRG, gene = topGene, intgroup="breeder", returnData = FALSE)

#################### Merging DE genes with original DB 21699 genes. ################

# extracting the count genes
filt_norm_DRG_cou <- as.data.frame(counts(Norm_filt_DEseq_data_DRG, normalized = TRUE))
filt_norm_DRG_cou$geneID <- rownames(filt_norm_DRG_cou)


# Extracting counts per sample. 
# Taco
all_DRG_genes_Taco <- filt_norm_DRG_cou[,1:6]
all_DRG_genes_Taco$geneID <- filt_norm_DRG_cou$geneID
all_DRG_genes_Taco$log2foldchange <- res_DRG_new_Pval$log2FoldChange
all_DRG_genes_Taco$padj <- res_DRG_new_Pval$padj

# Calculating CV after de testing.
# Taco
DE_DRG_taco_cou <- apply(all_DRG_genes_Taco[,1:6] ,1, function(x) ((sd(x)/mean(x))*100))
DE_DRG_taco_cou_melted <- melt(DE_DRG_taco_cou)
all_DRG_genes_Taco$CV <- DE_DRG_taco_cou_melted$value

# Calculating Means after de testing.
# Taco
DE_means_DRG_taco <- rowMeans(all_DRG_genes_Taco[,1:6], na.rm = TRUE, dims = 1)
DE_log2_means_DRG_taco <- log2(DE_means_DRG_taco)
DE_log2_means_DRG_taco <- as.data.frame(DE_log2_means_DRG_taco)
all_DRG_genes_Taco$means <- DE_log2_means_DRG_taco$DE_log2_means_DRG_taco

# Extracting counts per sample.
# Envi
all_DRG_genes_Envigo <- filt_norm_DRG_cou[,7:12]
all_DRG_genes_Envigo$geneID <- filt_norm_DRG_cou$geneID
all_DRG_genes_Envigo$log2foldchange <- res_DRG_new_Pval$log2FoldChange
all_DRG_genes_Envigo$padj <- res_DRG_new_Pval$padj

# Calculating CV after de testing.
# Envi
DE_DRG_envi_cou <- apply(all_DRG_genes_Envigo[,1:5] ,1, function(x) ((sd(x)/mean(x))*100))
DE_melted_data_DRG_envi <- melt(DE_DRG_envi_cou)
all_DRG_genes_Envigo$CV <- DE_melted_data_DRG_envi$value

# Calculating Means after de testing.
# envigo
DE_means_DRG_envi <- rowMeans(all_DRG_genes_Envigo[,1:5], na.rm = TRUE, dims = 1)
DE_log2_means_DRG_envi <- log2(DE_means_DRG_envi)
DE_log2_means_DRG_envi <- as.data.frame(DE_log2_means_DRG_envi)
all_DRG_genes_Envigo$means <- DE_log2_means_DRG_envi$DE_log2_means_DRG_envi


res_DRG_taco_envi <- as.data.frame(res_DRG_new_Pval)
res_DRG_taco_envi$gene_id <- rownames(res_DRG_taco_envi)
res_DRG_taco_envi$cv_taco <- all_DRG_genes_Taco$CV
res_DRG_taco_envi$mean_taco <- all_DRG_genes_Taco$means
res_DRG_taco_envi$cv_envi <- all_DRG_genes_Envigo$CV   
res_DRG_taco_envi$mean_envi <- all_DRG_genes_Envigo$means
res_DRG_taco_envi <- subset(res_DRG_taco_envi, res_DRG_taco_envi$padj < 0.05)
res_DRG_taco_envi <- res_DRG_taco_envi[(res_DRG_taco_envi$cv_taco < 75 & res_DRG_taco_envi$cv_envi < 75), ]
res_DRG_taco_envi <- res_DRG_taco_envi[which(res_DRG_taco_envi$mean_taco >5 | res_DRG_taco_envi$mean_envi >5), ]
res_DRG_taco_envi_up <- subset(res_DRG_taco_envi, res_DRG_taco_envi$log2FoldChange > 1)
res_DRG_taco_envi_down <- subset(res_DRG_taco_envi, res_DRG_taco_envi$log2FoldChange < -1)

dim(res_DRG_taco_envi)
dim(res_DRG_taco_envi_up)
dim(res_DRG_taco_envi_down)


#Identfying most significant up-reg genes:
sig_reg_drg_taco_envi_up <- as.data.frame(res_DRG_taco_envi_up$gene_id)
View(sig_reg_drg_taco_envi_up)

write.table(sig_reg_drg_taco_envi_up, file = "sig_reg_drg_up_taco-envi.txt", sep = "\t", col.names = FALSE , row.names = FALSE, quote = FALSE)

sig_reg_drg_names_up <- read.table("sig_drg_taco_envi_up.txt", header = TRUE, sep = ",")
colnames(sig_reg_drg_names_up)[1] <- "gene_id"
View(sig_reg_drg_names_up)

sig_reg_drg_up <- merge.data.frame(res_DRG_taco_envi_up, sig_reg_drg_names_up, by = "gene_id", all = TRUE)
dim.data.frame(sig_reg_drg_up)
View(sig_reg_drg_up)

sorted_sig_reg_drg_up <- sig_reg_drg_up[order(-sig_reg_drg_up$log2FoldChange),c(1,3,7,9,11,12,13,14)]
View(sorted_sig_reg_drg_up)
write.table(sorted_sig_reg_drg_up, file = "sig_drg_up_info.txt", sep = "\t", col.names = TRUE , row.names = FALSE, quote = FALSE)


#Identfying most significant down-reg genes:
sig_reg_drg_taco_envi_down <- as.data.frame(res_DRG_taco_envi_down$gene_id)
View(sig_reg_drg_taco_envi_down)

write.table(sig_reg_drg_taco_envi_down, file = "sig_reg_drg_down_taco-envi.txt", sep = "\t", col.names = FALSE , row.names = FALSE, quote = FALSE)

sig_reg_drg_names_down <- read.table("sig_drg_taco_envi_down.txt", header = TRUE, sep = ",")
colnames(sig_reg_drg_names_down)[1] <- "gene_id"

sig_reg_drg_down <- merge.data.frame(res_DRG_taco_envi_down, sig_reg_drg_names_down, by = "gene_id", all = TRUE)
dim.data.frame(sig_reg_drg_down)
View(sig_reg_drg_down)

sorted_sig_reg_drg_down <- sig_reg_drg_down[order(sig_reg_drg_down$log2FoldChange),c(1,3,7,9,11,12,13,14)]
View(sorted_sig_reg_drg_down)
write.table(sorted_sig_reg_drg_down, file = "sig_drg_down_info.txt", sep = "\t", col.names = TRUE , row.names = FALSE, quote = FALSE)



# volcano plot
#sig=ifelse(all_TG_genes_Taco$padj<0.05, "FDR<0.05","no_sig")
qplot(res_DRG_taco_envi$log2FoldChange  , res_DRG_taco_envi$cv_envi, xlab = "log2_fold_change", ylab = "CV", main = "CV versus log2(foldchange) plot of Envigo\n all DE (DRG)" )+ 
  annotate("text", x = 1, y = 200, label = paste0("",dim(res_DRG_taco_envi)[1]," DE genes"), color="red")
dev.off()

# Mean and cv distribution plot. 
qplot(res_DRG_taco_envi$mean_taco , res_DRG_taco_envi$log2FoldChange, xlab = "log2(means)", ylab = "CV", main = "CV versus log2(mean) plot of Taconic\n genes (DRG)") + geom_smooth() + 
  annotate("text", x = 16, y = 200, label = paste0("",dim(res_DRG_taco_envi)[1]," genes"), color="red")



# making the MA-plot.
# upregulated
plotMA(up_sig_gene_DRG, main="MA plot of\n up_regulated genes (DRG)", ylim=c(-2,2))

# down-regulated
plotMA(down_sig_gene_DRG, main="MA plot of\n down_regulated genes (DRG)")

idx <- identify(sig_gene_TG$baseMean)
rownames(sig_gene_TG)[idx]

# Plotting raw p-values. 
h_sig_DRG <- hist(m_DRG_genes$pvalue, col="lightblue", main="Histogram of raw P-values\n of all DE genes (DRG)", breaks=20, xlab="P-value")
dev.off()


# upregulated
h_up<- hist(up_sig_gene_DRG$pvalue, col="lightblue", main="Histogram of raw P-values\n up-regulated genes (DRG)", breaks=20, xlab="P-value")
dev.off()

# downregulated
h_down <- hist(down_sig_gene_DRG$pvalue, col="lightblue", main="Histogram of raw P-values\n down-regulated genes (DRG)", breaks=20, xlab="P-value")
dev.off()

# Plotting the volcano plots
# upregulated
volc_plot_up <- plot(x=up_sig_gene_DRG$log2FoldChange, y=-log10(up_sig_gene_DRG$padj),
                     xlab="log2(Fold-Change)", ylab="-log10(adjusted P-value",
                     col=ifelse(up_sig_gene_DRG$padj<=0.05, "red", "black"), main="Volcano plot\n upregulated-genes (DRG)")
dev.off()

# downregulated
volc_plot_down <- plot(x=down_sig_gene_DRG$log2FoldChange, y=-log10(down_sig_gene_DRG$padj),
                       xlab="log2(Fold-Change)", ylab="-log10(adjusted P-value",
                       col=ifelse(down_sig_gene_DRG$padj<=0.05, "red", "black"), main="Volcano plot\n downregulated-genes (DRG)")
dev.off()
