# Go-seq analysis and KEGG pathway analysis. 

# location of library in work-pc
library("goseq", lib.loc = "C:/Program Files/R/R-3.2.3/library")
library("GO.db", lib.loc="C:/Program Files/R/R-3.2.3/library")
library("org.Rn.eg.db", lib.loc="C:/Program Files/R/R-3.2.3/library")
library("KEGGgraph", lib.loc="C:/Program Files/R/R-3.2.3/library")
library("pathview", lib.loc="C:/Program Files/R/R-3.2.3/library")
library("biomaRt", lib.loc="C:/Program Files/R/R-3.2.3/library")

# The code is seprated accordingly by tissues TG and DRG. There are examples of what the result look like. 

################################# TG ##############################################

# Get the path of working directory. 
getwd()

# set the working directory for computerome. Please change the directory.
#setwd("/home/people/prasing/Master_Thesis/scripts/")

# Note: It is important to set the working directory same as where the files are generated by DEseq.

# loading the table Naive rats
# upregulated
up_sig_TG <- read.table("UPREG_filt_cv-mean_tg.txt", header = TRUE)
dim(up_sig_TG)
#View(up_sig_TG)
# Genes filtered with cv and means
# [1] 85 11


# downregulated
down_sig_TG <- read.table("DOWNREG_filt_cv-mean_tg.txt", header = TRUE)
dim(down_sig_TG)
#View(down_sig_TG)
# # Genes filtered with cv and means
# [1]  12 11

# Merge up and down-reg DBs. 
# Note: This will help to get more power. 
m_TG <- rbind.data.frame(up_sig_TG,down_sig_TG)
dim(m_TG)
# # Genes filtered with cv and means
# 97  11


# loading the dataset consisting of all the genes investigated for DE Analysis. 
load("dds_TG_Naive.R")

# extracting the genes that were investigated for DE-analysis. 
assayed_genes_TG <- rownames(nbinom_DEseq_data_TG)
length(assayed_genes_TG)
# 22102 genes (this is a character vector)


# All DE_Genes
de_genes_TG_all <- rownames(m_TG)
length(de_genes_TG_all)
# Genes filtered with cv and means
# 97


# extract upregulated DE-genes 
de_genes_TG_up <- rownames(up_sig_TG)
length(de_genes_TG_up)
# There are 85 upregulated. 


# exrtact downregulated DE-genes:
de_genes_TG_down <- rownames(down_sig_TG)
length(de_genes_TG_down)
# There are 12 downregulated.


#All
go_vector_TG_all <- as.integer(assayed_genes_TG %in% de_genes_TG_all)
names(go_vector_TG_all)=assayed_genes_TG
table(go_vector_TG_all)
# 0     1 
# 22005    97 

# upregulated
go_vector_TG_up <- as.integer(assayed_genes_TG %in% de_genes_TG_up)
names(go_vector_TG_up)=assayed_genes_TG
table(go_vector_TG_up)
# 0     1 
# 22017    85

# down-regulated
go_vector_TG_down <- as.integer(assayed_genes_TG %in% de_genes_TG_down)
names(go_vector_TG_down)=assayed_genes_TG
table(go_vector_TG_down)
# 0     1 
# 22090    12 

# All
pwf_TG_all <- nullp(go_vector_TG_all, "rn5", "ensGene")
dim(pwf_TG_all)
#[1] 22102     3

# upregulated
pwf_TG_up <- nullp(go_vector_TG_up, "rn5", "ensGene")
dim(pwf_TG_up)
#[1] 22102     3 


# down-regulated 
pwf_TG_down <- nullp(go_vector_TG_down, "rn5", "ensGene")
dim(pwf_TG_down)
#[1] 22102     3


# All
GO_TG_all <- goseq(pwf_TG_all, "rn5", "ensGene")
dim(GO_TG_all)
#[1] 20000     7


# upregulated
GO_TG_up <- goseq(pwf_TG_up, "rn5", "ensGene")
dim(GO_TG_up)
#[1] 20000     7


# down-regulated 
GO_TG_down <- goseq(pwf_TG_down, "rn5", "ensGene")
dim(GO_TG_down)
#[1] 20000     7


# All
enriched_GO_TG_all <- GO_TG_all$category[p.adjust(GO_TG_all$over_represented_pvalue, method = "BH")<0.05]
head(enriched_GO_TG_all)
length(enriched_GO_TG_all)
# 0 enriched GO-terms


# upregulated
enriched_GO_TG_up <- GO_TG_up$category[p.adjust(GO_TG_up$over_represented_pvalue, method = "BH")<0.05]
head(enriched_GO_TG_up)
length(enriched_GO_TG_up)
# @@@@ Up-regulated genes @@@@
# 0 enriched GO-terms



# down-regulated
enriched_GO_TG_down <- GO_TG_down$category[p.adjust(GO_TG_down$over_represented_pvalue, method = "BH")<0.05]
head(enriched_GO_TG_down)
length(enriched_GO_TG_down)
# @@@@ down-regulated genes @@@@
# 0 enriched GO-terms


#There were no over represented genes when filtered for both CV and means


# Printing the GO-terms: 
# TG
for (go in enriched_GO_TG_all[1:14]){  
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}


for (go in enriched_GO_TG_up[1:2]){  
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}


for (go in enriched_GO_TG_down[1:4]) {
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}


#KEGG Pathway Analysis. Currently used method for KEGG
# Get the mapping from ENSEMBL 2 Entrez
en2eg=as.list(org.Rn.egENSEMBL2EG)
# Get the mapping from Entrez 2 KEGG
eg2kegg <- as.list(org.Rn.egPATH)
# Define a function which gets all unique KEGG IDs, associated with a set of Entrez IDs.
grepKEGG=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
# Apply this function to every entry in the mapping from, ENSEMBL 2 Entrez to combine the two maps
kegg=lapply(en2eg,grepKEGG,eg2kegg)

class(kegg)
# list
head(kegg, n=10L)
#Length
length(kegg)
#[1] 19858 


#TG
#All
KEGG_all <- goseq(pwf_TG_all,gene2cat = kegg)
# Note: For 17075 genes, we could not find any categories. These genes will be excluded. This was after foldchange more than 1.
dim(KEGG_all)
#[1] 225   5


# upregulated
KEGG_up <- goseq(pwf_TG_up,gene2cat = kegg)
# Note: For 17075 genes, we could not find any categories. These genes will be excluded. This was after foldchange more than 1.
dim(KEGG_up)
#[1] 225   5


# down-regulated
KEGG_down <- goseq(pwf_TG_down, gene2cat = kegg)
# Note: For 17075 genes, we could not find any categories. These genes will be excluded.
dim(KEGG_down)
#[1] 225   5


# All
enriched_KEGG_all <-  KEGG_all$category[p.adjust(KEGG_all$over_represented_pvalue, method = "BH")<0.05]
enriched_KEGG_all
# [1] "04970"


# upregulated
enriched_KEGG_up<-  KEGG_up$category[p.adjust(KEGG_up$over_represented_pvalue, method = "BH")<0.05]
length(enriched_KEGG_up) 
# 0 enriched KEGG categories
# There were no over represented pathways when filtered for both CV and means also indivisually 
# for Cv and means. 


# down-regulated
enriched_KEGG_down <- KEGG_down$category[p.adjust(KEGG_down$over_represented_pvalue, method = "BH")<0.05]
length(enriched_KEGG_down)
# 0 enriched KEGG categories
# There were no over represented pathways when filtered for both CV and means also indivisually 
# for Cv and means. 


# The above code is used for getting gene names and gene symbols from biomart which are further used for making kegg pathways using 
# online pathview application. https://pathview.uncc.edu/home
# All genes TG
all_genes_TG <- read.table("all_genes_without_filtering_study_1.txt", header = TRUE, sep = "\t")
all_genes_TG$geneID <- rownames(all_genes_TG)
dim(all_genes_TG)
View(all_genes_TG)

all_de_genes_TG <- as.data.frame(all_genes_TG$geneID)
View(all_de_genes_TG)

write.table(all_de_genes_TG, file = "all_de_gene_ids_TG.txt", sep = "\t", col.names = FALSE , row.names = FALSE, quote = FALSE)

all_de_names_TG <- read.table("nofilt_de_genes_tg.txt", header = TRUE, sep = ",")
colnames(all_de_names_TG)[1] <- "geneID"
View(all_de_names_TG)

m_TG_de_all <- merge.data.frame(all_de_names_TG, all_genes_TG, by = "geneID")
dim.data.frame(m_TG_de_all)
View(m_TG_de_all)

path_names_TG_de_all <- as.data.frame(m_TG_de_all[,c(2,4)])
dim(path_names_TG_de_all)

write.table(path_names_TG_de_all, file = "Pathview_gene_names_TG_all_de.txt", sep = "\t", col.names = TRUE , row.names = FALSE , quote = FALSE)



################################# DRG ##############################################

# loading the table Naive rats
# upregulated
up_sig_DRG <- read.table("UPREG_cv-mean-filt_drg.txt", header = TRUE)
dim(up_sig_DRG)

# downregulated
down_sig_DRG <- read.table("DOWNREG_cv-mean-filt_drg.txt", header = TRUE)
dim(down_sig_DRG)

# All
m_DRG_all <- rbind.data.frame(down_sig_DRG, up_sig_DRG)
dim(m_DRG_all)
# [1] 638   9
# After filtering DE genes by Cv and mean:
# [1] 385  11


# loading the dataset consisting of all the genes investigated for DE Analysis. 
load("dds_DRG_Naive.R")

# extracting the genes that were investigated for DE-analysis. 
assayed_genes_DRG <- rownames(nbinom_DEseq_data_DRG)
length(assayed_genes_DRG)
# 21699 genes (this is a character vector)


# extract all DE-genes
de_genes_DRG_all <- rownames(m_DRG_all)
length(de_genes_DRG_all)
# There are 638 DE genes from the whole list of DE genes. 
# After filtering DE genes by Cv and mean:
# 385

# extract upregulated DE-genes 
de_genes_DRG_up <- rownames(up_sig_DRG)
length(de_genes_DRG_up)
# There are 140 upregulated from whole list of DE genes. 
# There are 161 upregulated from whole list of filtered CV DE genes. 
# After filtering DE genes by Cv and mean:
# [1] 146

# exrtact downregulated DE-genes:
de_genes_DRG_down <- rownames(down_sig_DRG)
class(de_genes_DRG_down)
# [1] "character"
length(de_genes_DRG_down)
# [1] 173
# There are 173 downregulated genes from the whole list of DE genes. 
# There are 418 downregulated genes from the whole list of filtered CV DE genes. 
# After filtering DE genes by Cv and mean:
# [1] 239

#Constructing the vector for GOseq:
# All genes
go_vector_DRG_all <- as.integer(assayed_genes_DRG %in% de_genes_DRG_all)
names(go_vector_DRG_all)=assayed_genes_DRG
table(go_vector_DRG_all)
# 0     1 
# 21061   638 
# After filtering DE genes by Cv and mean:
#  0     1 
# 21314   385


# upregulated
go_vector_DRG_up <- as.integer(assayed_genes_DRG %in% de_genes_DRG_up)
names(go_vector_DRG_up)=assayed_genes_DRG
table(go_vector_DRG_up)
# @@@@  After log2foldchange > 1 filtering @@@@
# 0     1 
# 21559   140 
#After filtering DE genes based on CV we get:
# 21538   161
# After filtering DE genes by Cv and mean:
#  0     1 
# 21553   146 


# down-regulated
go_vector_DRG_down <- as.integer(assayed_genes_DRG %in% de_genes_DRG_down)
names(go_vector_DRG_down)=assayed_genes_DRG
table(go_vector_DRG_down)
# @@@@ After log2foldchange > 1 filtering @@@@
# 0     1 
# 21526   173 
#After filtering DE genes based on CV we get:
#21281   418
# After filtering DE genes by Cv and mean:
#  0     1 
# 21460   239 


# calculate pwf(probablity weight function):
# All genes
pwf_DRG_all <- nullp(go_vector_DRG_all, "rn5", "ensGene")
dim(pwf_DRG_all)
#[1] 21699     3


# upregulated
pwf_DRG_up <- nullp(go_vector_DRG_up, "rn5", "ensGene")
class(pwf_DRG_up)
#[1] "data.frame"
dim(pwf_DRG_up)
#[1] 21699     3
#View(pwf_DRG_up)



# down-regulated 
pwf_DRG_down <- nullp(go_vector_DRG_down, "rn5", "ensGene")
#class(pwf_DRG_down)
#[1] "data.frame"
dim(pwf_DRG_down)
#[1] 21699     3
#View(pwf_DRG_down)


# GO category over representation in DE genes,using the Wallenius approximation.
# All genes
GO_DRG_all <- goseq(pwf_DRG_all, "rn5", "ensGene")
dim(GO_DRG_all)
#[1] 19961     7
# Note: For 6707 genes, we could not find any categories. These genes will be excluded.

# upregulated
GO_DRG_up <- goseq(pwf_DRG_up, "rn5", "ensGene")
dim(GO_DRG_up)
#[1] 19961     7
class(GO_DRG_up)
#[1] "data.frame"
# Note: For 6707 genes, we could not find any categories. These genes will be excluded.

# down-regulated
GO_DRG_down <- goseq(pwf_DRG_down, "rn5", "ensGene")
dim(GO_DRG_down)
#[1] 19961     7
class(GO_DRG_down)
# [1] "data.frame"
# Note: For 6707 genes, we could not find any categories. These genes will be excluded.


# Extract GO categories overenriched. Using FDR-corrected p-values. 
# All genes
enriched_GO_DRG_all <- GO_DRG_all$category[p.adjust(GO_DRG_all$over_represented_pvalue, method = "BH")<0.05]

head(enriched_GO_DRG_all)
#  "GO:0002376" "GO:0006955" "GO:0005833" "GO:0002684" "GO:0002682" "GO:0015669"

length(enriched_GO_DRG_all)
# 69 GO-terms in total
# After filtering DE genes by Cv and mean:
# 23 GO-terms in total

# ordering the dataset in ascending order. 
newdata_all <- GO_DRG_all[order(GO_DRG_all$over_represented_pvalue),] 
View(newdata_all)

# extracting the particular columns. 
temp_tab_all <- newdata_all[c(1:10),c(1,2,6,7,4,5)] 
View(temp_tab_all)

write.table(temp_tab_all, file = "GO-terms-combined-DRG.txt", sep = "\t", col.names = TRUE , row.names = FALSE, quote = FALSE)



# upregulated
enriched_GO_DRG_up<-GO_DRG_up$category[p.adjust(GO_DRG_up$over_represented_pvalue, method = "BH")<0.05]

class(enriched_GO_DRG_up)
# [1] "character"

head(enriched_GO_DRG_up)
#[1] "GO:0005886" "GO:0071944" "GO:0097285" "GO:0007610" "GO:0005887" "GO:0048699"

length(enriched_GO_DRG_up)
# @@@@ Up-regulated genes @@@@
# 74 enriched GO-terms
# After filtering DE genes based on CV we get:
# 60 enriched GO-terms
# After filtering DE genes by Cv and mean:
# 42 GO-terms in total

# ordering the dataset in ascending order. 
newdata_up <- GO_DRG_up[order(GO_DRG_up$over_represented_pvalue),] 
View(newdata_up)

# extracting the particular columns. 
temp_tab_up <- newdata_up[c(1:10),c(1,2,6,7,4,5)] 
View(temp_tab_up)

write.table(temp_tab_up, file = "GO-terms-upreg-DRG.txt", sep = "\t", col.names = TRUE , row.names = FALSE, quote = FALSE)


# down-regulated
enriched_GO_DRG_down <- GO_DRG_down$category[p.adjust(GO_DRG_down$over_represented_pvalue, method = "BH")<0.05]

head(enriched_GO_DRG_down)
#[1] "GO:0002376" "GO:0006955" "GO:0050776" "GO:0002682" "GO:0002684" "GO:0005833"

length(enriched_GO_DRG_down)
# @@@@ down-regulated genes @@@@
# 80 enriched GO-terms
# After filtering DE genes based on CV we get:
# 71 enriched GO-terms
# After filtering DE genes by Cv and mean:
# 48 GO-terms in total

# ordering the dataset in ascending order. 
newdata_down <- GO_DRG_down[order(GO_DRG_down$over_represented_pvalue),] 
View(newdata_down)

# extracting the particular columns. 
temp_tab_down <- newdata_down[c(1:10),c(1,2,6,7,4,5)] 
View(temp_tab_down)

write.table(temp_tab_down, file = "GO-terms-downreg-DRG.txt", sep = "\t", col.names = TRUE , row.names = FALSE, quote = FALSE)



# Printing the GO-terms: 
# TG
for (go in enriched_GO_DRG_up[1:10]) {
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}


for (go in enriched_GO_DRG_down[1:43]) {
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}

for (go in enriched_GO_DRG_all[1:13]) {
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}


#KEGG Pathway Analysis. Currently used method for KEGG
# Get the mapping from ENSEMBL 2 Entrez
en2eg=as.list(org.Rn.egENSEMBL2EG)
# Get the mapping from Entrez 2 KEGG
eg2kegg <- as.list(org.Rn.egPATH)
# Define a function which gets all unique KEGG IDs, associated with a set of Entrez IDs.
grepKEGG=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
# Apply this function to every entry in the mapping from, ENSEMBL 2 Entrez to combine the two maps
kegg=lapply(en2eg,grepKEGG,eg2kegg)

class(kegg)
# list
head(kegg, n=10L)
#Length
length(kegg)
#[1] 19858 


###########################################################
# All de_genes.
all_de_genes_DRG <- read.table("all_genes_without_filtering_study_1_DRG.txt", header = TRUE, sep = "\t")
all_de_genes_DRG$geneID <- rownames(all_de_genes_DRG)
class(all_de_genes_DRG)
dim(all_de_genes_DRG)
View(all_de_genes_DRG)

all_de_genes_id_DRG <- as.data.frame(all_de_genes_DRG$geneID)
View(all_de_genes_id_DRG)

write.table(all_de_genes_id_DRG, file = "all_de_gene_ids_DRG.txt", sep = "\t", col.names = FALSE , row.names = FALSE, quote = FALSE)

all_de_names_DRG <- read.table("all_de_gene_names_drg.txt", header = TRUE, sep = ",")
colnames(all_de_names_DRG)[1] <- "geneID"
View(all_de_names_DRG)

m_DRG_de_all <- merge.data.frame(all_de_names_DRG, all_de_genes_DRG, by = "geneID")
dim.data.frame(m_DRG_de_all)
View(m_DRG_de_all)

path_names_DRG_de_all <- as.data.frame(m_DRG_de_all[,c(2,4)])
dim(path_names_DRG_de_all)
View(path_names_DRG_de_all)

write.table(path_names_DRG_de_all, file = "Pathview_gene_names_DRG_nofilt_de.txt", sep = "\t", col.names = TRUE , row.names = FALSE , quote = FALSE)
