# Constructing the sample_table 
setwd("/home/jamshed/Desktop/Datasets/")

# Directory with HTSeq output
directory <- "/home/jamshed/Desktop/Datasets/study1_naive_diff/"

# extract files from directory
sample_files <- list.files(directory)

# load meta-data
meta_data <- read.csv("/home/jamshed/Desktop/Datasets/metadata_navie.csv")


# making the new-data frame with sample info. 
sample_table <- data.frame(sampleName = gsub(".txt", "", list.files(directory)), fileName = list.files(directory), stringsAsFactors=FALSE)


# library("compare")
# 
# compare(meta_data$HTSeqfile,sample_table$fileName)
# setdiff(meta_data$HTSeqfile,sample_table$fileName)

# construct the sample table
sample_table$id <- meta_data$ID[match(sample_table$sampleName, meta_data$Sample)]
sample_table$RIN <- meta_data$RIN[match(sample_table$sampleName, meta_data$Sample)]
sample_table$tissue <- meta_data$Tissue[match(sample_table$sampleName, meta_data$Sample)]
sample_table$breeder <- meta_data$Breeder[match(sample_table$sampleName, meta_data$Sample)]
sample_table$breed <- meta_data$Breed[match(sample_table$sampleName, meta_data$Sample)]
sample_table$pretreatment <- meta_data$Pretreatment[match(sample_table$sampleName, meta_data$Sample)]
sample_table$treatment <- meta_data$Treatment[match(sample_table$sampleName, meta_data$Sample)]
sample_table$age <- meta_data$Age_weeks[match(sample_table$sampleName, meta_data$Sample)]
sample_table$sex <- meta_data$Sex[match(sample_table$sampleName, meta_data$Sample)]
sample_table$weight_operation <- meta_data$Weight_operation[match(sample_table$sampleName, meta_data$Sample)]
sample_table$weight_treatment <- meta_data$weight_treatment[match(sample_table$sampleName, meta_data$Sample)]
sample_table$weight_loss <- meta_data$Weight_loss[match(sample_table$sampleName, meta_data$Sample)]
sample_table$treatment_date <- meta_data$treatment_date[match(sample_table$sampleName, meta_data$Sample)]
sample_table$operation_date <- meta_data$Operation_date[match(sample_table$sampleName, meta_data$Sample)]
sample_table$day_arr_treatment <- meta_data$Days_arr_operation[match(sample_table$sampleName, meta_data$Sample)]
sample_table$days_operation_treatment <- meta_data$Days_operation_treatment[match(sample_table$sampleName, meta_data$Sample)]
sample_table$days_arrive_treatment <- meta_data$Days_arrive_treatment[match(sample_table$sampleName, meta_data$Sample)]
sample_table$infusion_startime <- meta_data$Infusion_starttime[match(sample_table$sampleName, meta_data$Sample)]
sample_table$infusion_stoptime <- meta_data$Infustion_stoptime[match(sample_table$sampleName, meta_data$Sample)]
sample_table$infusion_timelength <- meta_data$Infusion_timelength[match(sample_table$sampleName, meta_data$Sample)]
sample_table$infusion_timelength_mins <- meta_data$Infusion_timelength_minutes[match(sample_table$sampleName, meta_data$Sample)]
sample_table$group <- meta_data$Group[match(sample_table$sampleName, meta_data$Sample)]
sample_table$original_condition <- meta_data$condition[match(sample_table$sampleName, meta_data$Sample)]
sample_table$arrival_date <- meta_data$Arrival_date[match(sample_table$sampleName, meta_data$Sample)]
sample_table$notes <- meta_data$Notes[match(sample_table$sampleName, meta_data$Sample)]






