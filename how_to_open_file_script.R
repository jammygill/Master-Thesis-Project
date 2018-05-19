getwd()
setwd("H:/Master_Thesis/Datasets/study1_naive_diff")
library(readxl)

# file <- read.csv(file = "metadata.csv", header = TRUE, sep = ",")
# file
# file2 <- read_excel(path = "metadata_naiverats_envigo_without_semi.xlsx", sheet = 1, col_names = TRUE)

filelist = list.files(pattern = ".*.txt")
datalist = lapply(filelist, function(x)read.table(x, header=F))
datafr = do.call("cbind", datalist)

dim.data.frame(datafr)

View(datafr)

class(datafr)
