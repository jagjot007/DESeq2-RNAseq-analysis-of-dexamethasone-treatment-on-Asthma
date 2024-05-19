#working directory swt
getwd()
setwd("~/Desktop/data science /DEseq2 project")

#installing and loading the packages
install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

if (!"tidyverse" %in% rownames(installed.packages())) {
  message("Package 'tidyverse' is not installed.")
} else {
  message("Package 'tidyverse' is already installed.")
}
install.packages("tidyverse")
library(tidyverse)
library(readr)

BiocManager::install("airway")
library(airway)

#looking at the data
data(airway)
str(airway)
View(airway)
airway 

#preparing the sample data 
sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]    #we only need column 2 and 3

sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
sample_info$dexamethasone <- as.factor(sample_info$dexamethasone)  #make them categoriacal so that helps later in DESeq

#wrangling the assay data (raw gene expression couunts)
countsData <- assay(airway)
countsdata <- as.data.frame(countsData)    #easier to work with data frames
head(countsdata)

#cleaning the environment
rm(airway)
rm(countsData)

#Checking if the row names in ColumnData match the column names in counts_data
if (all(colnames(countsdata) %in% rownames(sample_info))) {
  print("yes") 
} else { print("No") }

all(colnames(countsdata) == rownames(sample_info))

#constructing a DESeq dataset object
dds <- DESeqDataSetFromMatrix(countData = countsdata,
                       colData = sample_info,
                       design = ~ dexamethasone)
dds

#pre-processing step :- removing low read count data genes (keeping read count min limit as 10)
keep <- rowSums(counts(dds)) >= 10       #counts is a function here
keep
dds <- dds[keep,]
dds        #dimensions decreased

#setting up the factor/reference level (telling DESeq that untreated in design is our reference level)
#otherwise it chooses ref level alphabetically by default
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")
dds$dexamethasone

#NOTE: if there are technical replicates then we need to collapse them before performing DESeq (never collapse biological replicates)
#running DESeq
dds <- DESeq(dds)

#exploring results
res <- results(dds)
res
summary(res)    #it is using p value of 0.1 which is high so we reduce it
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

#contrasts
resultsNames(dds)

#for multiple design factors we can compare each category individually
#results(dds, contrast = c("dexamethasone","treatment1","untreated"))
#results(dds, contrast = c("dexamethasone","treatment2","untreated"))
#results(dds, contrast = c("dexamethasone","treatment2","treatment1"))
#we already defined a reference in our case therefore this is not needed.

# MA plot
plotMA(res0.01)



