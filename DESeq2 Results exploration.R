getwd()
setwd("~/Desktop/data science /DEseq2 project")

#After DESeq2 RNAseq analysis

#finding out the top 10 most differentially expresssed genes on treatment
#sorting the results data
sorted_res0.01 <- res0.01[order(res0.01$padj), ]
#top 10 genes
top10_genes <- head(sorted_res0.01, 10)
#extracting gene names
top10_gene_names <- rownames(top10_genes)
top10_gene_names

#we got the Ensemble IDs of the genes but we want the gene names
BiocManager::install("biomaRt")
library(biomaRt)
#connecting to Ensembl database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#getting gene names
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = top10_gene_names,
                   mart = mart)
print(gene_info)

#there are some genes of interest :- DPP10, PHF11, ADAM33, NPSR1, HLA-G
#these genes are commonly associated with Asthma
#we want to check if there is any differential expression in our genes of interest
#making list of genes of interest
genes_of_interest <- c("DPP10", "PHF11", "ADAM33", "NPSR1", "HLA-G")
#converting gene names to Ensembl IDs
gene_ids <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                  filters = "external_gene_name",
                  values = genes_of_interest,
                  mart = mart)
print(gene_ids)
ensembl_ids_of_interest <- gene_ids$ensembl_gene_id
#checking results of our Ensembl IDs
res_of_interest <- res0.01[rownames(res0.01) %in% ensembl_ids_of_interest, ]
print(res_of_interest)

#3 genes were found in the data :- PHF11, ADAM33 and HLA-G
#none of them were significant


