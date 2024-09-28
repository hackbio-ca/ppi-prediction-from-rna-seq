library(tidyverse)
library(biomaRt)
#####huvec
huvec <- read.table("GSE221514_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t", header = T)
df <- huvec
  
cl_control <- as.matrix(df[,1:7])

# if (!requireNamespace("devtools", quietly = TRUE))
#   install.packages("devtools")
# devtools::install_github("gillislab/MetaMarkers")
# library(MetaMarkers)

cl_cpm <- convert_to_cpm(cl_control[,2:7])
cl_cpm <- cbind(cl_control[,1], cl_cpm)
cl_cpm_fil <- cl_cpm[ rowSums(cl_cpm[,-1]) > 0, ]

entrez_id_cl <- as.data.frame(cl_cpm_fil[,1])
colnames(entrez_id_cl) <- c("V1")


ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")                                                                         
filters = listFilters(ensembl)                                                                                                        
entrezgene = (entrez_id_cl$V1)             
genes <- getBM(filters="entrezgene_id", attributes=c("ensembl_gene_id","entrezgene_id", "hgnc_symbol"), values=entrezgene, mart=ensembl)                                                                                                                 


entrez_id_cl <- left_join(as.data.frame(cl_cpm_fil), genes, by = c("V1" = "entrezgene_id"))
entrez_id_cl_full <- entrez_id_cl[complete.cases(entrez_id_cl), ]

entrez_id_cl_full = entrez_id_cl_full[!duplicated(entrez_id_cl_full$ensembl_gene_id),]

rownames(entrez_id_cl_full) <- entrez_id_cl_full[,8]
entrez_id_cl_full <- entrez_id_cl_full[,2:7]
cl_cpm_fil_transposed <- t(entrez_id_cl_full)
cl_cpm_fil_t_cor <- cor(cl_cpm_fil_transposed)

write.table(cl_cpm_fil_t_cor, "huvec_corr_mat.csv", sep = ",", quote = F)
