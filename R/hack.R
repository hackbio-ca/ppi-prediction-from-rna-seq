library(tidyverse)
library(biomaRt)
library(MetaMarkers)
#####huvec
hek293t <- read.table("GSE233957_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t", header = T)
df <- hek293t

cl_control <- as.matrix(df[,c(1, 2:4, 23:25)])

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
entrez_id_cl_full <-  entrez_id_cl[!(is.na(entrez_id_cl$hgnc_symbol) | entrez_id_cl$hgnc_symbol==""), ]


entrez_id_cl_full = entrez_id_cl_full[!duplicated(entrez_id_cl_full$hgnc_symbol),]

rownames(entrez_id_cl_full) <- entrez_id_cl_full[,9]
entrez_id_cl_full <- entrez_id_cl_full[,2:7]
cl_cpm_fil_transposed <- t(entrez_id_cl_full)
cl_cpm_fil_t_cor <- cor(cl_cpm_fil_transposed)

write.table(cl_cpm_fil_t_cor, "hek293t_corr_mat.csv", sep = ",", quote = F)
cl_cpm_fil_t_cor_1 <- as.data.frame(cl_cpm_fil_t_cor[1:10000,])

write.table(cl_cpm_fil_t_cor_1, file = "hek293t_corr_mat.csv", sep = ",", col.names = F)


cl_cpm_fil_t_cor_small <- cl_cpm_fil_t_cor[1:1000, 1:1000]
write.table(cl_cpm_fil_t_cor_small, "small_hek_corr_mat.csv", sep = ",", quote = F)
genes_1k <- colnames(cl_cpm_fil_t_cor_small)
write.table(genes_1k, "genes_1k_jurkat.tsv", sep = "\t", quote = F)


genes_23k <- colnames(jurkat_cpm_fil_t_cor)
write.table(genes_23k, "genes_23k.tsv", sep = "\t", quote = F)

#match to protein protein
read.table("GSE233957_raw_counts_GRCh38.p13_NCBI.tsv", sep = "\t", header = T)