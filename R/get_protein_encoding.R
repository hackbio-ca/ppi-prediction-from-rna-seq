library(tidyverse)
library(biomaRt)
###get protein encoding only genes from expression matrices

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")                                                                         
filters = listFilters(ensembl)                                                                                                        
myResult <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "hgnc_symbol"), mart=ensembl)
protein_coding <- filter(myResult, gene_biotype == "protein_coding")
protein_coding_list <- protein_coding$hgnc_symbol
#where cl_cpm_fil_t_cor is a matrix with col and row names of the gene symbols
all_hek293t_subset <- cl_cpm_fil_t_cor[rownames(cl_cpm_fil_t_cor)%in%protein_coding_list,colnames(cl_cpm_fil_t_cor)%in%protein_coding_list]

cl_cpm_fil_t_cor_small_hek293t <- all_hek293t_subset[1:100, 1:100]
write.table(all_hek293t_subset, "../output/all_hek293t_subset.csv", sep = ",", quote = F)

#for input k562
k562 <- read.table("../data/perturb_rbp_coexp.csv/perturb_rbp_coexp.csv", sep = ",", header = T)
k562_genes <- colnames(k562)
k562$genes <- k562_genes
rownames(k562) <- k562[,17893]

all_k562_subset <- k562[rownames(k562)%in%protein_coding_list,colnames(k562)%in%protein_coding_list]

all_k562_subset_small <- all_k562_subset[1:10, 1:10]
write.table(all_k562_subset, "../output/all_k562_subset.csv", sep = ",", quote = F)