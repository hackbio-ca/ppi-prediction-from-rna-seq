suppressPackageStartupMessages({
  library(dplyr)
})

input_path <- "ppi-prediction-from-rna-seq/data/ppi/PROPER_v1.csv"
output_path <- "ppi-prediction-from-rna-seq/output/johnson2021.preprocessed_PPIs.tsv"

proper <- read.csv(input_path)

# Remove interactions that are not significant in BH corrected p values
proper_sig <- proper %>%
    filter(BH.corrected.p.value < 0.05) %>%
    rename(gene1 = "Gene1", gene2 = "Gene2", pvalue = "BH.corrected.p.value") %>%
    mutate(source = "Johnson2021")

proper_sig_hek <- proper_sig %>%
    filter(Cell.line.specificity == "HEK293" | Cell.line.specificity == "shared") %>%
    mutate(cell_line = "HEK293T") %>%
    select(gene1, gene2, pvalue, cell_line, source)

proper_sig_jurkat <- proper_sig %>%
    filter(Cell.line.specificity == "Jurkat" | Cell.line.specificity == "shared") %>%
    mutate(cell_line = "Jurkat") %>%
    select(gene1, gene2, pvalue, cell_line, source)

proper_sig_huvec <- proper_sig %>%
    filter(Cell.line.specificity == "HUVEC" | Cell.line.specificity == "shared") %>%
    mutate(cell_line = "HUVEC") %>%
    select(gene1, gene2, pvalue, cell_line, source)

# Decided to merge all df together to store in a single file
proper_sig_all <- rbind(proper_sig_hek, proper_sig_jurkat, proper_sig_huvec)

# save to file
write.table(proper_sig_all, output_path, row.names = FALSE, quote = F, sep = "\t", row.names = F)
