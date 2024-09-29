library(dplyr)

johnson_path <- "ppi-prediction-from-rna-seq/output/johnson2021.preprocessed_PPIs.tsv"
khoroshkin_path <- "ppi-prediction-from-rna-seq/output/khoroshkin2024.preprocessed_PPIs.tsv"
huttlin_path <- "ppi-prediction-from-rna-seq/output/huttlin2021.preprocessed_PPIs.tsv"
goos_path <- "ppi-prediction-from-rna-seq/output/goos2022.preprocessed_PPIs.tsv"

johnson <- read.csv(johnson_path, sep = "\t")
khoroshkin <- read.csv(khoroshkin_path, sep = "\t")
huttlin <- read.csv(huttlin_path, sep = "\t")
goos <- read.csv(goos_path, sep = "\t")

output_johnson <- "ppi-prediction-from-rna-seq/output/johnson2021.preprocessed_PPIs.standardized.tsv"
output_khoroshkin <- "ppi-prediction-from-rna-seq/output/khoroshkin2024.preprocessed_PPIs.standardized.tsv"
output_huttlin <- "ppi-prediction-from-rna-seq/output/huttlin2021.preprocessed_PPIs.standardized.tsv"
output_goos <- "ppi-prediction-from-rna-seq/output/goos2022.preprocessed_PPIs.standardized.tsv"

# transform p values to -logp
johnson_df <- johnson %>%
  mutate(pval_fixed = ifelse(pvalue == 0, 0.0000001, pvalue)) %>%
  mutate(logp = -log10(pval_fixed))

# Step 1: Rank the values (ignoring NA values, and using the average for ties)
ranked_johnson_temp <- rank(johnson_df$logp, na.last = "keep", ties.method = "average")
ranked_khoroshkin_temp <- rank(khoroshkin$conf_score, na.last = "keep", ties.method = "average")
ranked_huttlin_temp <- rank(huttlin$pscore, na.last = "keep", ties.method = "average")
ranked_goos_temp <- rank(goos$conf_score, na.last = "keep", ties.method = "average")


# Step 2: Standardize the ranked values by dividing by the maximum rank (excluding NA)
standardized_ranked_johnson <- ranked_johnson_temp / max(ranked_johnson_temp, na.rm=T)
standardized_ranked_khoroshkin <- ranked_khoroshkin_temp / max(ranked_khoroshkin_temp, na.rm=T)
standardized_ranked_huttlin <- ranked_huttlin_temp / max(ranked_huttlin_temp, na.rm=T)
standardized_ranked_goos <- ranked_goos_temp / max(ranked_goos_temp, na.rm=T)

# Step 3: cbind babey
johnson_ranks <- cbind(johnson_df, standardized_ranked_johnson) %>% rename(standardized_rank = "standardized_ranked_johnson")
khoroshkin_ranks <- cbind(khoroshkin, standardized_ranked_khoroshkin) %>% rename(standardized_rank = "standardized_ranked_khoroshkin")
huttlin_ranks <- cbind(huttlin, standardized_ranked_huttlin) %>% rename(standardized_rank = "standardized_ranked_huttlin")
goos_ranks <- cbind(goos, standardized_ranked_goos) %>% rename(standardized_rank = "standardized_ranked_goos")

# write to file
write.table(johnson_ranks, output_johnson, quote = F, sep = "\t", row.names = F)
write.table(khoroshkin_ranks, output_khoroshkin, quote = F, sep = "\t", row.names = F)
write.table(huttlin_ranks, output_huttlin, quote = F, sep = "\t", row.names = F)
write.table(goos_ranks, output_goos, quote = F, sep = "\t", row.names = F)

