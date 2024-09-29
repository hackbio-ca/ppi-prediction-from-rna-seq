suppressPackageStartupMessages({
  library(dplyr)
})

# define input/output paths and read
johnson_path_hek <- "ppi-prediction-from-rna-seq/output/johnson2021.preprocessed_PPIs.HEK293T.tsv"
johnson_path_jurkat <- "ppi-prediction-from-rna-seq/output/johnson2021.preprocessed_PPIs.jurkat.tsv"
johnson_path_huvec <- "ppi-prediction-from-rna-seq/output/johnson2021.preprocessed_PPIs.HUVEC.tsv"
khoroshkin_path <- "ppi-prediction-from-rna-seq/output/khoroshkin2024.preprocessed_PPIs.K562.tsv"
huttlin_path <- "ppi-prediction-from-rna-seq/output/huttlin2021.preprocessed_PPIs.HEK293T.tsv"
goos_path <- "ppi-prediction-from-rna-seq/output/goos2022.preprocessed_PPIs.HEK293.tsv"

johnson_hek <- read.csv(johnson_path_hek, sep = "\t")
johnson_jurkat <- read.csv(johnson_path_jurkat, sep = "\t")
johnson_huvec <- read.csv(johnson_path_huvec, sep = "\t")
khoroshkin <- read.csv(khoroshkin_path, sep = "\t")
huttlin <- read.csv(huttlin_path, sep = "\t")
goos <- read.csv(goos_path, sep = "\t")

output_johnson_hek <- "ppi-prediction-from-rna-seq/output/johnson2021.preprocessed_PPIs.HEK293T.standardized.tsv"
output_johnson_jurkat <- "ppi-prediction-from-rna-seq/output/johnson2021.preprocessed_PPIs.jurkat.standardized.tsv"
output_johnson_huvec <- "ppi-prediction-from-rna-seq/output/johnson2021.preprocessed_PPIs.HUVEC.standardized.tsv"
output_khoroshkin <- "ppi-prediction-from-rna-seq/output/khoroshkin2024.preprocessed_PPIs.K562.standardized.tsv"
output_huttlin <- "ppi-prediction-from-rna-seq/output/huttlin2021.preprocessed_PPIs.HEK293T.standardized.tsv"
output_goos <- "ppi-prediction-from-rna-seq/output/goos2022.preprocessed_PPIs.HEK293.standardized.tsv"

# transform p values to -logp
johnson_df_hek <- johnson_hek %>%
  mutate(pval_fixed = ifelse(pvalue == 0, .Machine$double.xmin, pvalue)) %>%
  mutate(logp = -log10(pval_fixed))
johnson_df_jurkat <- johnson_jurkat %>%
  mutate(pval_fixed = ifelse(pvalue == 0, .Machine$double.xmin, pvalue)) %>%
  mutate(logp = -log10(pval_fixed))
johnson_df_huvec <- johnson_huvec %>%
  mutate(pval_fixed = ifelse(pvalue == 0, .Machine$double.xmin, pvalue)) %>%
  mutate(logp = -log10(pval_fixed))

# Step 1: Rank the values (ignoring NA values, and using the average for ties)
ranked_johnson_temp_hek <- rank(johnson_df_hek$logp, na.last = "keep", ties.method = "average")
ranked_johnson_temp_jurkat <- rank(johnson_df_jurkat$logp, na.last = "keep", ties.method = "average")
ranked_johnson_temp_huvec <- rank(johnson_df_huvec$logp, na.last = "keep", ties.method = "average")
ranked_khoroshkin_temp <- rank(khoroshkin$conf_score, na.last = "keep", ties.method = "average")
ranked_huttlin_temp <- rank(huttlin$pscore, na.last = "keep", ties.method = "average")
ranked_goos_temp <- rank(goos$conf_score, na.last = "keep", ties.method = "average")

# Step 2: Standardize the ranked values by dividing by the maximum rank (excluding NA)
standardized_ranked_johnson_hek <- ranked_johnson_temp_hek / max(ranked_johnson_temp_hek, na.rm=T)
standardized_ranked_johnson_jurkat <- ranked_johnson_temp_jurkat / max(ranked_johnson_temp_jurkat, na.rm=T)
standardized_ranked_johnson_huvec <- ranked_johnson_temp_huvec / max(ranked_johnson_temp_huvec, na.rm=T)
standardized_ranked_khoroshkin <- ranked_khoroshkin_temp / max(ranked_khoroshkin_temp, na.rm=T)
standardized_ranked_huttlin <- ranked_huttlin_temp / max(ranked_huttlin_temp, na.rm=T)
standardized_ranked_goos <- ranked_goos_temp / max(ranked_goos_temp, na.rm=T)

# Step 3: cbind babey
johnson_ranks_hek <- cbind(johnson_df_hek, standardized_ranked_johnson_hek) %>% rename(standardized_rank = "standardized_ranked_johnson_hek")
johnson_ranks_jurkat <- cbind(johnson_df_jurkat, standardized_ranked_johnson_jurkat) %>% rename(standardized_rank = "standardized_ranked_johnson_jurkat")
johnson_ranks_huvec <- cbind(johnson_df_huvec, ranked_johnson_temp_huvec) %>% rename(standardized_rank = "ranked_johnson_temp_huvec")
khoroshkin_ranks <- cbind(khoroshkin, standardized_ranked_khoroshkin) %>% rename(standardized_rank = "standardized_ranked_khoroshkin")
huttlin_ranks <- cbind(huttlin, standardized_ranked_huttlin) %>% rename(standardized_rank = "standardized_ranked_huttlin")
goos_ranks <- cbind(goos, standardized_ranked_goos) %>% rename(standardized_rank = "standardized_ranked_goos")

# write to file
write.table(johnson_ranks_hek, output_johnson_hek, quote = F, sep = "\t", row.names = F)
write.table(johnson_ranks_huvec, output_johnson_jurkat, quote = F, sep = "\t", row.names = F)
write.table(johnson_ranks_huvec, output_johnson_huvec, quote = F, sep = "\t", row.names = F)
write.table(khoroshkin_ranks, output_khoroshkin, quote = F, sep = "\t", row.names = F)
write.table(huttlin_ranks, output_huttlin, quote = F, sep = "\t", row.names = F)
write.table(goos_ranks, output_goos, quote = F, sep = "\t", row.names = F)

