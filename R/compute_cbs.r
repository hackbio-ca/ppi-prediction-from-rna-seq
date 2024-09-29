# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(Matrix)
  library(pROC)
})

# Compute CBS from gene-gene and protein-protein matrices
compute_cbs <- function(coexpr_matrix, bioid_matrix, top_n, line,
                        fix_format = TRUE, ppi_stat = "score") {

  if (fix_format) {
    coexpr_matrix <- fix_data(read_csv(coexpr_matrix))
  } else {
    coexpr_matrix <- read_csv(coexpr_matrix) %>%
      as.matrix()
    rownames(coexpr_matrix) <- colnames(coexpr_matrix)
  }

  if (ppi_stat == "score") {
    bioid_matrix <- read_csv(bioid_matrix) %>%
      column_to_rownames("gene2") %>%
      as.matrix()
  } else {
    bioid_matrix <- read_csv(bioid_matrix) %>%
      column_to_rownames("gene2") %>%
      mutate(across(everything(), ~ -log10(. + .Machine$double.xmin))) %>% 
      as.matrix()
  }

  # Initialize result storage
  predictive_scores <- data.frame(bait_gene = character(),
                                  auroc_score = numeric(),
                                  stringsAsFactors = FALSE)

  # Loop through each bait gene (which is also the bait protein)
  for (bait_gene in rownames(coexpr_matrix)) {

    # Check if bait gene exists in the BioID matrix as a column (bait)
    if (!bait_gene %in% colnames(bioid_matrix)) {
      next  # Skip if bait gene is not in BioID matrix
    }

    # Step 1: Rank coexpression values for the bait gene
    coexpr_values <- coexpr_matrix[bait_gene, ]
    coexpr_ranks <- rank(-coexpr_values, ties.method = "average")

    # Step 2: Rank prey proteins by confidence scores
    prey_scores <- bioid_matrix[, bait_gene]
    prey_ranks <- rank(-prey_scores, ties.method = "average")
    top_n_prey_genes <- names(prey_ranks[prey_ranks <= top_n])

    # Step 3: Create a binary label vector: 1 if gene is in the top N preys,
    # 0 otherwise
    labels <- ifelse(rownames(coexpr_matrix) %in% top_n_prey_genes, 1, 0)

    # Ensure that labels and coexpression ranks are aligned
    if (!all(names(coexpr_ranks) == rownames(coexpr_matrix))) {
      stop("Mismatch between coexpression ranks and gene names")
    }

    # Step 4: Compute AUROC-like score based on how well top N
    # coexpression partners predict top N interaction partners
    roc_result <- roc(labels, coexpr_ranks, quiet = TRUE)
    auroc_score <- as.numeric(auc(roc_result))

    # Store the result
    predictive_scores <- rbind(predictive_scores,
                               data.frame(bait_gene = bait_gene,
                                          auroc_score = auroc_score,
                                          stringsAsFactors = FALSE))
  }
  #write_csv(predictive_scores, paste0("../output/", cell_line, "_cbs.csv"))
  write_csv(predictive_scores, paste0("final_cbs/", line, "_cbs.csv"))
  return(predictive_scores)
}

# Resolve coexpression matrix format issue
fix_data <- function(data) {
  num <- length(colnames(data)) - 1
  keep <- colnames(data)[1:num - 1]
  data <- data[1:(num - 1), 2:num]
  data <- as.matrix(data)
  colnames(data) <- keep
  rownames(data) <- keep
  return(data)
}

## note: we have hek293t in both Johnson and Huttlin datasets

# HEK293T (Huttlin)
hek293T_h_cbs <- compute_cbs(paste0("RNA/all_hek293t_subset.csv"),
                             paste0("PPI/hek293T_h_ppi_matrix.csv"),
                             top_n = 100,
                             line = "hek293T_h")

# HEK293
hek293_cbs <- compute_cbs(paste0("RNA/all_hek293_subset.csv"),
                          paste0("PPI/hek293_ppi_matrix.csv"),
                          top_n = 100,
                          line = "hek293")

# K562
k562_cbs <- compute_cbs(paste0(paste0("RNA/all_k562_subset.csv")),
                        paste0("PPI/k562_ppi_matrix.csv"),
                        top_n = 10,
                        line = "k562")

# johnson datasets
## HEK293T
compute_cbs(paste0("RNA/all_hek293t_subset.csv"),
            paste0("PPI/hek293T_j_ppi_matrix.csv"),
            top_n = 10,
            line = "hek293T_j",
            ppi_stat = "pval")

## HUVEC, JURKAT
lapply(c("huvec", "jurkat"), function(line) {
  compute_cbs(paste0("RNA/all_", line, "_subset.csv"),
              paste0("PPI/", line, "_ppi_matrix.csv"),
              top_n = 10,
              line = line,
              ppi_stat = "pval")
})