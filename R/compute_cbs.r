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
compute_cbs <- function(coexpr_matrix, bioid_matrix, top_n, cell_line) {

  coexpr_matrix <- read_csv(coexpr_matrix) %>%
    as.matrix()
  rownames(coexpr_matrix) <- colnames(coexpr_matrix)

  bioid_matrix <- read_csv(bioid_matrix) %>%
    column_to_rownames("gene2") %>%
    as.matrix()

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
    top_n_coexpr_genes <- names(coexpr_ranks[coexpr_ranks <= top_n])

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
  write(predictive_scores, paste0("../output/", cell_line, "_cbs.csv"))
  return(predictive_scores)
}

cell_lines <- c("hek293", "hek293t", "k562", "jurkat", "huvec")
## note: we have hek293t in both Johnson and Huttlin datasets

lapply(cell_lines, function(cell_line) {
  compute_cbs(paste0("RNA/all_", cell_line, "_subset.csv"),
              paste0("PPI/", cell_line, "_ppi_matrix.csv"),
              top_n = 100,
              cell_line = cell_line)
})
