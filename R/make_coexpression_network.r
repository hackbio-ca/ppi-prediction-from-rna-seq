# Packages
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(SingleCellExperiment)
  library(MetaMarkers)
})

# Load the data
sample_to_rbp_mapping <- read_tsv("perturb_rbp/GSM7056650_bc_to_sgrna_mapping_3_06.csv")
## perturb_rbp is a directory containing the data downloaded from GSE225807
perturb_rbps <- Read10X(data.dir = "perturb_rbp") %>%
  CreateSeuratObject()

colnames(perturb_rbps) <- str_remove(colnames(perturb_rbps), "-1")
perturb_rbps_sub <- subset(perturb_rbps, cells = sample_to_rbp_mapping$barcode)

sample_to_rbp_mapping_sub <- sample_to_rbp_mapping %>%
  filter(barcode %in% colnames(perturb_rbps_sub)) %>%
  .[match(colnames(perturb_rbps_sub), .$barcode), ]

perturb_rbps_sub <- AddMetaData(perturb_rbps_sub, sample_to_rbp_mapping_sub)

# Get pseudo-bulk
get_cpm_pseudo_bulk <- function(scdata, new_labels, nbds) {
  sc_sub2 <- matrix(NA, nrow = dim(scdata)[1], ncol = length(nbds))

  for (i in seq_along(nbds)){
    cols <- which(new_labels == nbds[i])

    if(length(cols) > 1) {
      sc_sub2[, i] <- rowSums(cpm(scdata)[, cols]) / length(cols)
    } else if (length(cols) == 1) {
      sc_sub2[, i] <- cpm(scdata)[, cols]
    } else {
      sc_sub2[, i] <- NA
    }
  }
  rownames(sc_sub2) <- rownames(scdata)
  colnames(sc_sub2) <- nbds
  return(sc_sub2)
}

perturb_rbps_sce <- as.SingleCellExperiment(perturb_rbps_sub)

assay(perturb_rbps_sce, "cpm") <- MetaMarkers::convert_to_cpm(assay(perturb_rbps_sce))
perturb_rbps_bulk <- get_cpm_pseudo_bulk(perturb_rbps_sce, perturb_rbps_sce$RBP, 
                                         unique(perturb_rbps_sce$RBP))

# Remove genes with no counts
perturb_rbps_bulk_filtered <- perturb_rbps_bulk[rowSums(perturb_rbps_bulk) > 0, ]

# Coexpression network
perturb_rbps_coexp <- cor(t(perturb_rbps_bulk_filtered))

write_csv(as.data.frame(perturb_rbps_coexp), "perturb_rbp_coexp.csv")