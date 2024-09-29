args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
})

input_path <- args[1]
# extract author and year name from the input and use for output naming
author_year <- strsplit(basename(input_path), "_")[[1]][1]
output_path <- paste0("ppi-prediction-from-rna-seq/output/", author_year, ".preprocessed_PPIs.csv")

# Lazy case handling of getting the cell lines... because we know what datasets 
# were dealing with
if (author_year == "khoroshkin2024") {
    current_cell_line <- "K562"
    } else if (author_year == "goos2022") {
    current_cell_line <- "HEK293"
    } else {
    print("You messed this up somehow and aren't using the right input files.")
    stop()
    }

mitab_raw <- read.csv(input_path, sep = "\t")

# Custom fx to get the gene name from the Alias.es..interactor.A and B columns
get_geneSymbol <- function(alias) {
    if(is.na(alias)) {
        return(NA)
    }
    divided_aliases <- strsplit(alias, "\\|")[[1]]

    uniprot_gene_name <- divided_aliases[grep("(gene name)", divided_aliases, fixed = T)]
    if(length(uniprot_gene_name) == 0) {
        return(NA)
    }
    # grab the gene name that is between ":" and "("
    gene_name_noUniprot <- gsub(".*:", "", uniprot_gene_name)
    gene_name <- gsub("\\(.*", "", gene_name_noUniprot)

    return(gene_name)
}

# Extract gene name and organism from columns 
# gene name will be found in Alias.es..interactor.A and B
mitab_formatted <- mitab_raw %>%
    mutate(gene1 = lapply(Alias.es..interactor.A, FUN=get_geneSymbol),
           gene2 = lapply(Alias.es..interactor.B, FUN=get_geneSymbol)) %>%
    rename(conf_score = "Confidence.value.s.") %>%
    mutate(conf_score = as.numeric(gsub("[^0-9.]", "", conf_score))) %>%
    select(gene1, gene2, conf_score) %>%
    mutate(cell_line = current_cell_line,
           source = author_year)

mitab_final <- data.frame(lapply(mitab_formatted, as.character), stringsAsFactors=FALSE)

write.csv(mitab_final, output_path, row.names = F)
