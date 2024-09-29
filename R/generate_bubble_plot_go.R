library(ggplot2)
library(tidyverse)
library(gprofiler2)

#general pipeline for inputting cba analysis to do low hanging go enrichment analysis

cba_plot <- function(input_file, cell_line){
df <- read.table(input_file, sep = ",", header = T)
df_cutoff <- filter(df ,auroc_score > .80)
df_cutoff_genes <- df_cutoff$bait_gene
gostres_cell_line <- gost(query = df_cutoff_genes,
                organism = "hsapiens")

gostres_cell_line_res <- gostres_cell_line$result
gostres_cell_line_plot <- dplyr::select(gostres_cell_line_res, p_value, intersection_size, term_name)
gostres_cell_line_plot$cell_line <- c(cell_line)

return(gostres_cell_line_plot)}

#append all the outputs, colored by cell 
cell_line <- c("k62")
k562_for_graph <- cba_plot("k562_cbs.csv", "k562")

all_cell_lines <- cbind()

p <- ggplot(all_cell_lines, aes(x=term_name , y=-log10(p_value), size = intersection_size, color = cell_line)) +
  geom_point(alpha=0.5) +
  scale_size(range = c(1, 10), name="interaction size") +
  theme_bw()
