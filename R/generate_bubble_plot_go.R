library(ggplot2)
library(tidyverse)
library(gprofiler2)

#general pipeline for inputting cba analysis to do low hanging go enrichment analysis

cba_plot <- function(input_file, cell_line){
df <- read.table(input_file, sep = ",", header = T)
df_cutoff <- filter(df ,auroc_score > .7)
df_cutoff_genes <- df_cutoff$bait_gene
gostres_cell_line <- gost(query = df_cutoff_genes,
                organism = "hsapiens")

gostres_cell_line_res <- gostres_cell_line$result
gostres_cell_line_plot <- dplyr::select(gostres_cell_line_res, p_value, intersection_size, term_name)
gostres_cell_line_plot$cell_line <- c(cell_line)

return(gostres_cell_line_plot)}

#append all the outputs, colored by cell 
k562_for_graph <- cba_plot("final_cbs/k562_cbs.csv", "k562")
hek293_for_graph <- cba_plot("final_cbs/hek293_cbs.csv", "hek293")
hek293t_for_graph <- cba_plot("final_cbs/hek293T_j_cbs.csv", "hek293t")
jurkat_for_graph <- cba_plot("final_cbs/jurkat_cbs.csv", "jurkat")
huvec_for_graph <- cba_plot("final_cbs/huvec_cbs.csv", "huvec")

all_cell_lines <- rbind(k562_for_graph, hek293_for_graph, hek293t_for_graph, jurkat_for_graph, huvec_for_graph)

p <- ggplot(all_cell_lines, aes(x=term_name , y=-log10(p_value), size = intersection_size, color = cell_line)) +
  geom_point(alpha=0.2) +
  scale_size(range = c(1, 10), name="interaction size") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p

ggsave("../output/corr.png", width = 50, height = 50, units = "cm", limitsize = F)

k562 <- read.table("../data/final_cbs/k562_cbs.csv", sep = ",", header = T)
hek293 <- read.table("../data/final_cbs/hek293_cbs.csv", sep = ",", , header = T)
hek293t<- read.table("../data/final_cbs/hek293T_j_cbs.csv", sep = ",", , header = T)
jurkat <- read.table("../data/final_cbs/jurkat_cbs.csv", sep = ",", , header = T)
huvec<- read.table("../data/final_cbs/huvec_cbs.csv", sep = ",", , header = T)

all_cell_lines_hist <- rbind(k562, hek293, hek293t, jurkat, huvec)


ggplot(all_cell_lines_hist, aes(x = auroc_score)) +    
  geom_histogram(alpha = 0.5, fill = "#69b3a2") + theme_bw()

ggsave("../output/hist.png", width = 50, height = 50, units = "cm", limitsize = F)

