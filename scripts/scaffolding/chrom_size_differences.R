library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

# Read in genome index
CGC2 <- readr::read_tsv("../../processed_data/genomes/c_briggsae.CGC2.hifi.ONT.HiC.Feb2026.genome.fa.fai", col_names = c("chrom","CGC2","x1","x2","x3")) %>%
  dplyr::select(chrom,CGC2) %>%
  dplyr::filter(chrom != "AC186293.1")

# Read in genome index
AF16 <- readr::read_tsv("../../processed_data/genomes/c_briggsae.AF16.PRJNA10731.WS276.genome.fa.fai", col_names = c("chrom","AF16","x1","x2","x3")) %>%
  dplyr::select(chrom,AF16) %>%
  dplyr::filter(!grepl("cb",chrom)) 

both_chrom_sizes <- AF16 %>%
  dplyr::left_join(CGC2, by = 'chrom') %>%
  pivot_longer(
    cols = c(AF16, CGC2),
    names_to = "identity",
    values_to = "length")

# Calculating the difference in chromosome size of CGC2 compared to AF16
size_diff <- AF16 %>%
  dplyr::left_join(CGC2, by = 'chrom') %>%
  dplyr::mutate(diff_size_kb = round((CGC2 - AF16) / 1000, 2))

# Plotting 
chrom_size_differences <- ggplot(both_chrom_sizes, aes(x = chrom, y = length / 1e6, fill = identity)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(data = size_diff, aes(x = chrom, y = pmax(AF16, CGC2) / 1e6 + 0.3, label = paste0(round((CGC2 - AF16) / 1000, 1), " kb"), fill = NULL), size = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.025))) +
  scale_fill_manual(values = c("AF16" = "#E6B3B3", "CGC2" = "#7570B3")) +
  theme(
    axis.title = element_text(size = 12, color = 'black'),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, color = 'black'),
    panel.background = element_blank(),
    panel.grid.major= element_line(color = 'gray80'),
    legend.text = element_text(size = 11, color = 'black'),
    legend.title = element_text(size = 11, color = 'black'),
    legend.box.background = element_rect(color = "black", fill = NA),
    legend.position = "inside",
    legend.position.inside = c(0.15,0.9),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 10, color = 'black'),
    panel.border = element_rect(fill = NA)) +
  labs(y = "Chromosome size (Mb)", fill = "Genome assembly")
chrom_size_differences


# Save plot
ggsave("../../figures/supplementary/chrom_size_differences.png", chrom_size_differences, width = 7, height = 7, dpi = 600)
