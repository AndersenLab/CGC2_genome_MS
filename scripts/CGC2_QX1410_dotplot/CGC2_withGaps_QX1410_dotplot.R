library(ggplot2)
library(readr)
library(dplyr)
library(cowplot)

# Load in BED file containing coordinates of assembly gaps
gaps <- readr::read_tsv("../../processed_data/scaffolding/CGC2_gaps.bed", col_names = c("CGC2_chrom", "start", "end"))

# Load in genome-genome alignments of CGC2 with gaps present against QX1410
alignment_withGaps <- readr::read_tsv("../../processed_data/genome_genome_alignments/CGC2_withGaps_QX1410.transformed.tsv", col_names = c("QXS","QXE","CGS","CGE","L1","L2","IDY","LENR","LENQ","QX_chrom","CGC2_chrom"))

# Gaps on V
CG_telomeres <- readr::read_tsv("../../processed_data/genomes/CGC2_withGaps_telomeres_binned_1kb.bed", col_names = c("chrom","start","end","count")) %>%
  dplyr::filter(chrom == "V" & count > 150) %>% # Represents 1 kb segments that are entirely telomeric repeats
  dplyr::filter(start < 150000) %>%
  dplyr::mutate(start = min(start), end = max(end)) %>%
  dplyr::distinct(chrom, start, end)

chromV_gaps <- ggplot(alignment_withGaps %>% dplyr::filter(CGC2_chrom == "V" & QX_chrom == "V")) +
  geom_rect(data = gaps %>% dplyr::filter(CGC2_chrom == "V") %>% dplyr::rename(QX_chrom = CGC2_chrom), aes(xmin = -Inf, xmax = Inf, ymin = start/1e6 - 0.001, ymax = end/1e6 + 0.001), fill = '#CC79A7')+
  geom_rect(data = CG_telomeres, aes(xmin = -Inf, xmax = Inf, ymin = start / 1e6, ymax = end / 1e6), fill = '#009E73') +
  geom_segment(aes(x = QXS / 1e6, xend = QXE / 1e6, y = CGS / 1e6, yend = CGE / 1e6), color = 'black', linewidth = 1) +
  geom_vline(xintercept = 0.068, color = '#E69F00') +
  scale_color_gradient(low = "gold", high = "black") +
  facet_wrap(~QX_chrom, scales = "free") +
  theme(
    panel.background = element_blank(),
    legend.position = 'none',
    axis.title = element_text(size = 11, color = 'black'),
    strip.text = element_text(size = 12, color = 'black'),
    panel.grid = element_blank(),
    axis.text = element_text(size = 11, color = 'black'),
    panel.border = element_rect(fill = NA)) +
  labs(y = "CGC2 genome coordinates (Mb)", x = "QX1410 genome coordinates (Mb)") +
  coord_cartesian(xlim = c(0 ,0.3), ylim = c(0,0.5))
chromV_gaps


# After contig reordering and gap closing
gapfree_alignment <- readr::read_tsv("../../processed_data/genome_genome_alignments/CGC2_QX1410.transformed.tsv", col_names = c("QXS","QXE","CGS","CGE","L1","L2","IDY","LENR","LENQ","QX_chrom","CGC2_chrom")) 

# Telomeres after contig reordeding and gap closing
final_tels <- readr::read_tsv("../../processed_data/genomes/CGC2_telomeres_binned_1kb.bed", col_names = c("chrom","start","end","count")) %>%
  dplyr::filter(chrom == "V" & count > 150) %>% # Represents 1 kb segments that are entirely telomeric repeats
  dplyr::filter(start < 150000) %>%
  dplyr::mutate(start = min(start), end = max(end)) %>%
  dplyr::distinct(chrom, start, end)

noGaps <- ggplot(gapfree_alignment %>% dplyr::filter(CGC2_chrom == "V" & QX_chrom == "V")) +
  geom_rect(data = final_tels, aes(xmin = -Inf, xmax = Inf, ymin = start / 1e6, ymax = end / 1e6), fill = '#009E73') +
  geom_segment(aes(x = QXS / 1e6, xend = QXE / 1e6, y = CGS / 1e6, yend = CGE / 1e6), color = 'black', linewidth = 1) +
  geom_hline(yintercept = 0.350 ) +
  geom_vline(xintercept = 0.068, color = '#E69F00') +
  scale_color_gradient(low = "gold", high = "black") +
  facet_wrap(~QX_chrom, scales = "free") +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size = 11, color = 'black'),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 12, color = 'black'),
    legend.position = 'none',
    legend.box.background = element_rect(color = 'black'),
    panel.grid = element_blank(),
    axis.text = element_text(size = 11, color = 'black'),
    panel.border = element_rect(fill = NA)) +
  labs(y = "CGC2 genome coordinates (Mb)", x = "QX1410 genome coordinates (Mb)") +
  coord_cartesian(xlim = c(0 ,0.3), ylim = c(0,0.5))
noGaps


# Concatenating plots together
final_plot <- cowplot::plot_grid(
  chromV_gaps, noGaps,
  labels = c("a","b"),
  nrow = 1)
final_plot

# Saving final plot:
ggsave("../../figures/supplementary/chromV_gaps.png", final_plot, width = 7, height = 5, dpi = 600)

