library(ggplot2)
library(dplyr)
library(readr)

# Chromosome ordering for plot
chrom_levels <- c("I","II","III","IV","V","X")
inv_chrom_levels <- c("X","V","IV","III","II","I")

# Loading in nucmer alignment data of AF16 to QX1410
af_qx <- readr::read_tsv("../../processed_data/genome_genome_alignments/AF16_QX1410.transformed.tsv", col_names = c("QXS","QXE","AFS","AFE","L1","L2","IDY","LENR","LENQ","QX_chrom","AF_chrom")) %>%
  dplyr::filter(!grepl("cb25", AF_chrom)) %>%
  dplyr::mutate(
    inv = ifelse(AFE < AFS, "yes", "no"),
    # Cap alignments boundaries to the max and min of each chromosome to remove spurious inter-chromosomal alignments 
    QXS_c = pmin(pmax(QXS, 0), LENR),
    QXE_c = pmin(pmax(QXE, 0), LENR),
    AFS_c = pmin(pmax(AFS, 0), LENQ),
    AFE_c = pmin(pmax(AFE, 0), LENQ),
    QX_chrom = factor(QX_chrom, levels = chrom_levels),
    AF_chrom = factor(AF_chrom, levels = inv_chrom_levels))

# Plotting alignment dotplot faceted by each strains chromosomes to visualzie co-linearity
af_qx_plt <- ggplot(af_qx) +
  geom_blank(data = af_qx, aes(x = 0, y = 0)) +
  geom_blank(data = af_qx, aes(x = LENR/1e6, y = 0)) +
  geom_segment(aes(x = QXS_c/1e6, xend = QXE_c/1e6, y = AFS_c/1e6, yend = AFE_c/1e6, color = inv), linewidth = 1) +
  facet_grid(AF_chrom ~ QX_chrom, scales = "free", space = "free") +
  scale_color_manual(values = c(no = "black", yes = "red")) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    strip.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(b = 0),
    legend.position = "none"
  ) +
  labs(y = "AF16 genome coordinates (Mb)") +
  scale_x_continuous(expand = c(0,0), breaks = seq(5,20,5)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(5,20,5))
af_qx_plt

# Filtering alignments to only unplaced AF16 contigs
af_unplaced <- readr::read_tsv("../../processed_data/genome_genome_alignments/AF16_QX1410.transformed.tsv", col_names = c("QXS","QXE","AFS","AFE","L1","L2","IDY","LENR","LENQ","QX_chrom","AF_chrom")) %>%
  dplyr::filter(grepl("cb25", AF_chrom))

# Plotting where unplaced AF16 contigs align to the QX1410 genome
af_unplaced_plt <- ggplot(af_unplaced) +
  geom_blank(data = af_qx, aes(x = 0, y = 0)) +
  geom_blank(data = af_qx, aes(x = LENR/1e6, y = 0)) +
  geom_segment(aes(x = QXS / 1e6, xend = QXE / 1e6, y = AFS / 1e6, yend = AFE / 1e6), color = "blue", linewidth = 1) +
  facet_grid(~QX_chrom, scales="free_x", space="free_x", drop = FALSE) +
  theme(
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.title = element_text(size = 14, color = 'black'),
    axis.text = element_text(size = 10, color = 'black'),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 0),
    panel.border = element_rect(fill = NA)) +
  labs(x = "QX1410 genome coordinates (Mb)") +
  scale_x_continuous(expand = c(0,0), breaks = seq(5,20,5)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,0.08,0.04))
af_unplaced_plt

# Ensuring each x-axis facet (QX1410 chromosome) is aligned correctly
aligned <- cowplot::align_plots(af_qx_plt, af_unplaced_plt, align = "v", axis = "lr")

# Concatenating plots together
final_af_qx <- cowplot::plot_grid(
  aligned[[1]], aligned[[2]],
  ncol = 1,
  rel_heights = c(0.9, 0.1))
final_af_qx


# Loading in heatmap for panels B and C
load("../../processed_data/AF16_ancestry/HEATMAP1.Rda")
load("../../processed_data/AF16_ancestry/HEATMAP2.Rda")

hm1 <- HM1 +
  ggplot2::theme(axis.text.y = element_text(size = 20, color = 'black'))
hm1

hm2 <- HM2 +
  ggplot2::theme()
hm2

HM2



dotplot <- cowplot::plot_grid(
  aligned[[1]], aligned[[2]],
  ncol = 1,
  rel_heights = c(0.9, 0.1))
dotplot

heatmaps <- cowplot::plot_grid()
heatmaps

final_plot <- cowplot::plot_grid()
final_plot



# Saving plot
ggsave("../../figures/AF16_QX1410_dotplot/AF16_QX1410_dotplot.png", final_af_qx, width = 7, height = 10, dpi = 600)

