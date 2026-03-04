library(ggplot2)
library(readr)
library(dplyr)
                  
# Chromosome ordering for plot
chrom_levels <- c("I","II","III","IV","V","X")
inv_chrom_levels <- c("X","V","IV","III","II","I")

# Loading in nucmer alignment data of AF16 to QX1410
cg_qx <- readr::read_tsv("../../processed_data/genome_genome_alignments/CGC2_QX1410.transformed.tsv", col_names = c("QXS","QXE","CGS","CGE","L1","L2","IDY","LENR","LENQ","QX_chrom","CGC2_chrom")) %>%
  dplyr::filter(!grepl(">AC186293.1", CGC2_chrom)) %>% # removing MtDNA
  dplyr::mutate(
    inv = ifelse(CGE < CGS, "yes", "no"),
    # Cap alignments boundaries to the max and min of each chromosome to remove spurious inter-chromosomal alignments 
    QXS_c = pmin(pmax(QXS, 0), LENR),
    QXE_c = pmin(pmax(QXE, 0), LENR),
    CGS_c = pmin(pmax(CGS, 0), LENQ),
    CGE_c = pmin(pmax(CGE, 0), LENQ),
    QX_chrom = factor(QX_chrom, levels = chrom_levels),
    CGC2_chrom = factor(CGC2_chrom, levels = inv_chrom_levels))

# Plotting alignment dotplot faceted by each strains chromosomes to visualzie co-linearity
cg_qx_plt <- ggplot(cg_qx) +
  geom_blank(data = cg_qx, aes(x = 0, y = 0)) +
  geom_blank(data = cg_qx, aes(x = LENR/1e6, y = 0)) +
  geom_segment(aes(x = QXS_c/1e6, xend = QXE_c/1e6, y = CGS_c/1e6, yend = CGE_c/1e6, color = inv), linewidth = 1) +
  facet_grid(CGC2_chrom ~ QX_chrom, scales = "free", space = "free") +
  scale_color_manual(values = c(no = "black", yes = "red")) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    strip.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    plot.margin = margin(b = 0),
    legend.position = "none"
  ) +
  labs(y = "CGC2 genome coordinates (Mb)", x = "QX1410 genome coordinates (Mb)") +
  scale_x_continuous(expand = c(0,0), breaks = seq(5,20,5)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(5,20,5))
cg_qx_plt

# Saving plot
ggsave("../../figures/CGC2_QX1410_dotplot/CGC2_QX1410_dotplot.png", cg_qx_plt, width = 7.5, height = 7.5, dpi = 500)