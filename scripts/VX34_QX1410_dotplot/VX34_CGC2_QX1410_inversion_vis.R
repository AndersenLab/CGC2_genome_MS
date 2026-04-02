library(ggplot2)
library(readr)
library(dplyr)
                  

# Loading in nucmer alignment data of VX34 to QX1410
vx_qx <- readr::read_tsv("../../processed_data/genome_genome_alignments/VX34_QX1410.transformed.tsv", col_names = c("QXS","QXE","VXS","VXE","L1","L2","IDY","LENR","LENQ","QX_chrom","VX34_chrom")) %>%
  dplyr::mutate(inv = ifelse(VXE < VXS, "yes", "no")) %>%
  dplyr::filter(QX_chrom == "V" & VX34_chrom == "V")

# Plotting alignment of inversion on chromosome V in VX34
vx_qx_plt <- ggplot(vx_qx) +
  geom_segment(aes(x = QXS/1e6, xend = QXE/1e6, y = VXS/1e6, yend = VXE/1e6, color = inv), linewidth = 1) +
  facet_grid(~QX_chrom) +
  scale_color_manual(values = c(no = "black", yes = "red")) +
  geom_vline(xintercept = 14.934651, linetype = "dashed", linewidth = 0.75, color = "blue") +
  geom_vline(xintercept = 15.559394, linetype = "dashed", linewidth = 0.75, color = "blue") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    strip.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  labs(y = "VX34 genome coordinates (Mb)", x = "QX1410 genome coordinates (Mb)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(13.5,17), ylim = c(13.5,17))
vx_qx_plt


# Loading in nucmer alignment data of CGC2 to QX1410
cg_qx <- readr::read_tsv("../../processed_data/genome_genome_alignments/CGC2_QX1410.transformed.tsv", col_names = c("QXS","QXE","CGS","CGE","L1","L2","IDY","LENR","LENQ","QX_chrom","CGC2_chrom")) %>%
  dplyr::mutate(inv = ifelse(CGE < CGS, "yes", "no")) %>%
  dplyr::filter(QX_chrom == "V" & CGC2_chrom == "V")

filtered <- cg_qx %>% dplyr::filter(inv == "yes" & L2 > 1000)

# Plotting alignment of inversion on chromosome V in CGC2
cg_qx_plt <- ggplot(cg_qx) +
  geom_segment(aes(x = QXS/1e6, xend = QXE/1e6, y = CGS/1e6, yend = CGE/1e6, color = inv), linewidth = 1) +
  facet_grid(~QX_chrom) +
  scale_color_manual(values = c(no = "black", yes = "red")) +
  geom_vline(xintercept = 14.934651, linetype = "dashed", linewidth = 0.75, color = "blue") +
  geom_vline(xintercept = 15.559394, linetype = "dashed", linewidth = 0.75, color = "blue") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    strip.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  labs(y = "CGC2 genome coordinates (Mb)", x = "QX1410 genome coordinates (Mb)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(13.5,17), ylim = c(13.5,17))
cg_qx_plt

# Creating final plot:
final_plot <- cowplot::plot_grid(
  vx_qx_plt, cg_qx_plt,
  nrow = 1,
  align = "h",
  axis = "tb",
  labels = c("a", "b")
)
final_plot


# Saving plot
ggsave("../../figures/supplementary/VX34_CGC2_QX1410_inversion.png", final_plot, width = 7, height = 3.5, dpi = 600)

