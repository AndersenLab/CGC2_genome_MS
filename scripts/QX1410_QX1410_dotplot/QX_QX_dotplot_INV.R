library(ggplot2)
library(readr)
library(dplyr)

# Loading in nucmer alignment data of QX1410 to QX1410
qxqx <- readr::read_tsv("../../processed_data/genome_genome_alignments/QX_QX.transformed.tsv", col_names = c("QXS1","QXE1","QXS2","QXE2","L1","L2","IDY","LENR","LENQ","QX_chrom","QX_chrom2")) %>%
  dplyr::mutate(inv = ifelse(QXE2 < QXS2, "yes", "no")) %>%
  dplyr::filter(QX_chrom == "V" & QX_chrom2 == "V")

# Plotting alignment of inversion locus in QX1410
qx_qx_plt_main <- ggplot(qxqx) +
  geom_segment(aes(x = QXS1 / 1e6, xend = QXE1 / 1e6, y = QXS2 / 1e6, yend = QXE2 / 1e6, color = inv), linewidth = 0.75) +
  facet_grid(~QX_chrom2) +
  scale_color_manual(values = c(no = "black", yes = "red")) +
  geom_vline(xintercept = 14.934651, linetype = "dashed", linewidth = 0.75, color = "blue") +
  geom_vline(xintercept = 15.559394, linetype = "dashed", linewidth = 0.75, color = "blue") +
  geom_hline(yintercept = 14.934651, linetype = "dashed", linewidth = 0.75, color = "blue") +
  geom_hline(yintercept = 15.559394, linetype = "dashed", linewidth = 0.75, color = "blue") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    strip.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 11, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    plot.margin = margin(t = 15, r = 15, l = 15, b = 15),
    legend.position = "none"
  ) +
  labs(y = "QX1410 genome coordinates (Mb)", x = "QX1410 genome coordinates (Mb)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(14.5,16), ylim = c(14.5,16))

# Save the main plot
ggsave("../../figures/supplementary/qx_qx_INV_visualization.png", qx_qx_plt_main, height = 5, width = 5, dpi = 600)

# Looking at multi-alignment at left inversion breakpoint
left_inset <- ggplot(qxqx) +
  geom_segment(aes(x = QXS1 / 1e6, xend = QXE1 / 1e6, y = QXS2 / 1e6, yend = QXE2 / 1e6, color = inv), linewidth = 0.4) +
  facet_grid(~QX_chrom2) +
  scale_color_manual(values = c(no = "black", yes = "red")) +
  geom_vline(xintercept = 14.934651, linetype = "dashed", linewidth = 0.3, color = "blue") +
  geom_vline(xintercept = 15.559394, linetype = "dashed", linewidth = 0.3, color = "blue") +
  geom_hline(yintercept = 14.934651, linetype = "dashed", linewidth = 0.3, color = "blue") +
  geom_hline(yintercept = 15.559394, linetype = "dashed", linewidth = 0.3, color = "blue") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    strip.text = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) +
  labs(y = "QX1410 genome coordinates (Mb)", x = "QX1410 genome coordinates (Mb)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(14.9245,14.9445), ylim = c(15.55,15.57))

# Save the left inset
ggsave("../../figures/supplementary/qx_qx_INV_visualization_left_inset.png", left_inset, height = 2, width = 2, dpi = 600)

# Looking at multi-alignment at right inversion breakpoint
right_inset <- ggplot(qxqx) +
  geom_segment(aes(x = QXS1 / 1e6, xend = QXE1 / 1e6, y = QXS2 / 1e6, yend = QXE2 / 1e6, color = inv), linewidth = 0.4) +
  facet_grid(~QX_chrom2) +
  scale_color_manual(values = c(no = "black", yes = "red")) +
  geom_vline(xintercept = 14.934651, linetype = "dashed", linewidth = 0.3, color = "blue") +
  geom_vline(xintercept = 15.559394, linetype = "dashed", linewidth = 0.3, color = "blue") +
  geom_hline(yintercept = 14.934651, linetype = "dashed", linewidth = 0.3, color = "blue") +
  geom_hline(yintercept = 15.559394, linetype = "dashed", linewidth = 0.3, color = "blue") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    strip.text = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) +
  labs(y = "QX1410 genome coordinates (Mb)", x = "QX1410 genome coordinates (Mb)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(15.5495,15.5695), ylim = c(14.9245,14.9445))

# Save the right inset
ggsave("../../figures/supplementary/qx_qx_INV_visualization_right_inset.png", right_inset, height = 2, width = 2, dpi = 600)


