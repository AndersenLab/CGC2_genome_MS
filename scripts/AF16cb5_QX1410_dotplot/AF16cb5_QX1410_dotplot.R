library(ggplot2)
library(readr)
library(dplyr)
                  
# Chromosome ordering for plot
chrom_levels <- c("I","II","III","IV","V","X")

# Loading in nucmer alignment data of AF16 CB5 to QX1410
afcb_qx <- readr::read_tsv("../../processed_data/genome_genome_alignments/AF16cb5_QX1410.transformed.tsv", 
                         col_names = c("QXS","QXE","AFS","AFE","L1","L2","IDY","LENR","LENQ","QX1410","AF16cb")) %>%
  dplyr::mutate(class = ifelse(AFS > AFE, "INV", "NOINV")) %>%
  dplyr::filter(QX1410 == AF16cb)
  
# Plotting alignment dotplot faceted by each strains chromosomes to visualzie co-linearity
afcb_qx_main <- ggplot(afcb_qx) + 
  geom_segment(aes(x = QXS / 1e6, xend = QXE / 1e6, y = AFS / 1e6, yend = AFE / 1e6, color = class), linewidth = 0.75, alpha = 0.7) +
  scale_color_manual(values = c("INV" = "red", "NOINV" = "black")) +
  facet_wrap(~QX1410, scales = "free", ncol = 3) +
  scale_x_continuous(breaks = seq(0, 20, 5)) +
  scale_y_continuous(breaks = seq(0, 20, 5)) +
  labs(x = "QX1410 genome coordinates (Mb)", y = "AF16 cb5 genome coordinates (Mb)") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position = "none",
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12, color = 'black'),
    panel.background = element_blank(),
    axis.text = element_text(size = 10, color = 'black'),
    axis.title = element_text(size = 10, color = 'black'))

# Saving plot
ggsave("../../figures/supplementary/AF16cb5_QX1410_dotplot.png", afcb_qx_main, width = 7, height = 5, dpi = 600)