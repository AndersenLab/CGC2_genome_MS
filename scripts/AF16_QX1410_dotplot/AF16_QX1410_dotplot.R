library(ggplot2)
library(dplyr)
library(readr)
library(ComplexHeatmap)
library(ggplotify)

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
  labels = c('a'), 
  # label_size = 12,
  ncol = 1,
  rel_heights = c(0.9, 0.1))
final_af_qx


# Loading in heatmap for panels B and C
gtcheck <- readr::read_tsv("../../processed_data/AF16_ancestry/gtcheck.allpairs.SNV.nohead.txt",col_names = c("DCv2","S1","S2","discordance","eval","total_sites","matches")) %>%
  dplyr::mutate(concordance=matches/total_sites)

heatmap_data <- gtcheck %>% dplyr::filter(S1!="QX1410" & S2!="QX1410")

samples <- sort(unique(c(heatmap_data$S1, heatmap_data$S2)))

pair_grid <- tidyr::expand_grid(
  S1 = samples,
  S2 = samples) %>%
  dplyr::mutate(concordance = dplyr::case_when(S1 == S2 ~ 1, TRUE ~ NA_real_))

hm_full <- pair_grid %>%
  dplyr::left_join(heatmap_data %>% dplyr::select(S1, S2, concordance), by = c("S1", "S2"), suffix = c("", ".from_data")) %>%
  dplyr::left_join(heatmap_data %>% dplyr::select(S1, S2, concordance) %>% dplyr::rename(S1_rev = S2, S2_rev = S1, concordance_rev = concordance), by = c("S1" = "S1_rev", "S2" = "S2_rev")) %>%
  dplyr::mutate(concordance = dplyr::coalesce(concordance.from_data, concordance_rev, concordance)) %>%
  dplyr::select(S1, S2, concordance)

mat <- hm_full %>%
  tidyr::pivot_wider(names_from = S1, values_from = concordance) %>%
  tibble::column_to_rownames("S2") %>%
  as.matrix()

HM1 <- ComplexHeatmap::Heatmap(
  mat,
  name = "Proportion of\nidentical SNV\nalleles",
  col = viridisLite::viridis(100),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  rect_gp = grid::gpar(col = NA),
  column_names_gp = grid::gpar(fontsize = 8),
  row_names_gp = grid::gpar(fontsize = 8),
  row_names_side = "left",
  column_names_rot = 45,
  column_names_side = "top",
  heatmap_legend_param = list(
    title = "Proportion of\nidentical SNV\nalleles",
    title_gp = gpar(fontface = "plain", fontsize = 8), 
    labels_gp = gpar(fontface = "plain", fontsize = 8) 
  ))
HM1

HM1_grob <- grid::grid.grabExpr(
  ComplexHeatmap::draw(HM1,padding = grid::unit(c(5, 10, 5, 5), "mm")))

HM1_gg <- ggplotify::as.ggplot(HM1_grob)

target <- "QX1410"

hm2_df <- gtcheck %>%
  dplyr::filter(S1 == target | S2 == target) %>%
  dplyr::mutate(sample = ifelse(S1 == target, S2, S1)) %>%
  dplyr::select(sample, concordance) %>%
  dplyr::arrange(sample)

mat2 <- matrix(
  hm2_df$concordance,
  nrow = 1,
  dimnames = list(target, hm2_df$sample))

HM2 <- ComplexHeatmap::Heatmap(
  t(mat2),
  name = "Proportion of\nidentical SNV\nalleles",
  col = viridisLite::viridis(100),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_gp = grid::gpar(fontsize = 8),
  row_names_gp = grid::gpar(fontsize = 8),
  rect_gp = grid::gpar(col = NA),
  column_names_side = "top",
  column_names_rot = 45,
  row_names_side = "left",
  heatmap_legend_param = list(
    title = "Proportion of\nidentical SNV\nalleles",
    title_gp = gpar(fontface = "plain", fontsize = 8), 
    labels_gp = gpar(fontface = "plain", fontsize = 8) 
  ))

HM2_grob <- grid::grid.grabExpr(
  ComplexHeatmap::draw(HM2,padding = grid::unit(c(5, 10, 6, 5), "mm")))

HM2_gg <- ggplotify::as.ggplot(HM2_grob)

heatmaps <- cowplot::plot_grid(HM2_gg,HM1_gg,
                           nrow=1,
                           align="h",
                           axis = "tb",
                           rel_widths = c(1,2.2),
                           labels = c("b","c"))
heatmaps

final_plot <- cowplot::plot_grid(
  final_af_qx, heatmaps,
  nrow = 2,
  rel_heights = c(0.7,0.3))
final_plot



# Saving plot
ggsave("../../figures/AF16_QX1410_dotplot/AF16_QX1410_dotplot.png", final_plot, width = 7, height = 10, dpi = 600, bg = 'white')

