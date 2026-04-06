library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)

# Chromosome order for plotting
chr_levels <- c("I", "II", "III", "IV", "V", "X")

# Read in nucmer genome-genome alignments
aln <- readr::read_tsv("../../processed_data/genome_genome_alignments/AF16_QX1410.transformed.tsv",
                       col_names = c("QXS","QXE","AFS","AFE","L1","L2","IDY","LENR","LENQ","QX1410","AF16")) 

# Ordering chromosomes and flagging unplaced scaffold alignments
aln_filt <- aln %>%
  dplyr::mutate(QX1410 = factor(QX1410, levels = chr_levels), is_cb = stringr::str_detect(AF16, "^cb")) 

# Filtering to only intra-chromosomal alignments of AF16 nuclear scaffolds
aln_main <- aln_filt %>%
  dplyr::filter(!is_cb) %>%
  dplyr::mutate(class = ifelse(AFS > AFE, "INV", "NOINV")) %>% 
  dplyr::filter(AF16 == QX1410)

# Filtering to only unplaced scaffold alignments
aln_unplaced <- aln_filt %>%
  dplyr::filter(is_cb)

# Plotting genome-genome alignments of AF16 nuclear scaffolds to QX1410
af_qx_main <- ggplot(aln_main) + 
  geom_segment(aes(x = QXS / 1e6, xend = QXE / 1e6, y = AFS / 1e6, yend = AFE / 1e6, color = class), linewidth = 1, alpha = 0.7) +
  scale_color_manual(values = c("INV" = "red", "NOINV" = "black")) +
  facet_wrap(~QX1410, scales = "free", ncol = 3) +
  scale_x_continuous(breaks = seq(5,20,5), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(5,20,5), expand = c(0,0)) +
  labs(x = NULL, y = "AF16 genome coordinates (Mb)") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position = "none",
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12, color = 'black'),
    panel.background = element_blank(),
    axis.text = element_text(size = 10, color = 'black'),
    axis.title.y = element_text(size = 10, color = 'black', vjust = 0))

# Reading in QX1410 chromosome sizes
qx_chrom_sizes <- read.table("../../processed_data/genomes/c_briggsae.QX1410.nanopore.Feb2020.genome.fa.fai", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(qx_chrom_sizes) <- c("QX1410", "chr_len", "offset", "line_bases", "line_width")

# Setting a window size
win_size <- 1e5

# Creating 100 kb windows of QX1410
windows <- qx_chrom_sizes %>%
  dplyr::filter(QX1410 %in% chr_levels) %>%
  dplyr::mutate(QX1410 = factor(QX1410, levels = chr_levels)) %>%
  purrr::pmap_dfr( function(QX1410, chr_len, offset, line_bases, line_width) {
    tibble::tibble(
      QX1410 = QX1410,
      win_start = seq(1, chr_len, by = win_size)) %>%
      dplyr::mutate(
        win_end = pmin(win_start + win_size - 1, chr_len),
        win_mid = (win_start + win_end) / 2,
        win_width = win_end - win_start + 1) })

# Looking at alignments over each 100 kb window and the proportion that are from unplaced AF16 scaffolds
cb_overlap <- windows %>%
  dplyr::left_join(aln_unplaced %>% dplyr::transmute(
    QX1410,
    aln_start = pmin(QXS, QXE),
    aln_end = pmax(QXS, QXE)), by = "QX1410", relationship = "many-to-many") %>%
  dplyr::mutate(overlap_bp = pmax(0, pmin(win_end, aln_end) - pmax(win_start, aln_start) + 1)) %>%
  dplyr::group_by(QX1410, win_start, win_end, win_mid, win_width) %>%
  dplyr::summarise(
    cb_bp = sum(overlap_bp, na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::mutate(prop_cb = pmin(cb_bp / win_width, 1))

# Plotting the proportion of unplaced scaffold alignments for each 100 kb window across the QX1410 genome
af_qx_unplaced <- ggplot(cb_overlap) +
  geom_col(aes(x = win_mid / 1e6, y = prop_cb), width = win_size / 1e6, fill = 'black') +
  facet_wrap(~QX1410, scales = "free", ncol = 3) +
  scale_x_continuous(breaks = seq(5,20,5), expand = c(0,0), name = "QX1410 genome coordinates (Mb)") +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0), labels = scales::percent_format(accuracy = 1),name = "Proportion of sequence\nunplaced in AF16") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position = "none",
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12, color = 'black'),
    panel.background = element_blank(),
    axis.text = element_text(size = 10, color = 'black'),
    axis.title = element_text(size = 10, color = 'black'))

# Concatenating together the two plots
final_af_qx <- cowplot::plot_grid(
  af_qx_main, af_qx_unplaced,
  ncol = 1,
  labels = c("a","b"),
  rel_heights = c(1,0.75),
  align = "v")
final_af_qx

# Saving plot
ggsave("../../figures/AF16_QX1410_dotplot.png", final_af_qx, width = 7, height = 7, dpi = 600)
