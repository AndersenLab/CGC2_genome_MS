library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)

# Chromosome order for plotting
chr_levels <- c("I", "II", "III", "IV", "V", "X")

# Read in nucmer genome-genome alignments
aln <- readr::read_tsv("../../processed_data/genome_genome_alignments/AF16_CGC2.transformed.tsv",
                       col_names = c("CGS","CGE","AFS","AFE","L1","L2","IDY","LENR","LENQ","CGC2","AF16")) 

# Ordering chromosomes and flagging unplaced scaffold alignments
aln_filt <- aln %>%
  dplyr::mutate(CGC2 = factor(CGC2, levels = chr_levels), is_cb = stringr::str_detect(AF16, "^cb")) 

# Filtering to only intra-chromosomal alignments of AF16 nuclear scaffolds
aln_main <- aln_filt %>%
  dplyr::filter(!is_cb) %>%
  dplyr::mutate(class = ifelse(AFS > AFE, "INV", "NOINV")) %>% 
  dplyr::filter(AF16 == CGC2)

# Filtering to only unplaced scaffold alignments
aln_unplaced <- aln_filt %>%
  dplyr::filter(is_cb)

# Plotting genome-genome alignments of AF16 nuclear scaffolds to CGC2
af_CG_main <- ggplot(aln_main) + 
  geom_segment(aes(x = CGS / 1e6, xend = CGE / 1e6, y = AFS / 1e6, yend = AFE / 1e6, color = class), linewidth = 0.75, alpha = 0.7) +
  scale_color_manual(values = c("INV" = "red", "NOINV" = "black")) +
  facet_wrap(~CGC2, scales = "free", ncol = 3) +
  scale_x_continuous(breaks = seq(0,20,5)) +
  scale_y_continuous(breaks = seq(0,20,5)) +
  labs(x = NULL, y = "AF16 genome coordinates (Mb)") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position = "none",
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12, color = 'black'),
    panel.background = element_blank(),
    axis.text = element_text(size = 10, color = 'black'),
    axis.title.y = element_text(size = 10, color = 'black', vjust = 0))

# Reading in CGC2 chromosome sizes
CG_chrom_sizes <- read.table("../../processed_data/genomes/c_briggsae.CGC2.hifi.ONT.HiC.Feb2026.genome.fa.fai", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(CG_chrom_sizes) <- c("CGC2", "chr_len", "offset", "line_bases", "line_width")

# Setting a window size
win_size <- 1e5

# Creating 100 kb windows of CGC2
windows <- CG_chrom_sizes %>%
  dplyr::filter(CGC2 %in% chr_levels) %>%
  dplyr::mutate(CGC2 = factor(CGC2, levels = chr_levels)) %>%
  purrr::pmap_dfr( function(CGC2, chr_len, offset, line_bases, line_width) {
    tibble::tibble(
      CGC2 = CGC2,
      win_start = seq(1, chr_len, by = win_size)) %>%
      dplyr::mutate(
        win_end = pmin(win_start + win_size - 1, chr_len),
        win_mid = (win_start + win_end) / 2,
        win_width = win_end - win_start + 1) })

# Looking at alignments over each 100 kb window and the proportion that are from unplaced AF16 scaffolds
cb_overlap <- windows %>%
  dplyr::left_join(aln_unplaced %>% dplyr::transmute(
    CGC2,
    aln_start = pmin(CGS, CGE),
    aln_end = pmax(CGS, CGE)), by = "CGC2", relationship = "many-to-many") %>%
  dplyr::mutate(overlap_bp = pmax(0, pmin(win_end, aln_end) - pmax(win_start, aln_start) + 1)) %>%
  dplyr::group_by(CGC2, win_start, win_end, win_mid, win_width) %>%
  dplyr::summarise(
    cb_bp = sum(overlap_bp, na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::mutate(prop_cb = pmin(cb_bp / win_width, 1))

# Plotting the proportion of unplaced scaffold alignments for each 100 kb window across the CGC2 genome
af_CG_unplaced <- ggplot(cb_overlap) +
  geom_col(aes(x = win_mid / 1e6, y = prop_cb), width = win_size / 1e6, fill = 'black') +
  facet_wrap(~CGC2, scales = "free", ncol = 3) +
  scale_x_continuous(breaks = seq(5,20,5), expand = c(0,0), name = "CGC2 genome coordinates (Mb)") +
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
final_af_CG <- cowplot::plot_grid(
  af_CG_main, af_CG_unplaced,
  ncol = 1,
  labels = c("a","b"),
  rel_heights = c(1,0.75),
  align = "v")
final_af_CG

# Saving plot
ggsave("../../figures/supplementary/AF16_CGC2_dotplot.png", final_af_CG, width = 7, height = 7, dpi = 600)




















af_CG_main <- ggplot(aln_main) + 
  geom_segment(aes(x = CGS / 1e6, xend = CGE / 1e6, y = AFS / 1e6, yend = AFE / 1e6, color = class), linewidth = 0.75, alpha = 0.7) +
  scale_color_manual(values = c("INV" = "red", "NOINV" = "black")) +
  facet_wrap(~CGC2, scales = "free", ncol = 3) +
  scale_x_continuous(breaks = seq(0,20,5)) +
  scale_y_continuous(breaks = seq(0,20,5)) +
  labs(x = "CGC2 genome coordinates (Mb)", y = "AF16 genome coordinates (Mb)") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position = "none",
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 22, color = 'black'),
    panel.background = element_blank(),
    axis.text = element_text(size = 16, color = 'black'),
    axis.title = element_text(size = 18, color = 'black'))


# Plotting the proportion of unplaced scaffold alignments for each 100 kb window across the CGC2 genome
af_CG_unplaced <- ggplot(cb_overlap) +
  geom_col(aes(x = win_mid / 1e6, y = prop_cb), width = win_size / 1e6, fill = 'black') +
  facet_wrap(~CGC2, scales = "free", ncol = 3) +
  scale_x_continuous(breaks = seq(5,20,5), expand = c(0,0), name = "CGC2 genome coordinates (Mb)") +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0), labels = scales::percent_format(accuracy = 1),name = "Proportion of sequence\nunplaced in AF16") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position = "none",
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 20, color = 'black'),
    panel.background = element_blank(),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.text.x = element_text(size = 16, color = 'black'),
    axis.title = element_text(size = 18, color = 'black'))

# Concatenating together the two plots
cowplot::plot_grid(
  af_CG_main, af_CG_unplaced,
  ncol = 1,
  labels = c("a","b"),
  label_size = 20,
  rel_heights = c(1,0.75),
  align = "v")

