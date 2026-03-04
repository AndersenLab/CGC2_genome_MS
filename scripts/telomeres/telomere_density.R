library(ggplot2)
library(readr)
library(dplyr)

# Desired chromosome ordering for plotting
chrom_order <- c("I","II","III","IV","V","X")

# Laod in CGC2 telomere counts
tel_ct_CGC2 <- readr::read_tsv("../../processed_data/genomes/CGC2_telomeres_binned_1kb.bed", col_names = c("chrom","start","end","count")) %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order), mid_Mb = (start + end) / 2 / 1e6) %>%
  dplyr::filter(!grepl(">AC186293.1","chrom"))

# Laod in AF16 telomere counts
tel_ct_AF16 <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_HiC/telomeres/af_telomeres_binned_1kb.bed", col_names = c("chrom","start","end","count")) %>%
  dplyr::filter(chrom != "MtDNA" & !grepl("cb25", chrom)) %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order), mid_Mb = (start + end) / 2 / 1e6)

# Laod in QX1410 telomere counts
tel_ct_QX1410 <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_HiC/telomeres/QX_telomeres_binned_1kb.bed", col_names = c("chrom","start","end","count")) %>%
  dplyr::filter(chrom != "MtDNA") %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order), mid_Mb = (start + end) / 2 / 1e6)

# Laod in VX34 telomere counts
tel_ct_VX34 <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_HiC/telomeres/VX_telomeres_binned_1kb.bed", col_names = c("chrom","start","end","count")) %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order), mid_Mb = (start + end) / 2 / 1e6)

# Plotting telomere density
tel_count_plt <- cowplot::plot_grid(
  ggplot(data = tel_ct_CGC2) +
    geom_col(aes(x = mid_Mb, y = count), color = 'black', width = 0.0001) +
    facet_wrap(~chrom, nrow = 1, scales = "free") +
    theme_bw() +
    theme(
      panel.border = element_rect(color = 'black', fill = NA),
      strip.text = element_text(size = 11, color = 'black'),
      axis.text = element_text(size = 10, color = 'black'),
      axis.title = element_text(size = 11, color = 'black')
    ) +
    labs(x = "CGC2 genome position (Mb)", y = "Counts of TTAGGC per kb") +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    coord_cartesian(ylim = c(0,180)),
  
  ggplot(data = tel_ct_AF16) +
    geom_col(aes(x = mid_Mb, y = count), color = 'blue', width = 0.0001) +
    facet_wrap(~chrom, nrow = 1, scales = "free_x") +
    theme_bw() +
    theme(
      panel.border = element_rect(color = 'black', fill = NA),
      strip.text = element_text(size = 11, color = 'black'),
      axis.text = element_text(size = 10, color = 'black'),
      axis.title = element_text(size = 11, color = 'black')
    ) +
    labs(x = "AF16 genome position (Mb)", y = "Counts of TTAGGC per kb") +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    coord_cartesian(ylim = c(0,180)),
  
  ggplot(data = tel_ct_QX1410) +
    geom_col(aes(x = mid_Mb, y = count), color = 'forestgreen', width = 0.0001) +
    facet_wrap(~chrom, nrow = 1, scales = "free_x") +
    theme_bw() +
    theme(
      panel.border = element_rect(color = 'black', fill = NA),
      strip.text = element_text(size = 11, color = 'black'),
      axis.text = element_text(size = 10, color = 'black'),
      axis.title = element_text(size = 11, color = 'black')
    ) +
    labs(x = "QX1410 genome position (Mb)", y = "Counts of TTAGGC per kb") +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    coord_cartesian(ylim = c(0,180)),
  
  ggplot(data = tel_ct_VX34) +
    geom_col(aes(x = mid_Mb, y = count), color = 'purple', width = 0.0001) +
    facet_wrap(~chrom, nrow = 1, scales = "free_x") +
    theme_bw() +
    theme(
      panel.border = element_rect(color = 'black', fill = NA),
      strip.text = element_text(size = 11, color = 'black'),
      axis.text = element_text(size = 10, color = 'black'),
      axis.title = element_text(size = 11, color = 'black')
    ) +
    labs(x = "VX34 genome position (Mb)", y = "Counts of TTAGGC per kb") +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    coord_cartesian(ylim = c(0,180)),
  
  
  nrow = 4)

tel_count_plt

# Save plot
ggsave("../../figures/supplementary/telomere_density.png", tel_count_plt, width = 6, height = 7, dpi = 600)

