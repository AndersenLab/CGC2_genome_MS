library(ggplot2)
library(readr)
library(dplyr)

telomere_count <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_HiC/telomeres/CGC2_telomeres_binned_1kb.bed", col_names = c("chrom","start","end","count"))

chrom_order <- c("I","II","III","IV","V","X")

tel_plt <- telomere_count %>% dplyr::mutate(chrom = factor(chrom, levels = chrom_order), mid_Mb = (start + end) / 2 / 1e6)

# ggplot(data = tel_plt) +
#   geom_col(aes(x = mid_Mb, y = count), color = 'black', width = 0.0001) +
#   facet_wrap(~chrom, nrow = 1, scales = "free") +
#   theme_bw() +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     strip.text = element_text(size = 14, color = 'black'),
#     axis.text = element_text(size = 12, color = 'black'),
#     axis.title = element_text(size = 14, color = 'black')
#   ) +
#   labs(x = "CGC2 genome position (Mb)", y = "Counts of TTAGGC per kb") +
#   scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
#   coord_cartesian(ylim = c(0,180))


vx_telomere_count <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_HiC/telomeres/VX_telomeres_binned_1kb.bed", col_names = c("chrom","start","end","count"))

tel_plt_vx <- vx_telomere_count %>% dplyr::mutate(chrom = factor(chrom, levels = chrom_order), mid_Mb = (start + end) / 2 / 1e6)

# ggplot(data = tel_plt_vx) +
#   geom_col(aes(x = mid_Mb, y = count), color = 'red', width = 0.0001) +
#   facet_wrap(~chrom, nrow = 1, scales = "free_x") +
#   theme_bw() +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     strip.text = element_text(size = 14, color = 'black'),
#     axis.text = element_text(size = 12, color = 'black'),
#     axis.title = element_text(size = 14, color = 'black')
#   ) +
#   labs(x = "VX34 genome position (Mb)", y = "Counts of TTAGGC per kb") +
#   scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
#   coord_cartesian(ylim = c(0,180))



qx_telomere_count <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_HiC/telomeres/QX_telomeres_binned_1kb.bed", col_names = c("chrom","start","end","count")) %>%
  dplyr::filter(chrom != "MtDNA")

tel_plt_qx <- qx_telomere_count %>% dplyr::mutate(chrom = factor(chrom, levels = chrom_order), mid_Mb = (start + end) / 2 / 1e6)

# ggplot(data = tel_plt_qx) +
#   geom_col(aes(x = mid_Mb, y = count), color = 'forestgreen', width = 0.0001) +
#   facet_wrap(~chrom, nrow = 1, scales = "free_x") +
#   theme_bw() +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     strip.text = element_text(size = 14, color = 'black'),
#     axis.text = element_text(size = 12, color = 'black'),
#     axis.title = element_text(size = 14, color = 'black')
#   ) +
#   labs(x = "QX1410 genome position (Mb)", y = "Counts of TTAGGC per kb") +
#   scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
#   coord_cartesian(ylim = c(0,180))



af_telomere_count <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_HiC/telomeres/af_telomeres_binned_1kb.bed", col_names = c("chrom","start","end","count")) %>%
  dplyr::filter(chrom != "MtDNA" & !grepl("cb25", chrom))

tel_plt_af <- af_telomere_count %>% dplyr::mutate(chrom = factor(chrom, levels = chrom_order), mid_Mb = (start + end) / 2 / 1e6)

# ggplot(data = tel_plt_af) +
#   geom_col(aes(x = mid_Mb, y = count), color = 'seagreen', width = 0.0001) +
#   facet_wrap(~chrom, nrow = 1, scales = "free_x") +
#   theme_bw() +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     strip.text = element_text(size = 14, color = 'black'),
#     axis.text = element_text(size = 12, color = 'black'),
#     axis.title = element_text(size = 14, color = 'black')
#   ) +
#   labs(x = "AF16 genome position (Mb)", y = "Counts of TTAGGC per kb") +
#   scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
#   coord_cartesian(ylim = c(0,180))








cowplot::plot_grid(
  ggplot(data = tel_plt) +
    geom_col(aes(x = mid_Mb, y = count), color = 'black', width = 0.0001) +
    facet_wrap(~chrom, nrow = 1, scales = "free") +
    theme_bw() +
    theme(
      panel.border = element_rect(color = 'black', fill = NA),
      strip.text = element_text(size = 14, color = 'black'),
      axis.text = element_text(size = 12, color = 'black'),
      axis.title = element_text(size = 14, color = 'black')
    ) +
    labs(x = "CGC2 genome position (Mb)", y = "Counts of TTAGGC per kb") +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    coord_cartesian(ylim = c(0,180)),
  
  ggplot(data = tel_plt_af) +
    geom_col(aes(x = mid_Mb, y = count), color = 'blue', width = 0.0001) +
    facet_wrap(~chrom, nrow = 1, scales = "free_x") +
    theme_bw() +
    theme(
      panel.border = element_rect(color = 'black', fill = NA),
      strip.text = element_text(size = 14, color = 'black'),
      axis.text = element_text(size = 12, color = 'black'),
      axis.title = element_text(size = 14, color = 'black')
    ) +
    labs(x = "AF16 genome position (Mb)", y = "Counts of TTAGGC per kb") +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    coord_cartesian(ylim = c(0,180)),
  
  ggplot(data = tel_plt_qx) +
    geom_col(aes(x = mid_Mb, y = count), color = 'forestgreen', width = 0.0001) +
    facet_wrap(~chrom, nrow = 1, scales = "free_x") +
    theme_bw() +
    theme(
      panel.border = element_rect(color = 'black', fill = NA),
      strip.text = element_text(size = 14, color = 'black'),
      axis.text = element_text(size = 12, color = 'black'),
      axis.title = element_text(size = 14, color = 'black')
    ) +
    labs(x = "QX1410 genome position (Mb)", y = "Counts of TTAGGC per kb") +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    coord_cartesian(ylim = c(0,180)),
  
  ggplot(data = tel_plt_vx) +
    geom_col(aes(x = mid_Mb, y = count), color = 'purple', width = 0.0001) +
    facet_wrap(~chrom, nrow = 1, scales = "free_x") +
    theme_bw() +
    theme(
      panel.border = element_rect(color = 'black', fill = NA),
      strip.text = element_text(size = 14, color = 'black'),
      axis.text = element_text(size = 12, color = 'black'),
      axis.title = element_text(size = 14, color = 'black')
    ) +
    labs(x = "VX34 genome position (Mb)", y = "Counts of TTAGGC per kb") +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    coord_cartesian(ylim = c(0,180)),
  
  
  nrow = 4
)
