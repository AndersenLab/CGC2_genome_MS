library(ggplot2)
library(readr)
library(dplyr)
library(cowplot)

# Desired chromosome ordering for plotting
chrom_order <- c("I","II","III","IV","V","X")

# Laod in CGC2 telomere counts
tel_ct_CGC2 <- readr::read_tsv("../../processed_data/genomes/CGC2_telomeres_binned_1kb.bed", col_names = c("chrom","start","end","count")) %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order), mid_Mb = (start + end) / 2 / 1e6) %>%
  dplyr::filter(!grepl(">AC186293.1","chrom"))

# Load in AF16 cb5 telomere counts
tel_ct_AF16cb5 <- readr::read_tsv("../../processed_data/genomes/AF16_cb5_telomeres_binned_1kb.bed", col_names = c("chrom","start","end","count")) %>%
  dplyr::filter(chrom != "MtDNA" & !grepl("cb25", chrom)) %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order), mid_Mb = (start + end) / 2 / 1e6)

# Laod in AF16 cb4 telomere counts
tel_ct_AF16 <- readr::read_tsv("../../processed_data/genomes/AF16_telomeres_binned_1kb.bed", col_names = c("chrom","start","end","count")) %>%
  dplyr::filter(chrom != "MtDNA" & !grepl("cb25", chrom)) %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order), mid_Mb = (start + end) / 2 / 1e6)

# Laod in QX1410 telomere counts
tel_ct_QX1410 <- readr::read_tsv("../../processed_data/genomes/QX1410_telomeres_binned_1kb.bed", col_names = c("chrom","start","end","count")) %>%
  dplyr::filter(chrom != "MtDNA") %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order), mid_Mb = (start + end) / 2 / 1e6)

# Laod in VX34 telomere counts
tel_ct_VX34 <- readr::read_tsv("../../processed_data/genomes/VX34_telomeres_binned_1kb.bed", col_names = c("chrom","start","end","count")) %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order), mid_Mb = (start + end) / 2 / 1e6)


# Plotting telomere density
tel_count_plt <- 
  cowplot::plot_grid(
      ggplot(data = tel_ct_CGC2) +
        geom_col(aes(x = mid_Mb, y = count), color = '#7570B3', width = 0.0001) +
        facet_wrap(~chrom, nrow = 1, scales = "free_x") +
        theme_bw() +
        theme(
          panel.border = element_rect(color = 'black', fill = NA),
          strip.text = element_text(size = 11, color = 'black'),
          axis.text = element_text(size = 9, color = 'black'),
          axis.title = element_text(size = 11, color = 'black'),
          plot.margin = margin(l = 20, r = 2, t = 2, b = 2),
          axis.title.y = element_blank()
        ) +
        labs(x = "CGC2 genome coordinates (Mb)") +
        scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
        coord_cartesian(ylim = c(0,180)),
      
      ggplot(data = tel_ct_AF16cb5) +
        geom_col(aes(x = mid_Mb, y = count), color = '#1F78B4', width = 0.0001) +
        facet_wrap(~chrom, nrow = 1, scales = "free_x") +
        theme_bw() +
        theme(
          panel.border = element_rect(color = 'black', fill = NA),
          strip.text = element_text(size = 11, color = 'black'),
          axis.text = element_text(size = 9, color = 'black'),
          axis.title = element_text(size = 11, color = 'black'),
          plot.margin = margin(l = 20, r = 2, t = 2, b = 2),
          axis.title.y = element_blank()
        ) +
        labs(x = "AF16 cb5 genome coordinates (Mb)") +
        scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
        coord_cartesian(ylim = c(0,180)),
      
      ggplot(data = tel_ct_AF16) +
        geom_col(aes(x = mid_Mb, y = count), color = '#E6B3B3', width = 0.0001) +
        facet_wrap(~chrom, nrow = 1, scales = "free_x") +
        theme_bw() +
        theme(
          panel.border = element_rect(color = 'black', fill = NA),
          strip.text = element_text(size = 11, color = 'black'),
          axis.text = element_text(size = 9, color = 'black'),
          axis.title = element_text(size = 11, color = 'black'),
          plot.margin = margin(l = 20, r = 2, t = 2, b = 2),
          axis.title.y = element_blank()
        ) +
        labs(x = "AF16 cb4 genome coordinates (Mb)") +
        scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
        coord_cartesian(ylim = c(0,180)),
      
      ggplot(data = tel_ct_QX1410) +
        geom_col(aes(x = mid_Mb, y = count), color = '#53886C', width = 0.0001) +
        facet_wrap(~chrom, nrow = 1, scales = "free_x") +
        theme_bw() +
        theme(
          panel.border = element_rect(color = 'black', fill = NA),
          strip.text = element_text(size = 11, color = 'black'),
          axis.text = element_text(size = 9, color = 'black'),
          axis.title = element_text(size = 11, color = 'black'),
          plot.margin = margin(l = 20, r = 2, t = 2, b = 2),
          axis.title.y = element_blank()
        ) +
        labs(x = "QX1410 genome coordinates (Mb)") +
        scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
        coord_cartesian(ylim = c(0,180)),
      
      ggplot(data = tel_ct_VX34) +
        geom_col(aes(x = mid_Mb, y = count), color = '#CC799D', width = 0.0001) +
        facet_wrap(~chrom, nrow = 1, scales = "free_x") +
        theme_bw() +
        theme(
          panel.border = element_rect(color = 'black', fill = NA),
          strip.text = element_text(size = 11, color = 'black'),
          axis.text = element_text(size = 9, color = 'black'),
          plot.margin = margin(l = 20, r = 2, t = 2, b = 2),
          axis.title = element_text(size = 11, color = 'black'),
          axis.title.y = element_blank()
        ) +
        labs(x = "VX34 genome coordinates (Mb)") +
        scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
        coord_cartesian(ylim = c(0,180)),
      
      nrow = 5) +
  draw_label("Counts of TTAGGC per kb", x=0.005, y=0.5, vjust= 1.5, angle=90, size = 11)

# tel_count_plt

# Save plot
ggsave("../../figures/supplementary/telomere_density.png", tel_count_plt, width = 7, height = 9, dpi = 600)

