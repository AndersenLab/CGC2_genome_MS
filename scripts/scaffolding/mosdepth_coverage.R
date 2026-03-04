library(ggplot2)
library(dplyr)
library(readr)

hifi_cov <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_HiC/alignment/mosdepth_cov/hifi_chromV_1kb_cov.regions.mosdepth.bed", col_names = c("chrom", "start", "end", "cov")) %>%
  dplyr::mutate(mid_point = (start + end) / 2)

per_base_hifi <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_HiC/alignment/mosdepth_cov/hifi_per_base_cov.bed", col_names = c("chrom", "start", "end", "cov")) %>%
  dplyr::filter(end <= 335000) %>%
  dplyr::summarize(rDNA_cov = mean(cov))

mean_cov_hifi <- hifi_cov %>%
  dplyr::summarize(mean_cov = mean(cov),
                   mean_cov_first335 = )

ont_cov <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_HiC/alignment/mosdepth_cov/ONT_chromV_1kb_cov.regions.mosdepth.bed", col_names = c("chrom", "start", "end", "cov")) %>%
  dplyr::mutate(mid_point = (start + end) / 2) %>%
  dplyr::mutate(cov_20 = cov / 20)

per_base_ONT <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/CGC2_HiC/alignment/mosdepth_cov/ONT_per_base_cov.bed", col_names = c("chrom", "start", "end", "cov")) %>%
  dplyr::filter(end <= 335000) %>%
  dplyr::summarize(rDNA_cov = mean(cov))


mean_cov_ont <- ont_cov %>%
  dplyr::summarize(mean_cov = mean(cov))



rdna_cov <- ggplot() + 
  geom_line(data = ont_cov, aes(x = mid_point / 1000, y = cov_20, color = "ONT coverage / 20")) +
  geom_point(data = ont_cov, aes(x = mid_point / 1000, y = cov_20, color = "ONT coverage / 20")) + 
  geom_line(data = hifi_cov, aes(x = mid_point / 1000, y = cov, color = "HiFi coverage")) +
  geom_point(data = hifi_cov, aes(x = mid_point / 1000, y = cov, color = "HiFi coverage")) +
  scale_color_manual(values = c("HiFi coverage" = "deeppink", "ONT coverage / 20" = "steelblue")) +
  geom_hline(yintercept = mean_cov_hifi$mean_cov, linetype = "dashed") +
  geom_vline(xintercept = 335) +
  theme_bw() + 
  theme(
    axis.title = element_text(size = 14, color = 'black'),
    axis.text = element_text(size = 12, color = 'black'),
    legend.position = "inside",
    legend.position.inside = c(0.75,0.75),
    legend.text = element_text(size = 12, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    legend.box.background = element_rect(color = "black")
  ) +
  labs(x = "CGC2 chromosome V position (kb)", y = "Mean coverage / kb", color = "Seq technology") +
  coord_cartesian(xlim = c(0, 750)) +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.001)), breaks = seq(0, 1000, 100)) +
  scale_y_continuous(expand = expansion(mult = c(0.005, 0.01)))
rdna_cov

ggsave("/vast/eande106/projects/Lance/THESIS_WORK/manuscript_repos/CGC2_genome_MS/figures/supplementary/rDNA_coverage.png", rdna_cov, width = 7.5, height = 5, dpi = 600)
