library(dplyr)
library(readr)
library(ggplot2)
library(cowplot)

# Hi-C contig coverage:
cov <- readr::read_tsv("../../processed_data/scaffolding/HiC_contig_coverage.tsv") %>%
  dplyr::select(contig, coverage, meandepth) 

# Loading in scaffolds and corresponding contigs
chrom_cont <- readr::read_tsv("../../processed_data/scaffolding/scaffold_contig_IDs.tsv", col_names = c("scaffold","contig")) %>%
  dplyr::mutate(chrom = ifelse(scaffold == "scaffold_1","X", # assigning scaffolds to their chromosome ID
                               ifelse(scaffold == "scaffold_2","V",
                                      ifelse(scaffold == "scaffold_3","IV",
                                             ifelse(scaffold == "scaffold_4","II",
                                                    ifelse(scaffold == "scaffold_5","I",
                                                           ifelse(scaffold == "scaffold_6","III","scaffold_7"))))))) %>%
  dplyr::left_join(cov, by = "contig") %>% # adding Hi-C coverage per contig information to data frame
  dplyr::mutate(contig = factor(contig, levels = c("ptg000009l","ptg000001l","ptg000004l","ptg000007l",
                                                   "ptg000005l", "ptg000003l","ptg000018l", "ptg000023l", "ptg000016l", "ptg000031l",
                                                   "ptg000024l", "ptg000002l", "ptg000011l"))) %>%
  mutate(chrom = factor(chrom, levels = c("I","II","III","IV","V","X","scaffold_7"))) 

# Plotting Hi-C coverage per 
hic_coverage <- ggplot(chrom_cont, aes(x = contig, y = meandepth, fill = chrom)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(
    breaks = seq(0, 800, 50),
    minor_breaks = seq(0, 800, 25),
    expand = c(0, 5)) +
  coord_cartesian(ylim = c(0,750)) +
  scale_fill_manual(values = c("I" = "seagreen", "II" = "blue", "III" = "red", "IV" = "green", "V" = "orange","X" = "purple", "scaffold_7" = "gray")) +
  theme(
    axis.title = element_text(size = 10, color = 'black'),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10, color = 'black', angle = 60, hjust = 1.05),
    panel.background = element_blank(),
    panel.grid.major= element_line(color = 'gray80'),
    legend.text = element_text(size = 10, color = 'black'),
    legend.title = element_text(size = 10, color = 'black'),
    legend.position = 'top',
    # legend.justification = "left", 
    legend.justification = c(0, 0.5),
    legend.margin = margin(0, 0, 0, 0),     
    legend.box.margin = margin(0, 0, 0, 0), 
    legend.box.just = "left",
    plot.margin = margin(b = 10),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 10, color = 'black'),
    panel.border = element_rect(fill = NA)) +
  labs(y = "Mean sequencing depth", fill = "Chromosome \nassignment") 
hic_coverage


# Scaffold 7 aligned to scaffolded CGC2 genome
scaff_7 <- readr::read_tsv("../../processed_data/genome_genome_alignments/scaffold7_CGC2.transformed.tsv", col_names = c("s7S","s7E","CGS","CGE","L1","L2","IDY","LENR","LENQ","scaffold_7","CGC2_scaffold")) %>%
  dplyr::mutate(CGC2_scaffold = ifelse(CGC2_scaffold == "scaffold_1","scaffold_1 (Chromosome X)", 
                                       ifelse(CGC2_scaffold == "scaffold_2","scaffold_2 (Chromosome V)",
                                              ifelse(CGC2_scaffold == "scaffold_3","scaffold_3 (Chromosome IV)",
                                                     ifelse(CGC2_scaffold == "scaffold_4","scaffold_4 (Chromosome II)",
                                                            ifelse(CGC2_scaffold == "scaffold_5","scaffold_5 (Chromosome I)",
                                                                   ifelse(CGC2_scaffold == "scaffold_6","scaffold_6 (Chromosome III)", CGC2_scaffold))))))) 

# Plotting the alingment(s) of scaffold 7 to CGC2 genome
sc_7_aln <- ggplot(scaff_7 %>% dplyr::filter(L1 > 30000)) + # removing short alignments of scaffold 7
  geom_segment(aes(x = CGS / 1e6, xend = CGE / 1e6, y = s7S / 1e6, yend = s7E / 1e6), color = "black", linewidth = 1) +
  facet_wrap(~CGC2_scaffold, scales = "free") +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size = 10, color = 'black'),
    panel.grid = element_blank(),
    strip.text = element_text(size = 10, color = 'black'),
    axis.text = element_text(size = 10, color = 'black'),
    panel.border = element_rect(fill = NA)) +
  labs(y = "Scaffold 7 coordinates (Mb)", x = "CGC2 genome coordinates (Mb)")
sc_7_aln


# Create final plot:
final_plot <- cowplot::plot_grid(
  hic_coverage, sc_7_aln,
  nrow = 1,
  align = "h",
  axis = "tb"
)
final_plot


# Save the plot
ggsave("../../figures/supplementary/duplicate_haplotig.png", final_plot, width = 7, height = 5, dpi = 600)

