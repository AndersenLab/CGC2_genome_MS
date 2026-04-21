library(readr)
library(dplyr)
library(tidyr)
library(ape)
library(ggtree)
library(ggplot2)

geno <- readr::read_tsv("../../processed_data/derivative_analysis/genotype_matrix.tsv") %>% dplyr::filter(QX1410!=2)

geno_simp <- setdiff(names(geno), c("Site", "QX1410"))

counts <- geno %>%
  drop_na() %>%
  dplyr::mutate(
    n_2_NA = rowSums(dplyr::across(dplyr::all_of(geno_simp), ~ is.na(.)), na.rm = TRUE),
    n_2 = rowSums(dplyr::across(dplyr::all_of(geno_simp), ~ . == 2), na.rm = TRUE)
  )

summary_counts <- counts %>%
  dplyr::count(n_2)

h2 <- ggplot(summary_counts, aes(x = as.character(n_2), y = n)) +
  geom_col() +
  geom_text(aes(label = n), vjust = 0.5,angle=90,hjust=-0.1) +
  labs(
    x = "Number of AF16 derivatives\nwith ALT allele",
    y = "Number of SNV sites",
    title = ""
  ) +
  theme_bw() +
  expand_limits(y = max(summary_counts$n) * 1.2)

#cowplot::plot_grid(h1,h2,ncol=2)

privates <- geno %>%
  dplyr::filter(
    rowSums(across(-QX1410, ~ . == 2), na.rm = TRUE) == 1
  ) %>%
  drop_na() %>%
  summarise(across(-Site, ~ sum(. == 2, na.rm = TRUE)))


privates_long <- pivot_longer(privates,
                        cols = everything(),
                        names_to = "Strain",
                        values_to = "Count")

# classic ggplot invocation
h3 <- ggplot(data = privates_long %>% dplyr::filter(Strain!="QX1410"), aes(x = Strain, y = Count)) +
  geom_col() +
  geom_text(aes(label = Count), vjust = 0.5,angle=90,hjust=-0.1) +
  labs(title = "",
       x = "",
       y = "Number of private SNV sites") +
  theme_bw() +
  expand_limits(y = max(privates_long$Count) * 1.1)


snvcounts <-  cowplot::plot_grid(h2 +
                                   theme(panel.grid = element_blank(),
                                         axis.text = element_text(color="black")) +
                                   scale_y_continuous(expand = c(0,0))
                                   ,
                                 h3+
                                   theme(panel.grid = element_blank(),
                                         axis.text = element_text(angle=45,hjust=1,color="black")) +
                                   scale_y_continuous(expand = c(0,0)),
                                 align = "h",
                                 axis="tb",
                                 labels = c("a","b"))

ggsave(snvcounts,filename = "../../figures/FigureS1_SNV_ct.png",width = 7.5,height = 5.5,dpi = 600,device = 'png')
write.table(geno %>% dplyr::select(-QX1410),"../../tables/TableS3_SNV_counts_derivatives.tsv",sep = "\t",quote = F,row.names = F)


