library(dplyr)
library(tidyr)
library(readr)
library(phylotools)

tsv <- readr::read_tsv("../../processed_data/marker_liftovers/markers/AFHK_indels.tsv")
colnames(tsv) <- c("marker","forward","reverse","ref","note")

tsv_clean <- tsv %>%
  dplyr::mutate(forward=toupper(forward),reverse=toupper(reverse))

write.table(tsv_clean,"../../processed_data/marker_liftovers/markers/AFHK_indels.clean.tsv",quote = F,sep="\t",row.names = F)
write.table(tsv_clean,"../../tables/TableS7_AFHK_indels_info.tsv",quote = F,sep="\t",row.names = F)

forward <- tsv_clean %>%
  dplyr::select(marker,forward) %>%
  dplyr::mutate(marker=paste0(marker,"_F")) %>%
  dplyr::rename(seq.text=forward,seq.name=marker)

reverse <- tsv_clean %>%
  dplyr::select(marker,reverse) %>%
  dplyr::mutate(marker=paste0(marker,"_R")) %>%
  dplyr::rename(seq.text=reverse,seq.name=marker)

primers <- rbind(forward,reverse)
dat2fasta(primers, outfile = "../../processed_data/marker_liftovers/markers/AFHK.primers.fasta")
