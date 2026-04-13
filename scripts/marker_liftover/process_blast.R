library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

indels <- readr::read_tsv("../../processed_data/marker_liftover/markers/AFHK_indels.clean.tsv") %>%
  dplyr::mutate(flen=nchar(forward),rlen=nchar(reverse))

indels_forward <- indels %>% dplyr::select(marker,flen) %>% dplyr::mutate(marker=paste0(marker,"_F")) %>% dplyr::rename(qlen=flen)
indels_reverse <- indels %>% dplyr::select(marker,rlen) %>% dplyr::mutate(marker=paste0(marker,"_R")) %>% dplyr::rename(qlen=rlen)
indel_lengths <- rbind(indels_forward,indels_reverse)

af_index <- readr::read_tsv("../../processed_data/marker_liftover/libraries/AF16/c_briggsae.PRJNA10731.WS280.genomic.fa.fai", col_names = c("name","length","offset","linebases","linewidth","qualoffset"))
cgc2_index <- readr::read_tsv("../../processed_data/marker_liftover/libraries/CGC2/c_briggsae.CGC2.hifi.ONT.HiC.Feb2026.genome.fa.fai", col_names = c("name","length","offset","linebases","linewidth","qualoffset"))

af_hits <- readr::read_tsv("../../processed_data/marker_liftover/blastn_results/AF16_marker_hits.tsv", col_names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","strand","eval","bitscore"))
cgc2_hits <- readr::read_tsv("../../processed_data/marker_liftover/blastn_results/CGC2_marker_hits.tsv", col_names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","strand","eval","bitscore"))


cgc2_best <- cgc2_hits %>%
  dplyr::left_join(indel_lengths,by=c("qseqid"="marker")) %>%
  dplyr::arrange(qseqid) %>%
  dplyr::group_by(qseqid) %>%
  dplyr::filter(qend >= qlen-1) %>%
  #dplyr::filter(eval==min(eval) | bitscore == max(bitscore)) %>%
  dplyr::filter(length==max(length) | length >=max(qend)) %>%
  dplyr::filter(mismatch==min(mismatch)) %>%
  dplyr::filter(gapopen==min(gapopen)) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ymin=3,ymax=3.5)

af_best <- af_hits %>%
  dplyr::left_join(indel_lengths,by=c("qseqid"="marker")) %>%
  dplyr::arrange(qseqid) %>%
  dplyr::group_by(qseqid) %>%
  dplyr::filter(qend >= qlen-1) %>%
  #dplyr::filter(eval==min(eval) | bitscore == max(bitscore)) %>%
  dplyr::filter(length==max(length) | length >=max(qend)) %>%
  dplyr::filter(mismatch==min(mismatch)) %>%
  dplyr::filter(gapopen==min(gapopen)) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ymin=0,ymax=0.5)

af_single <- af_best %>% dplyr::filter(gsize == 1) %>%dplyr::rename(name=sseqid)
cgc2_single <- cgc2_best %>% dplyr::filter(gsize == 1 ) %>% dplyr::rename(name=sseqid)

af_multi <- af_best %>% dplyr::filter(gsize > 1) %>% dplyr::rename(name=sseqid)
cgc2_multi <- cgc2_best %>% dplyr::filter(gsize > 1) %>% dplyr::rename(name=sseqid)

multi_all <- af_multi %>% dplyr::left_join(cgc2_multi,by="qseqid",relationship = "many-to-many")

primers_mapped <- af_single %>% dplyr::left_join(cgc2_single,by="qseqid") %>%
  tidyr::separate(qseqid, into=c("primerID","direction"),sep = "_",remove = F) %>%
  dplyr::mutate(is_sameChrom=ifelse(name.x==name.y,T,F)) %>%
  dplyr::group_by(primerID) %>%
  dplyr::mutate(prim_set_num=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(primer_set=ifelse(prim_set_num==1,"mult.mapped pair","paired"))

primers_complete <- primers_mapped %>%
  dplyr::filter(prim_set_num>1) %>%
  dplyr::group_by(primerID) %>%
  #dplyr::mutate(anyNAg=ifelse(any(is.na())))
  dplyr::mutate(newpos.x=min(sstart.x)+(max(send.x)-min(sstart.x))/2,newpos.y=min(sstart.y)+(max(send.y)-min(sstart.y))/2) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(lentest=abs(newpos.x-newpos.y))

af_counts <- af_multi %>%
  dplyr::filter(!grepl("cb25", name)) %>%
  dplyr::count(qseqid, name) %>%
  dplyr::group_by(qseqid) %>%
  dplyr::mutate(total = sum(n)) %>%
  dplyr::arrange(dplyr::desc(n), name, .by_group = TRUE) %>%
  dplyr::mutate(local_order = dplyr::row_number(), n_bars = dplyr::n()) %>%
  dplyr::ungroup()

af_qseqid_layout <- af_counts %>%
  dplyr::distinct(qseqid, total, n_bars) %>%
  dplyr::arrange(dplyr::desc(total)) %>%
  dplyr::mutate(start = cumsum(dplyr::lag(n_bars + 1, default = 0)), center = start + (n_bars - 1) / 2)

af_counts <- af_counts %>%
  dplyr::inner_join(dplyr::select(af_qseqid_layout, qseqid, start), by = "qseqid") %>%
  dplyr::mutate(x = start + local_order - 1)

cgc2_counts <- cgc2_multi %>%
  #dplyr::filter(!grepl("cb25", name)) %>%
  dplyr::count(qseqid, name) %>%
  dplyr::group_by(qseqid) %>%
  dplyr::mutate(total = sum(n)) %>%
  dplyr::arrange(dplyr::desc(n), name, .by_group = TRUE) %>%
  dplyr::mutate(local_order = dplyr::row_number(), n_bars = dplyr::n()) %>%
  dplyr::ungroup()

cgc2_qseqid_layout <- cgc2_counts %>%
  dplyr::distinct(qseqid, total, n_bars) %>%
  dplyr::arrange(dplyr::desc(total)) %>%
  dplyr::mutate(start = cumsum(dplyr::lag(n_bars + 1, default = 0)), center = start + (n_bars - 1) / 2)

cgc2_counts <- cgc2_counts %>%
  dplyr::inner_join(dplyr::select(cgc2_qseqid_layout, qseqid, start), by = "qseqid") %>%
  dplyr::mutate(x = start + local_order - 1)

primers_unmap <- primers_complete %>% 
  dplyr::filter(is.na(lentest)) %>%
  dplyr::filter(direction=="R") %>%
  dplyr::mutate(primer_set="CGC2 failed to map")

LIFT <- ggplot() + 
  geom_segment(data=primers_mapped %>% dplyr::filter(is_sameChrom==T & prim_set_num==1) %>% dplyr::rename(name=name.x),aes(x=sstart.x/1e6,xend=sstart.y/1e6,y=ymax.x,yend=ymin.y),color="black",linetype="11") +
  geom_segment(data=primers_complete %>% dplyr::filter(is_sameChrom==T) %>% dplyr::rename(name=name.x),aes(x=newpos.x/1e6,xend=newpos.y/1e6,y=ymax.x,yend=ymin.y),color="black") +
  geom_segment(data=primers_unmap %>% dplyr::rename(name=name.x),aes(x=sstart.x/1e6,xend=send.x/1e6,y=ymin.x,yend=ymax.x),color= "black")+
  geom_rect(data=af_index %>% dplyr::filter(!grepl("cb",name) & !grepl("MtDNA",name)),aes(xmin = 1/1e6,xmax=length/1e6,ymin=0,ymax=0.5,fill="AF16")) + 
  geom_rect(data=cgc2_index %>% dplyr::filter(!grepl("MtDNA",name)),aes(xmin = 1/1e6,xmax=length/1e6,ymin=3,ymax=3.5,fill="CGC2")) + 
  facet_wrap(~name,nrow=6,ncol=1,scales = 'free_x',strip.position = 'left') +
  scale_x_continuous(expand = c(0.005, 0.005)) +
  scale_fill_manual(values = c("AF16"="#E6B3B3","CGC2"="#7570B3"), name="strain") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x = element_line(),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        strip.text = element_text(size=12),
        axis.text.x=element_text(color="black")) +
  xlab("Genome coordinates (Mb)")

unique_unpaired <- primers_mapped %>% 
  dplyr::filter(is_sameChrom==T & prim_set_num==1) %>%
  dplyr::select(-primerID,-direction,-eval.x,-bitscore.x,-eval.y,-bitscore.y,-qend.y,-ymin.x,-ymax.x,-ymin.y,-ymax.y,-qstart.x,-qend.x,-qstart.y,-qend.y,-pident.x,-pident.y,-gsize.y,-is_sameChrom,-prim_set_num,-primer_set,-gsize.x,-qlen.y)
colnames(unique_unpaired) <- c("primer","AF16_chrom","AF16_match_length","AF16_n_mismatch","AF16_n_gap","AF16_start","AF16_end","AF16_strand","primer_length","CGC2_chrom","CGC2_match_length","CGC2_n_mismatch","CGC2_n_gap","CGC2_start","CGC2_end","CGC2_strand")

tbl2 <- primers_complete %>% 
  dplyr::filter(is_sameChrom==T & prim_set_num>1) %>%
  dplyr::select(-primerID,-direction,-eval.x,-bitscore.x,-eval.y,-bitscore.y,-qend.y,-ymin.x,-ymax.x,-ymin.y,-ymax.y,-qstart.x,-qend.x,-qstart.y,-qend.y,-pident.x,-pident.y,-gsize.y,-is_sameChrom,-prim_set_num,-primer_set,-gsize.x,-qlen.y,-newpos.x,-newpos.y,-lentest)
colnames(tbl2) <- c("primer","AF16_chrom","AF16_match_length","AF16_n_mismatch","AF16_n_gap","AF16_start","AF16_end","AF16_strand","primer_length","CGC2_chrom","CGC2_match_length","CGC2_n_mismatch","CGC2_n_gap","CGC2_start","CGC2_end","CGC2_strand")

tbl3 <- cgc2_counts %>%
  dplyr::select(qseqid,total) %>%
  dplyr::distinct(qseqid,.keep_all = T)

tbl4 <- af_counts %>%
  dplyr::select(qseqid,total) %>%
  dplyr::distinct(qseqid,.keep_all = T)

unique_paired <- tbl2 %>%
  dplyr::arrange(primer) %>%
  tidyr::separate(primer,into=c("indel_id","orientation"),sep = "_",remove=F) %>%
  dplyr::group_by(indel_id) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(AF16_mismatches=AF16_n_mismatch+AF16_n_gap,CGC2_mismatches=CGC2_n_mismatch+CGC2_n_gap) %>%
  dplyr::select(primer,indel_id,orientation,AF16_chrom,AF16_start,AF16_end,AF16_mismatches,CGC2_chrom,CGC2_start,CGC2_end,CGC2_mismatches)

perfect_matches <- unique_paired %>% dplyr::filter(CGC2_mismatches==AF16_mismatches & AF16_mismatches==0) %>%
  dplyr::group_by(indel_id) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(gsize==2)

single_perf_match <- unique_paired %>% dplyr::filter(CGC2_mismatches==AF16_mismatches & AF16_mismatches==0) %>%
  dplyr::group_by(indel_id) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(gsize==1)

single_perf_match_pair <- unique_paired %>% dplyr::filter(CGC2_mismatches==AF16_mismatches & AF16_mismatches==1)

near_perf_matches <- rbind(single_perf_match %>% dplyr::select(-gsize),single_perf_match_pair)

multi_mapping_summ <- bind_rows(unique_unpaired%>%
            tidyr::separate(primer,into=c("indel_id","orientation"),sep = "_",remove=F) %>%
            dplyr::mutate(AF16_mismatches=AF16_n_mismatch+AF16_n_gap,CGC2_mismatches=CGC2_n_mismatch+CGC2_n_gap) %>%
            dplyr::select(primer,indel_id,orientation,AF16_chrom,AF16_start,AF16_end,AF16_mismatches,CGC2_chrom,CGC2_start,CGC2_end,CGC2_mismatches),
          tbl4 %>%
            dplyr::rename(primer=qseqid,total_locations_AF16=total) %>%
            dplyr::left_join(tbl3 %>% dplyr::rename(primer=qseqid,total_locations_CGC2=total),by="primer")) %>%
  dplyr::mutate(total_locations_CGC2=ifelse(is.na(total_locations_CGC2),1,total_locations_CGC2),total_locations_AF16=ifelse(is.na(total_locations_AF16),1,total_locations_AF16)) %>%
  dplyr::arrange(primer)

write.table(rbind(perfect_matches,near_perf_matches),"../../tables/TableS8_markerLift_positions.tsv",quote = F,sep = "\t",row.names = F)
write.table(multi_mapping_summ,"../../tables/TableS9_markerLift_positions_multi_mapping_summ.tsv",quote = F,sep = "\t",row.names = F)
              
ggsave(LIFT + theme(strip.background = element_blank()),filename = "../figures/LIFT.png",width = 7,height = 4,dpi = 600,device = 'png')

multi_af_plot <- ggplot2::ggplot(af_counts, ggplot2::aes(x = x, y = n, fill = name)) +
  ggplot2::geom_col(width = 0.8) +
  ggplot2::geom_text(ggplot2::aes(label = n), hjust = -0.2, size = 3) +
  ggplot2::coord_flip() +
  ggplot2::scale_x_continuous(breaks = af_qseqid_layout$center, labels = af_qseqid_layout$qseqid) +
  ggplot2::labs(x = "qseqid", y = "Number of hits", fill = "Chromosome") +
  ggplot2::theme_bw()

multi_cgc2_plot <- ggplot2::ggplot(cgc2_counts, ggplot2::aes(x = x, y = n, fill = name)) +
  ggplot2::geom_col(width = 0.8) +
  ggplot2::geom_text(ggplot2::aes(label = n), hjust = -0.2, size = 3) +
  ggplot2::coord_flip() +
  ggplot2::scale_x_continuous(breaks = cgc2_qseqid_layout$center, labels = cgc2_qseqid_layout$qseqid) +
  ggplot2::labs(x = "qseqid", y = "Number of hits", fill = "Chromosome") +
  ggplot2::theme_bw() 

panel <- cowplot::plot_grid(
  multi_cgc2_plot + 
    ggplot2::theme(legend.position = "none", 
                   axis.title = element_blank(),
                   panel.grid = element_blank(),
                   axis.text=element_text(color="black")
                   ) + 
    expand_limits(y = max(cgc2_counts$n) * 1.1) +
    ggtitle("CGC2"),
  multi_af_plot + 
    ggplot2::theme(legend.position = "inside",
                   legend.position.inside = c(0.72,0.87),
                   axis.title = element_blank(),
                   panel.grid = element_blank(),
                   axis.text=element_text(color="black")
                   )+ 
    expand_limits(y = max(af_counts$n) * 1.1) +
    ggtitle("AF16"),
  labels = c("a","b"),
  rel_widths = c(1,1),
  nrow = 1,
  align = "h",
  axis = "tb")

MMAP <- cowplot::ggdraw(panel) +
  ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 20, l = 20)) +
  cowplot::draw_label("Primer ID", x = 0.005, y = 0.5, angle = 90, vjust = 0.5) +
  cowplot::draw_label("Number of BLASTn hits", x = 0.5, y = -0.005, vjust = 0.5)

ggsave(MMAP,filename = "../../figures/supplementary/FigureS13_MMAP.png",width = 7,height = 9,dpi = 600,device = 'png',bg="white")
