library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

chrom_order <- c("I", "II", "III", "IV", "V", "X")
rect_half_height <- 0.08
offset_spacer = 300000

af_gff<- readr::read_tsv("../../processed_data/orthofinder/gff/c_briggsae.PRJNA10731.WS280.protein_coding.longest.gff",col_names = c("seqid","source","type","start","end","score","strand","phase","attributes"))
qx_gff<- readr::read_tsv("../../processed_data/orthofinder/gff/c_briggsae.QX1410_20250929.csq.longest.gff",col_names = c("seqid","source","type","start","end","score","strand","phase","attributes"))
cgc2_gff<- readr::read_tsv("../../processed_data/orthofinder/gff/CGC2.softMasked.braker.longest.gff",col_names = c("seqid","source","type","start","end","score","strand","phase","attributes"))

qx_fai <- readr::read_tsv("../../processed_data/genomes/c_briggsae.QX1410.nanopore.Feb2020.genome.fa.fai", col_names = c("chrom", "length", "offset", "line_bases", "line_width"))
af_fai <- readr::read_tsv("../../processed_data/genomes/c_briggsae.AF16.PRJNA10731.WS276.genome.fa.fai", col_names = c("chrom", "length", "offset", "line_bases", "line_width")) %>%
  dplyr::filter(chrom %in% chrom_order) %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order)) %>%
  dplyr::arrange(chrom) %>%
  dplyr::mutate(chrom = as.character(chrom)) %>%
  dplyr::filter(!is.na(chrom))

cgc2_fai <- readr::read_tsv("../../processed_data/genomes/c_briggsae.CGC2.hifi.ONT.HiC.Feb2026.genome.fa.fai", col_names = c("chrom", "length", "offset", "line_bases", "line_width")) %>% 
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order)) %>%
  dplyr::arrange(chrom) %>%
  dplyr::mutate(chrom = as.character(chrom))

# orthos <- readr::read_tsv("../processed_data/orthofinder/OrthoFinder_nN2/Results_CGC2_nN2/Orthogroups/Orthogroups.tsv") %>% 
#   dplyr::rename(QX1410=c_briggsae.QX1410_20250929.csq.longest.protein,
#                 CGC2=CGC2.softMasked.braker.longest.protein,
#                 AF16=c_briggsae.PRJNA10731.WS280.protein_coding.longest.prot)

orthos <- readr::read_tsv("../../processed_data/orthofinder/orthofinder_results/Orthogroups/Orthogroups.tsv") %>% 
  dplyr::rename(QX1410=c_briggsae.QX1410_20250929.csq.longest.protein,
                CGC2=CGC2.softMasked.braker.longest.protein,
                AF16=c_briggsae.PRJNA10731.WS280.protein_coding.longest.prot,
                N2=N2.WBonly.WS283.PConly.prot)

orthos_counts <- orthos %>%
  dplyr::mutate(
    dplyr::across(
      -Orthogroup,
      ~ ifelse(
        is.na(.),
        NA_integer_,
        stringr::str_count(., ",") + 1),.names = "{.col}_counts")) %>%
  dplyr::select(-N2_counts) %>%
  dplyr::filter(dplyr::if_all(dplyr::ends_with("_counts"),~ . == 1 | is.na(.))) %>%
  dplyr::filter(!dplyr::if_all(dplyr::ends_with("_counts"),~ is.na(.)))   %>%
  dplyr::rowwise() %>%
  dplyr::filter(sum(CGC2_counts,AF16_counts,QX1410_counts,na.rm = T)>=2) %>%
  dplyr::ungroup()
      
single_complete <- orthos_counts %>%
  dplyr::filter(dplyr::if_all(dplyr::ends_with("_counts"),~ . == 1)) %>%
  dplyr::mutate(class="COMPLETE")

single_absent <- orthos_counts %>%
  dplyr::filter(dplyr::if_any(dplyr::ends_with("_counts"),~ is.na(.))) %>%
  dplyr::mutate(class=ifelse(is.na(AF16),"AF16_ABSENCE",ifelse(is.na(QX1410),"QX1410_ABSENCE",ifelse(is.na(CGC2),"CGC2_ABSENCE","PRIVATE"))))

single_class <- rbind(single_complete,single_absent)

cgc2_transcripts <- cgc2_gff %>%
  dplyr::filter(type=="mRNA") %>%
  dplyr::select(seqid,start,end,strand,attributes) %>%
  dplyr::mutate(attributes=gsub(";Parent=.*","",attributes)) %>%
  dplyr::mutate(attributes=gsub("ID=","",attributes)) %>%
  dplyr::rename(tranname=attributes,cgc2_chrom=seqid,cgc2_start=start,cgc2_end=end,cgc2_strand=strand) 

af_transcripts <- af_gff %>%
  dplyr::filter(type=="mRNA") %>%
  dplyr::select(seqid,start,end,strand,attributes) %>%
  dplyr::mutate(attributes=gsub(";Parent=.*","",attributes)) %>%
  dplyr::mutate(attributes=gsub("ID=Transcript:","Transcript_",attributes)) %>%
  dplyr::rename(tranname=attributes,af_chrom=seqid,af_start=start,af_end=end,af_strand=strand) 

qx_transcripts <- qx_gff %>%
  dplyr::filter(type=="mRNA") %>%
  dplyr::select(seqid,start,end,strand,attributes) %>%
  dplyr::mutate(attributes=gsub(";Parent=.*","",attributes)) %>%
  dplyr::mutate(attributes=gsub("ID=transcript:","transcript_",attributes)) %>%
  dplyr::rename(tranname=attributes,qx_chrom=seqid,qx_start=start,qx_end=end,qx_strand=strand) 

single_anno <- single_class %>%
  dplyr::left_join(cgc2_transcripts,by=c("CGC2"="tranname")) %>%
  dplyr::left_join(af_transcripts,by=c("AF16"="tranname")) %>%
  dplyr::left_join(qx_transcripts,by=c("QX1410"="tranname"))

cgc2_offsets <- cgc2_fai %>%
  dplyr::mutate(
    chrom = factor(chrom, levels = chrom_order)
  ) %>%
  dplyr::arrange(chrom) %>%
  dplyr::mutate(
    chrom_index = dplyr::row_number(),
    cgc2_offset = dplyr::lag(cumsum(length), default = 0) + (chrom_index - 1) * offset_spacer
  ) %>%
  dplyr::transmute(
    cgc2_chrom = as.character(chrom),
    cgc2_offset
  )

af_offsets <- af_fai %>%
  dplyr::mutate(
    chrom = factor(chrom, levels = chrom_order)
  ) %>%
  dplyr::arrange(chrom) %>%
  dplyr::mutate(
    chrom_index = dplyr::row_number(),
    af_offset = dplyr::lag(cumsum(length), default = 0) + (chrom_index - 1) * offset_spacer
  ) %>%
  dplyr::transmute(
    af_chrom = as.character(chrom),
    af_offset
  )

qx_offsets <- qx_fai %>%
  dplyr::mutate(
    chrom = factor(chrom, levels = chrom_order)
  ) %>%
  dplyr::arrange(chrom) %>%
  dplyr::mutate(
    chrom_index = dplyr::row_number(),
    qx_offset = dplyr::lag(cumsum(length), default = 0) + (chrom_index - 1) * offset_spacer
  ) %>%
  dplyr::transmute(
    qx_chrom = as.character(chrom),
    qx_offset
  )

single_anno_renamed <- single_anno
colnames(single_anno_renamed) <- colnames(single_anno_renamed) |>
  gsub("^cgc2_", "CGC2_", x = _) |>
  gsub("^qx_", "QX1410_", x = _) |>
  gsub("^af_", "AF16_", x = _)

write.table(single_anno_renamed,"../../tables/TableS11_single_copy_orthologs_3spp.tsv",sep = "\t",quote = F,row.names = F)

single_anno_shifted <- single_anno %>%
  dplyr::left_join(cgc2_offsets, by = "cgc2_chrom") %>%
  dplyr::left_join(af_offsets, by = "af_chrom") %>%
  dplyr::left_join(qx_offsets, by = "qx_chrom") %>%
  dplyr::mutate(cgc2_start_shift = cgc2_start + cgc2_offset,
                cgc2_end_shift   = cgc2_end   + cgc2_offset,
                af_start_shift   = af_start   + af_offset,
                af_end_shift     = af_end     + af_offset,
                qx_start_shift   = qx_start   + qx_offset,
                qx_end_shift     = qx_end     + qx_offset)
    

# build plotting data from shifted coordinates ---------------------------------

plot_df <- single_anno_shifted %>%
  dplyr::filter(
    is.na(qx_chrom)   | qx_chrom   %in% chrom_order,
    is.na(cgc2_chrom) | cgc2_chrom %in% chrom_order,
    is.na(af_chrom)   | af_chrom   %in% chrom_order
  ) %>%
  dplyr::mutate(
    qx_mid_shift   = qx_start_shift   + (qx_end_shift   - qx_start_shift) / 2,
    cgc2_mid_shift = cgc2_start_shift + (cgc2_end_shift - cgc2_start_shift) / 2,
    af_mid_shift   = af_start_shift   + (af_end_shift   - af_start_shift) / 2
  )

# genome rectangles from fai + cumulative offsets ------------------------------

qx_rects <- qx_fai %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order)) %>%
  dplyr::arrange(chrom) %>%
  dplyr::filter(!is.na(chrom)) %>%
  dplyr::left_join(qx_offsets, by = c("chrom" = "qx_chrom")) %>%
  dplyr::transmute(
    chrom = as.character(chrom),
    genome = "QX1410",
    y = 2,
    xmin = qx_offset,
    xmax = qx_offset + length,
    ymin = y - rect_half_height,
    ymax = y + rect_half_height)

cgc2_rects <- cgc2_fai %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order)) %>%
  dplyr::arrange(chrom) %>%
  dplyr::filter(!is.na(chrom)) %>%
  dplyr::left_join(cgc2_offsets, by = c("chrom" = "cgc2_chrom")) %>%
  dplyr::transmute(
    chrom = as.character(chrom),
    genome = "CGC2",
    y = 3,
    xmin = cgc2_offset,
    xmax = cgc2_offset + length,
    ymin = y - rect_half_height,
    ymax = y + rect_half_height
  )

cgc2_rects2 <- cgc2_fai %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order)) %>%
  dplyr::arrange(chrom) %>%
  dplyr::filter(!is.na(chrom)) %>%
  dplyr::left_join(cgc2_offsets, by = c("chrom" = "cgc2_chrom")) %>%
  dplyr::transmute(
    chrom = as.character(chrom),
    genome = "CGC2",
    y = 0,
    xmin = cgc2_offset,
    xmax = cgc2_offset + length,
    ymin = y - rect_half_height,
    ymax = y + rect_half_height
  )

af_rects <- af_fai %>%
  dplyr::mutate(chrom = factor(chrom, levels = chrom_order)) %>%
  dplyr::arrange(chrom) %>%
  dplyr::filter(!is.na(chrom)) %>%
  dplyr::left_join(af_offsets, by = c("chrom" = "af_chrom")) %>%
  dplyr::transmute(
    chrom = as.character(chrom),
    genome = "AF16",
    y = 1,
    xmin = af_offset,
    xmax = af_offset + length,
    ymin = y - rect_half_height,
    ymax = y + rect_half_height
  )

rect_df <- dplyr::bind_rows(cgc2_rects,cgc2_rects2,qx_rects,af_rects)
# optional chromosome labels centered on each block ----------------------------

label_df <- rect_df %>%
  dplyr::group_by(genome, chrom, y, ymin, ymax) %>%
  dplyr::summarise(
    x = (dplyr::first(xmin) + dplyr::first(xmax)) / 2,
    .groups = "drop"
  )

# connection segments: keep valid pairs only -----------------------------------

seg_qx_cgc2 <- plot_df %>%
  dplyr::filter(!is.na(qx_mid_shift), !is.na(cgc2_mid_shift)) %>%
  dplyr::transmute(
    class,
    type = dplyr::if_else(qx_chrom == cgc2_chrom, "intra", "inter"),
    x = cgc2_mid_shift,
    xend = qx_mid_shift,
    y = 3 - rect_half_height,
    yend = 2 + rect_half_height
  )

seg_cgc2_af <- plot_df %>%
  dplyr::filter(!is.na(af_mid_shift), !is.na(cgc2_mid_shift)) %>%
  dplyr::transmute(
    class,
    type = dplyr::if_else(cgc2_chrom == af_chrom, "intra", "inter"),
    x = af_mid_shift,
    xend = cgc2_mid_shift,
    y = 1 - rect_half_height,
    yend = 0 + rect_half_height
  )

seg_qx_af <- plot_df %>%
  dplyr::filter(!is.na(qx_mid_shift), !is.na(af_mid_shift)) %>%
  dplyr::transmute(
    class,
    type = dplyr::if_else(qx_chrom == af_chrom, "intra", "inter"),
    x = qx_mid_shift,
    xend = af_mid_shift,
    y = 2 - rect_half_height,
    yend = 1 + rect_half_height
  )

seg_df <- dplyr::bind_rows(seg_qx_cgc2, seg_cgc2_af,seg_qx_af)

# plot ------------------------------------------------------------------------
leg_scaff <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = seg_df,
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = class),
    alpha = 1,
    linewidth = 1
  ) +
  ggplot2::geom_rect(
    data = rect_df,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = genome),
    color = "black",
    linewidth = 0.3
  ) +
  ggplot2::geom_text(
    data = label_df %>% dplyr::filter(genome=="QX1410"),
    ggplot2::aes(x = x, y = ymax+0.06, label = chrom),
    size = 4
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "CGC2" = "#7570B3",
      "AF16" = "#E6B3B3",
      "QX1410" = "#53886C"
    ),
    breaks = c("QX1410", "CGC2", "AF16"),
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "COMPLETE" = "grey50",
      "AF16_ABSENCE" = "#E6B3B3",
      "QX1410_ABSENCE" = "#53886C",
      "CGC2_ABSENCE"= "#7570B3"
    ),
    breaks = c("COMPLETE", "AF16_ABSENCE", "QX1410_ABSENCE","CGC2_ABSENCE"),
    labels = c(
      "COMPLETE" = "Present (all)",
      "AF16_ABSENCE" = "Absent in AF16",
      "QX1410_ABSENCE" = "Absent in QX1410",
      "CGC2_ABSENCE"= "Absent in CGC2"
    )
  )+
  ggplot2::scale_y_continuous(
    breaks = c(0 ,1, 2, 3),
    labels = c("CGC2","AF16","QX1410","CGC2"),
    limits = c(-0.1, 3.3),
    expand = c(0.01, 0)
  ) +
  ggplot2::scale_x_continuous(
    labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
    name = "Genome coordinates (Mb)",
    expand = c(0.01, 0)
  ) +
  ggplot2::labs(
    x = "Genome coordinates (Mb)",
    y = NULL,
    color = "SC ortholog\nclass:",
    fill = "Genome:"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "white"),
    panel.grid = ggplot2::element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.ticks.y=element_blank(),
    axis.line.y=element_blank(),
    axis.text.y=element_blank()
  ) +
  ggplot2::guides(
    fill = guide_legend(order = 1), 
    color = guide_legend(order = 2,
                         override.aes = list(linewidth = 3))  
  )
leg <- cowplot::get_legend(leg_scaff)


p <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = seg_df %>% dplyr::filter(type=="intra" & class=="COMPLETE"),
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = class),
    alpha = 0.15,
    linewidth = 0.1
  ) +
  ggplot2::geom_rect(
    data = rect_df,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = genome),
    color = "black",
    linewidth = 0.3
  ) +
  ggplot2::geom_text(
    data = label_df %>% dplyr::filter(genome=="CGC2" & y==3),
    ggplot2::aes(x = x, y = ymax+0.15, label = chrom),
    size = 4
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "CGC2" = "#7570B3",
      "AF16" = "#E6B3B3",
      "QX1410" = "#53886C"
    )
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "COMPLETE" = "grey50",
      "AF16_ABSENCE" = "#E6B3B3",
      "QX1410_ABSENCE" = "#53886C",
      "CGC2_ABSENCE"= "#7570B3"
    )
  ) +
  ggplot2::scale_y_continuous(
    breaks = c(0 ,1, 2, 3),
    labels = c("CGC2","AF16","QX1410","CGC2"),
    limits = c(-0.1, 3.3),
    expand = c(0.01, 0)
  ) +
  ggplot2::scale_x_continuous(
    labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
    name = "Genome coordinates (Mb)",
    expand = c(0.01, 0)
  ) +
  ggplot2::labs(
    x = "Genome coordinates (Mb)",
    y = NULL,
    color = "",
    fill = "Genome"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "white"),
    panel.grid = ggplot2::element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.ticks.y=element_blank(),
    axis.line.y=element_blank(),
    axis.text.y=element_blank(),
    axis.text.x=element_text(color="black")
  )

p2 <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = seg_df %>% dplyr::filter(type=="intra" & class!="COMPLETE"),
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = class),
    alpha = 1,
    linewidth = 0.15
  ) +
  ggplot2::geom_rect(
    data = rect_df,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = genome),
    color = "black",
    linewidth = 0.3
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "CGC2" = "#7570B3",
      "AF16" = "#E6B3B3",
      "QX1410" = "#53886C"
    )
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "COMPLETE" = "grey50",
      "AF16_ABSENCE" = "#E6B3B3",
      "QX1410_ABSENCE" = "#53886C",
      "CGC2_ABSENCE"= "#7570B3"
    )
  ) +
  ggplot2::scale_y_continuous(
    breaks = c(0 ,1, 2, 3),
    labels = c("CGC2","AF16","QX1410","CGC2"),
    limits = c(-0.1, 3.3),
    expand = c(0.01, 0)
  ) +
  ggplot2::scale_x_continuous(
    labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
    name = "Genome coordinates (Mb)",
    expand = c(0.01, 0)
  ) +
  ggplot2::labs(
    x = "Genome coordinates (Mb)",
    y = NULL,
    color = "",
    fill = "Genome"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "white"),
    panel.grid = ggplot2::element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.ticks.y=element_blank(),
    axis.line.y=element_blank(),
    axis.text.y=element_blank(),
    axis.text.x=element_text(color="black")
  )

class_counts <- single_anno_shifted %>%
  count(class)

pbar <- ggplot(class_counts, aes(x = reorder(class, -n), y = n,fill=class)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.3,angle=60,hjust=0) +
  theme_classic() +
  labs(
    title = "",
    x = "",
    y = "Number of SC orthologs") +
  ggplot2::scale_fill_manual(values = c("COMPLETE" = "grey50",
                                        "AF16_ABSENCE" = "#E6B3B3",
                                        "QX1410_ABSENCE" = "#53886C",
                                        "CGC2_ABSENCE"= "#7570B3"  ),
                             breaks=c("COMPLETE","AF16_ABSENCE","QX1410_ABSENCE","CGC2_ABSENCE"))+
  scale_y_continuous(expand= expansion(mult = c(0, 0.29), add = c(1, 0))) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'none',
        axis.text=element_text(color="black"))
 


pmain <- cowplot::plot_grid(p + theme(axis.title = element_blank(), legend.position = 'none'),
                            p2 + theme(legend.position = 'none'),
                            nrow=2,labels = c("a","b"),align = "v",axis = "lr",rel_heights = c(1,1.1))
                   

leg_bar <- cowplot::plot_grid(leg,pbar,NULL,nrow=3,rel_heights = c(1,1.05,0.05),labels=c("","c",""))

pleg <- cowplot::plot_grid(pmain,leg_bar,rel_widths = c(1,0.3),align = "h",axis = "bt")

ggsave(pleg,filename = "../../figures/Figure4_ORTHOS.png",width = 7,height = 5.5,dpi = 600,device = 'png',bg = "white")

alt_seg_qx_cgc2 <- plot_df %>%
  dplyr::filter(!is.na(qx_mid_shift), !is.na(cgc2_mid_shift)) %>%
  dplyr::transmute(
    class,
    type = dplyr::if_else(cgc2_chrom == af_chrom & qx_chrom == cgc2_chrom, "intra", "inter"),
    typelab=dplyr::if_else(type=="inter" & cgc2_chrom != af_chrom & cgc2_chrom != qx_chrom,"inter_both",ifelse(type=="inter" & cgc2_chrom != af_chrom,"inter_AF16","inter_QX1410")),
    x = cgc2_mid_shift,
    xend = qx_mid_shift,
    y = 3 - rect_half_height,
    yend = 2 + rect_half_height
  )

alt_seg_cgc2_af <- plot_df %>%
  dplyr::filter(!is.na(cgc2_mid_shift), !is.na(af_mid_shift)) %>%
  dplyr::transmute(
    class,
    type = dplyr::if_else(cgc2_chrom == af_chrom & qx_chrom == cgc2_chrom, "intra", "inter"),
    typelab=dplyr::if_else(type=="inter" & cgc2_chrom != af_chrom & cgc2_chrom != qx_chrom,"inter_both",ifelse(type=="inter" & cgc2_chrom != af_chrom,"inter_AF16","inter_QX1410")),
    x = qx_mid_shift,
    xend = af_mid_shift,
    y = 2 - rect_half_height,
    yend = 1 + rect_half_height
  )

alt_seg_df <- dplyr::bind_rows(alt_seg_qx_cgc2, alt_seg_cgc2_af) 

pinter <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = alt_seg_df %>% dplyr::filter(type=="inter" & class=="COMPLETE"),
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = typelab),
    alpha = 1,
    linewidth = 0.5
  ) +
  ggplot2::geom_rect(
    data = rect_df %>% dplyr::filter(y!=0),
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = genome),
    color = "black",
    linewidth = 0.3
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "CGC2" = "#7570B3",
      "AF16" = "#E6B3B3",
      "QX1410" = "#53886C"
    )
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "inter_both" = "#7570B3",
      "inter_AF16" = "#E6B3B3",
      "inter_QX1410" = "#53886C"
    ),
    breaks = c("inter_both","inter_AF16", "inter_QX1410"),
    labels = c(
      "inter_both" = "Translocated in CGC2",
      "inter_AF16" = "Translocated in AF16",
      "inter_QX1410" = "Translocated in QX1410"
    )
  )+
  ggplot2::scale_y_continuous(
    breaks = c(1, 2, 3),
    labels = c("AF16", "QX1410", "CGC2"),
    limits = c(0.90, 3.3),
    expand = c(0.01, 0)
  ) +
  ggplot2::scale_x_continuous(
    labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
    name = "Genome coordinates (Mb)",
    expand = c(0.01, 0)
  ) +
  ggplot2::labs(
    x = "Genome coordinates (Mb)",
    y = NULL,
    color = "Translocation\nclassification:",
    fill = "Genome:"
  ) +
  ggplot2::guides(
    fill = guide_legend(order = 1), 
    color = guide_legend(order = 2,
                         override.aes = list(linewidth = 2))  
  )+
  ggplot2::theme_classic() +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "white"),
    panel.grid = ggplot2::element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(),
    axis.ticks.y=element_blank(),
    axis.line.y=element_blank(),
    axis.text.y=element_blank(),
    axis.text.x=element_text(color="black")
  ) #warning indicates that we omit double rectangles for CGC2 - intentional, for simplicity

ggsave(pinter,filename = "../../figures/supplementary/FigureS15_TRANSLOC.png",width = 7,height = 3,dpi = 600,device = 'png',bg = "white")

N2_gff <- readr::read_tsv("../../processed_data/orthofinder/gff/N2.WBonly.WS283.PConly.longest.gff",col_names = c("seqid","source","type","start","end","score","strand","phase","attributes"))

aacid <- c(
  QX1410 = "../../processed_data/orthofinder/protein/c_briggsae.QX1410_20250929.csq.longest.protein.fa",
  CGC2   = "../../processed_data/orthofinder/protein/CGC2.softMasked.braker.longest.protein.fa",
  AF16   = "../../processed_data/orthofinder/protein/c_briggsae.PRJNA10731.WS280.protein_coding.longest.prot.fa",
  N2     = "../../processed_data/orthofinder/protein/N2.WBonly.WS283.PConly.prot.fa"
)

seq_table <- purrr::imap_dfr(aacid, \(file, sample) {
  fa <- Biostrings::readAAStringSet(file)
  tibble::tibble(
    seqid = names(fa),
    length = BiocGenerics::width(fa),
    sample = sample
  )
}) %>%
  dplyr::mutate(seqid=gsub(":","_",seqid))

N2_genes <- N2_gff %>%
  dplyr::filter(type=="gene") %>%
  dplyr::select(attributes) %>%
  dplyr::mutate(attributes=gsub(";Name=.*","",attributes)) %>%
  tidyr::separate(attributes,into=c("gene","pre_alias"),sep = ";") %>%
  tidyr::separate(pre_alias,into=c("alias","alias2"),sep=",") %>%
  dplyr::mutate(alias=gsub("Alias=","",alias)) %>%
  dplyr::mutate(gene=gsub("ID=","",gene))

N2_transcripts <- N2_gff %>%
  dplyr::filter(type=="mRNA") %>%
  dplyr::select(seqid,start,end,strand,attributes) %>%
  dplyr::mutate(attributes=gsub(";Name=.*","",attributes)) %>%
  tidyr::separate(attributes,into=c("tranname","parent"),sep = ";") %>%
  dplyr::mutate(tranname=gsub("ID=Transcript:","Transcript_",tranname)) %>%
  dplyr::mutate(parent=gsub("Parent=","",parent)) %>%
  dplyr::rename(N2_chrom=seqid,N2_start=start,N2_end=end,N2_strand=strand)%>%
  dplyr::left_join(N2_genes,by=c("parent"="gene")) 

orthos_wlen <- orthos_counts %>%
  #dplyr::filter(!is.na(AF16) & !is.na(CGC2) & !is.na(QX1410) & !is.na(N2)) %>%
  dplyr::filter(!grepl(",",N2)) %>%
  dplyr::left_join(N2_transcripts,by=c("N2"="tranname")) %>%
  dplyr::left_join(seq_table %>% dplyr::filter(sample=="QX1410"),by=c("QX1410"="seqid")) %>%
  dplyr::rename(QX1410_length=length) %>%
  dplyr::select(-sample) %>%
  dplyr::left_join(seq_table %>% dplyr::filter(sample=="AF16"),by=c("AF16"="seqid")) %>%
  dplyr::rename(AF16_length=length) %>%
  dplyr::select(-sample) %>%
  dplyr::left_join(seq_table %>% dplyr::filter(sample=="CGC2"),by=c("CGC2"="seqid")) %>%
  dplyr::rename(CGC2_length=length) %>%
  dplyr::select(-sample) %>%
  dplyr::left_join(seq_table %>% dplyr::filter(sample=="N2"),by=c("N2"="seqid")) %>%
  dplyr::rename(N2_length=length) %>%
  dplyr::select(-sample) %>%
  dplyr::mutate(QX1410_acc=QX1410_length/N2_length,
                AF16_acc=AF16_length/N2_length,
                CGC2_acc=CGC2_length/N2_length) %>%
  dplyr::mutate(QX1410_deviation=ifelse(QX1410_acc>=1,QX1410_acc-1,1-QX1410_acc),
                AF16_deviation=ifelse(AF16_acc>=1,AF16_acc-1,1-AF16_acc),
                CGC2_deviation=ifelse(CGC2_acc>=1,CGC2_acc-1,1-CGC2_acc)) %>%
  # dplyr::mutate(class=ifelse(AF16_deviation==CGC2_deviation & AF16_deviation==1, "identical",
  #                            ifelse(AF16_deviation < 0.05 & CGC2_deviation < 0.05, "near identical",
  #                                   ifelse(AF16_deviation > 0.05 & CGC2_deviation > 0.05 & AF16_acc > 1 & CGC2_acc > 1, "positive deviation",
  #                                          ifelse(AF16_deviation > 0.05 & CGC2_deviation > 0.05 & AF16_acc < 1 & CGC2_acc < 1, "negative deviation","discordant deviation"))))) %>%
  dplyr::mutate(class=ifelse(CGC2_deviation==0  | AF16_deviation==0,"Identical",
                             ifelse(AF16_deviation < 0.05 | CGC2_deviation < 0.05,"Near-identical",
                                    ifelse(abs(AF16_acc-CGC2_acc) < 0.05,"Deviated (Concordant)","Deviated (Discordant)"))))


orthos_3way <- orthos_wlen %>% dplyr::filter(!is.na(AF16) & !is.na(CGC2) & !is.na(N2))
write.table(orthos_3way %>%
              dplyr::rename(N2_alias=alias,N2_alias2=alias2) %>%
              dplyr::select(-QX1410,-parent,-QX1410_counts,-QX1410_length,-QX1410_acc,-QX1410_deviation) %>%
              dplyr::arrange(N2_chrom,N2_start),"../../tables/supplementary/TableS12_sc_ortholog_protein_lengths_wN2.tsv",sep = "\t",quote = F,row.names = F)

plac1 <- ggplot(orthos_3way) + 
  geom_point(aes(x=CGC2_acc-1,y=AF16_acc-1,color=class),size=0.7) + 
  coord_cartesian(xlim=c(-1,1.05),ylim=c(-1,1.5)) + 
  geom_rect(
    xmin = -0.1, xmax = 0.1,
    ymin = -0.1, ymax = 0.1,
    fill = NA,
    color = "black",
    linewidth = 0.8,
    linetype="11"
  ) +
  scale_color_manual(
    values = c(
      "Identical" = "blue",              # cool teal
      "Near-identical" = "#66C2A5",         # lighter cool teal
      "Deviated (Concordant)" = "#FC8D62",  # warm orange
      "Deviated (Discordant)" = "#D73027"   # strong warm red
    ),breaks=c("Identical","Near-identical","Deviated (Concordant)","Deviated (Discordant)"),
    guide = guide_legend(override.aes = list(size = 4))
  )+
  theme_classic() + 
  scale_x_continuous(expand = c(0.01,0),labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(expand = c(0.01,0),labels = scales::percent_format(accuracy = 1)) +
  theme(legend.title=element_blank(),
        axis.text=element_text(color="black"))+
  xlab("CGC2 percent protein length\ndifference from N2") +
  ylab("AF16 percent protein length\ndifference from N2")
#pleg <- cowplot::get_legend(plac1)


plac2 <- ggplot(orthos_3way) + 
  geom_point(aes(x=CGC2_acc-1,y=AF16_acc-1,color=class),size=0.7) + 
  coord_cartesian(xlim=c(-0.1,0.1),ylim=c(-0.1,0.1)) + 
  theme_classic() + 
  scale_color_manual(
    values = c(
      "Identical" = "blue",              # cool teal
      "Near-identical" = "#66C2A5",         # lighter cool teal
      "Deviated (Concordant)" = "#FC8D62",  # warm orange
      "Deviated (Discordant)" = "#D73027"   # strong warm red
    ),breaks=c("Identical","Near-identical","Deviated (Concordant)","Deviated (Discordant)")
  )+
  scale_x_continuous(expand = c(0.01,0),labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(expand = c(0.01,0),labels = scales::percent_format(accuracy = 1)) +
  theme(axis.line = element_blank(),
        panel.border= element_rect(fill=NA,linetype = "22",linewidth = 1.5),
        plot.margin = margin(t = 5.5, r = 8, b = 5.5, l = 5.5),
        axis.text=element_text(color="black")) +
  xlab("CGC2 percent protein length\ndifference from N2") +
  ylab("AF16 percent protein length\ndifference from N2")

                   
class_counts <- orthos_3way %>%
  dplyr::transmute(
    AF16_class = case_when(
      AF16_deviation == 0 & CGC2_deviation > 0.05 ~ "Identical",
      AF16_deviation <= 0.05 & CGC2_deviation > 0.05 ~ "Near-identical",
      AF16_deviation > 0.05  ~ "Deviated",
      TRUE ~ NA_character_
    ),
    CGC2_class = case_when(
      CGC2_deviation == 0 & AF16_deviation > 0.05 ~ "Identical",
      CGC2_deviation <= 0.05 & AF16_deviation > 0.05  ~ "Near-identical",
      CGC2_deviation > 0.05 ~ "Deviated",
      TRUE ~ NA_character_
    )
  ) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "strain",
    values_to = "class"
  ) %>%
  dplyr::mutate(
    strain = recode(strain,
                    AF16_class = "AF16",
                    CGC2_class = "CGC2"
    ),
    class = factor(class, levels = c("Identical", "Near-identical", "Deviated"))
  ) %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::count(strain, class)

class_count_plt <- ggplot(class_counts %>% dplyr::filter(class!="Deviated"), aes(x = class, y = n, fill = class)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.3,size=2) +
  ggh4x::facet_wrap2(~ strain,
                     strip = ggh4x::strip_themed(
                       background_x = ggh4x::elem_list_rect(
                         fill = c("#E6B3B3","#7570B3")
                       )
                     )
  ) +
  scale_fill_manual(
    values = c(
      "Identical" = "blue",
      "Near-identical" = "#66C2A5",
      "Deviated" = "red"
    )
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Number of\nSC orthologs"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid = element_blank(),
    axis.text=element_text(color="black")
  ) +
  scale_y_continuous(expand = c(0,0))+
  expand_limits(y = max((class_counts %>% dplyr::filter(class!="Deviated"))$n) * 1.08)


deviation_type_counts <- orthos_3way %>%
  filter(AF16_deviation > 0.05, CGC2_deviation > 0.05) %>%
  mutate(
    deviation_type = case_when(
      AF16_deviation == CGC2_deviation ~ "Concordant",
      TRUE ~ "Discordant"
    ),
    deviation_type = factor(deviation_type, levels = c("Concordant", "Discordant"))
  ) %>%
  count(deviation_type)

dev_count_plt <- ggplot(deviation_type_counts %>% dplyr::mutate(banner="Deviated"), aes(x = deviation_type, y = n, fill = deviation_type)) +
  geom_col() +
  ggh4x::facet_wrap2(~ banner,
                     strip = ggh4x::strip_themed(
                       background_x = ggh4x::elem_list_rect(
                         fill = c("grey")
                       )
                     )
  ) +
  geom_text(aes(label = n), vjust = -0.3,size=2) +
  scale_fill_manual(
    values = c(
      "Concordant" = "#FC8D62",
      "Discordant" = "#D73027"
    )
  ) +
  theme_minimal() +
  labs(
    x = NULL,
    y = ""
  )  +
  theme_bw()+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid = element_blank(),
    axis.text=element_text(color="black")
  ) +
  scale_y_continuous(expand = c(0,0))+
  expand_limits(y = max(deviation_type_counts$n) * 1.08)

ct_plts <- cowplot::plot_grid(class_count_plt,dev_count_plt,nrow=1,rel_widths = c(1.8,1),align = "h",axis = "tb",labels = c("c","d"))
plac_zoom_pad <- cowplot::plot_grid(NULL,plac2+theme(legend.position = 'none'),labels = c("b",""),nrow=1,rel_widths = c(0.07,1))
ct_zoom_plt <- cowplot::plot_grid(plac_zoom_pad,ct_plts,align = "v",axis = "lr",nrow=2,rel_heights = c(1.2,1))
acc_plt <- cowplot::plot_grid(cowplot::plot_grid(NULL,
                                                 plac1+theme(legend.position = 'inside',legend.position.inside = c(0.35,0.92)),
                                                 nrow=2,
                                                 rel_heights = c(0.04,1)),
                              ct_zoom_plt,
                              nrow=1,
                              rel_widths = c(1,0.9),
                              labels = c("a",""))
#acc_plt

ggsave(acc_plt,filename = "../../figures/Figure5_ACC.png",width = 7,height = 4.9,dpi = 600,device = 'png',bg = "white")



dev_bar <- orthos_3way %>%
  tidyr::pivot_longer(
    cols = c(AF16_deviation, CGC2_deviation),
    names_to = "sample",
    values_to = "deviation"
  ) %>%
  dplyr::mutate(
    sample = dplyr::recode(
      sample,
      AF16_deviation   = "AF16",
      CGC2_deviation   = "CGC2"
    ),
    dev_bin = dplyr::case_when(
      deviation == 0 ~ "0",
      deviation > 0    & deviation <= 0.05 ~ ">0-5%",
      deviation > 0.05 & deviation <= 0.1  ~ ">5-10%",
      deviation > 0.1  & deviation <= 0.25 ~ ">10-25%",
      deviation > 0.25 ~ ">25%",
      TRUE ~ NA_character_
    ),
    dev_bin = factor(
      dev_bin,
      levels = c("0", ">0-5%", ">5-10%", ">10-25%", ">25%")
    )
  ) %>%
  dplyr::filter(!is.na(dev_bin)) %>%
  dplyr::count(sample, dev_bin)


counts1 <-ggplot2::ggplot(dev_bar, ggplot2::aes(x = dev_bin, y = n, fill = sample)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::labs(
    x = "Deviation bin",
    y = "Count",
    fill = "Sample"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = n),
    position = ggplot2::position_dodge(width = 0.9),
    vjust = -0.2,
    angle = 60,
    hjust = 0,
    size = 3
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "CGC2" = "#7570B3",
      "AF16" = "#E6B3B3"
    )
  ) +
  ggplot2::theme_classic() +
  scale_y_continuous(expand= expansion(mult = c(0, 0.1), add = c(1, 0))) +
  ggtitle("")+
  theme(axis.text.x = element_text(angle=60,hjust=1),
        axis.text=element_text(color="black"),
        legend.position="inside",
        legend.position.inside = c(0.9,0.9),
        legend.title = element_blank())+
  xlab("Absolute percent protein sequence\nlength difference from N2")+
  ylab("Number of genes")
                   
ggsave(counts1,filename = "../../figures/supplementary/FigureS16_ACC_ct.png",width = 7,height = 5,dpi = 600,device = 'png',bg = "white")                   


orthos_miss_AF16 <- orthos_wlen %>%
  dplyr::filter(is.na(AF16) & !is.na(CGC2) & !is.na(QX1410) & is.na(N2))

orthos_miss_CGC2 <- orthos_wlen %>%
  dplyr::filter(!is.na(AF16) & is.na(CGC2) & !is.na(QX1410) & is.na(N2))

orthos_miss_QX1410 <- orthos_wlen %>%
  dplyr::filter(!is.na(AF16) & !is.na(CGC2) & is.na(QX1410) & !is.na(N2))

strains <- c("CGC2", "AF16", "QX1410")

summary_table <- lapply(strains, function(target) {
  
  others <- setdiff(strains, target)
  
  # target missing, other two C.b. present
  base_subset <- orthos_counts |>
    dplyr::filter(
      is.na(.data[[target]]),
      !is.na(.data[[others[1]]]),
      !is.na(.data[[others[2]]])
    )
  
  # target present, N2 present
  target_present_subset <- orthos_counts |>
    dplyr::filter(
      !is.na(.data[[target]]),
      !is.na(N2)
    )
  
  tibble::tibble(
    mising_in = target,
    missing_with_N2_ortholog = sum(!is.na(base_subset$N2)),
    #N2_missing = sum(is.na(base_subset$N2)),
    total_missing = nrow(base_subset),
    total_with_N2_ortholog = nrow(target_present_subset),
    total_with_N2_single_copy_ortholog = sum(!stringr::str_detect(target_present_subset$N2, ","))
  )
  
}) |>
  dplyr::bind_rows()

strains <- c("CGC2", "AF16", "QX1410")

summary_table <- lapply(strains, function(target) {
  
  others <- setdiff(strains, target)
  
  base_subset <- orthos_counts |>
    dplyr::filter(
      is.na(.data[[target]]),
      !is.na(.data[[others[1]]]),
      !is.na(.data[[others[2]]])
    )
  
  target_present_subset <- orthos_counts |>
    dplyr::filter(
      !is.na(.data[[target]]),
      !is.na(N2)
    )
  
  tibble::tibble(
    absent_in = target,
    absent_with_N2_ortholog = sum(!is.na(base_subset$N2)),
    absent_with_N2_single_copy_ortholog = sum(!is.na(base_subset$N2) &
                              !stringr::str_detect(base_subset$N2, ",")),
    total_absent = nrow(base_subset),
    total_with_N2_ortholog = nrow(target_present_subset),
    total_with_N2_single_copy_ortholog = sum(!stringr::str_detect(target_present_subset$N2, ","))
  )
  
}) |>
  dplyr::bind_rows()

cgc2af <- readr::read_tsv("../../processed_data/liftoff/CGC2toAF16/c_briggsae.AF16.liftoff.gff",col_names = c("seqid","source","type","start","end","score","strand","phase","attributes"))

cgc2af_tran <- cgc2af %>%
  dplyr::filter(type=="mRNA") %>%
  dplyr::select(attributes) %>%
  dplyr::mutate(attributes=gsub(";Parent=.*","",attributes)) %>%
  dplyr::mutate(attributes=gsub("ID=","",attributes)) %>%
  dplyr::filter(attributes%in%orthos_miss_AF16$CGC2)

af2cgc <- readr::read_tsv("../../processed_data/liftoff/AF16toCGC2/c_briggsae.CGC2.hifi.ONT.HiC.Feb2026.liftoff.gff",col_names = c("seqid","source","type","start","end","score","strand","phase","attributes"))

af2cgc_tran <- af2cgc %>%
  dplyr::filter(type=="mRNA") %>%
  dplyr::select(attributes) %>%
  dplyr::mutate(attributes=gsub(";Parent=.*","",attributes)) %>%
  dplyr::mutate(attributes=gsub("ID=","",attributes)) %>%
  dplyr::mutate(attributes=gsub(":","_",attributes)) %>%
  dplyr::filter(attributes%in%orthos_miss_CGC2$AF16)
