#### GENERAL TABLES ####


#' ChIP samples table.
#'
#' Information for the PolII and Input ChIP samples.
#'
#' @format A data frame with named columns:
#' \describe{
#'   \item{sample_name}{style: siFFL_Input_CI623, all libraries}
#'   \item{used_input}{style: siFFL_Input_CI623, NA for Input samplenames}
#'   \item{siRNA}{...}
#'   \item{ab}{Input, PolII}
#'   \item{run}{CI62, CI63, MA72_z18, MA73_z18, MA72_rrp40, MA73_rrp40}
#'   \item{main_run}{CI siARS2, MA72 siZ18, MA siRRP40}
#' }
"chip_samples_table"



#### CHIP SCALINGS ####


#' Scaling using 2nd half of gene body
#'
#' Scaling for ChIP samples.
#'
#' @format A data frame with named columns:
#' \describe{
#'   \item{sample_name}{style: siFFL_Input_CI623, all libraries}
#'   \item{used_input}{style: siFFL_Input_CI623, NA for Input samplenames}
#'   \item{siRNA}{...}
#'   \item{ab}{Input, PolII, Tyr1P, Ser2P or Ser5P}
#'   \item{run}{CI62, CI63, MA72_z18, MA73_z18, MA72_rrp40, MA73_rrp40}
#'   \item{input_scale_factor}{multiply input with that number}
#'   \item{chip_scale_factor}{multiply (ChIP - Input*input_scale_factor) with that number}
#' }
"Refseq_pc_g2k_q1_second_half_body_scales"



#### ANNOTATIONS ####
## useless, should be removed ...

#' RefSeq GRCh37.3 Genes.
#'
#' From 'annotations/hg19/RefSeqGRCh37.3/ref_GRCh37.p5_top_level_genes.bed'
#'
#' @format A data frame with 6 named columns:
#' 'chr', 'start', 'end', 'name', 'score', 'strand'
#' \describe{
#'   \item{name}{style: LOC...; MIR1302-2; FAM138A etc}
#' }
"refseq_genes_anno"



#### GGPLOT2 THEMES ####


#' Theme for doing ggplot2 metagene lineplots.
#'
#' list(
#'   +         theme_bw(),
#'   +         theme(axis.text.x = element_text(angle=45, hjust=1),
#'                   +               panel.grid.major = element_blank(),
#'                   +               panel.grid.minor = element_blank(),
#'                   +               panel.background = element_blank()),
#'   +         scale_color_brewer(palette = 'Set1'),
#'   +         scale_fill_brewer(palette = 'Set1'),
#'   +         geom_vline(xintercept = 0,linetype=3)
#'   +     )
#'
"lineplot_theme"



#' Theme for doing ggplot2 RNAseq log2FC metagene lineplots.
#'
#' list(
#'   +     geom_vline(xintercept = 0, linetype=2),
#'   +     theme_bw(),
#'   +     scale_colour_manual(values = c(brewer.pal(6, "Set1")[c(1,4,3)])),
#'   +     scale_fill_manual(values = c(brewer.pal(6, "Set1")[c(1,4,3)])),
#'   +     theme(axis.text.x = element_text(angle=45, hjust=1, size=24),
#'               +           axis.text.y = element_text(size=24),
#'               +           panel.grid.minor.x = element_blank(),
#'               +           panel.grid.minor.y = element_blank(),
#'               +           legend.text = element_text(size=24),
#'               +           legend.title = element_text(size=24))
#'
"RNAseq_lineplot_theme"



#' Theme for doing ggplot2 TIFseq events lineplots.
#'
#' list(
#'   +     geom_vline(xintercept = 0, linetype=2),
#'   +     theme_bw(),
#'   +     scale_colour_manual(values = c(brewer.pal(5, "Set1")[c(5,2)])),
#'   +     theme(axis.text.x = element_text(angle=45, hjust=1, size=24),
#'               +           axis.text.y = element_text(size=24),
#'               +           panel.grid.minor.x = element_blank(),
#'               +           panel.grid.minor.y = element_blank(),
#'               +           legend.text = element_text(size=24),
#'               +           legend.title = element_text(size=24))
#'
"TIFseq_lineplot_theme"



#' Theme for doing ggplot2 RNAseq raw value heatmaps.
#'
#' heatmap_theme <- list(theme_bw(),
#' scale_fill_gradient2(low='white', high='navy'),
#' theme(panel.grid.major = element_blank(),
#' panel.grid.minor = element_blank(),
#' panel.background = element_blank(),
#' axis.text.x = element_text(angle=45, hjust=1),
#' axis.text.y = element_text(size=6)))
#'
"RNAseq_value_heatmap_theme"



#' Theme for doing ggplot2 RNAseq log2FC value heatmaps.
#'
#' heatmap_theme <- list(theme_bw(),
#' scale_fill_gradient2(low="firebrick3", mid='white', high='navy'),
#' theme(panel.grid.major = element_blank(),
#' panel.grid.minor = element_blank(),
#' panel.background = element_blank(),
#' axis.text.x = element_text(angle=45, hjust=1),
#' axis.title.y = element_blank(),
#' axis.text.y = element_text(size=6))),
#' scale_y_discrete(expand=c(0,0))
#'
"RNAseq_logFC_heatmap_theme"



#### RESULTS ####

#' DESeq2 Results for intronic counts.
#'
#' Introns from Refseq GRCH37.3
#'
#' @format A data frame with 4 named columns:
#' 'id', 'log2FoldChange', 'padj', 'main_run'
#' \describe{
#'   \item{main_run}{style: CI siARS2; MA siZ18; MA siRRP40}
#' }
"Refseq_introns_DESeq2_results"
