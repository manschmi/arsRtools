#' simplify_ChIP_name
#'
#' Simplifies names of ChIP samples used for Claudias analysis.
#'
#' @param name The name to simplify.
#' @param regex_list List of substrings to replace, (ie list of regular expression).
#' @param replace_list List of substrings that replace regex_list hits, usually empty.
#'
#' @details Simply uses sub(regex, replace) for all
#'
#' @return Simplified name
#'
#' @examples
#'
#'  name <- "/WIGfs_Hela-H9_WT_siARS2_Input_CI622_bw2align_hg19_mergedReads_elPairsAndEst150_AThr12_bin50"
#'  simplify_ChIP_name(name)
#'  #
#'
#' @export
simplify_ChIP_name <- function(name,
                                 regex_list=c('Pol_II_N20_',
                                              '_bw2align_hg19_mergedReads_elPairsAndEst.*_AThr.*_bin50',
                                              'WIGfs_Hela-H9_WT_',
                                              '_Scaled_BGSub'),
                                 replace_list=c('PolII_',
                                                '',
                                                '',
                                                '')) {
  for(i in seq_along(regex_list)) {
    name %<>% sub(regex_list[i], replace_list[i], .)
  }

  name
}



#' RefSeqGRCh37.3 Genes Bed file
#'
#' @export
refseq_genes_anno <- function() {
  fname <- '/Users/schmidm/Documents/genomewide_datasets/annotations/hg19/RefSeqGRCh37.3/ref_GRCh37.p5_top_level_genes.bed'
  print(paste0('loading Refseq_genes_anno from ', fname))
  read.table(fname, col.names=c('chr', 'start', 'end', 'name', 'pad', 'strand'))
}



#' Load ChIP
#'
#' Loads a deeptools matrix of ChIP (Pasha scaled or unscaled).
#'
#' @param matrix Filename for deeptools_matrix.
#' @param sample_table A dataframe mapping ChIP to Input.
#'
#' @details Loads a deeptools matrix. And adds info from sample_table
#'   and scales according to scales.
#'
#'   sample_table requires named columns:
#'      . sample_name
#'      . siRNA (siFFL, siARS2, siZ18 or siRRP40)
#'      . ab (usually PolII or Input)
#'      . run (CI62, CI63, MA72_rrp40, MA72_z18, MA73_rrp40, MA73_z18)
#'      . used_input (a sample_name, blank for Inputs, ie siFFL_Input_CI623)
#'
#'
#' @return A tidy tibble with columns id, siRNA, ab, run, rel_pos, value
#'
#' @examples
#'
#'  fname <- system.file("extdata", "PolII_snRNA_TES_matrix.gz", package = "arsRtools")
#'  data(chip_samples_table, package='arsRtools')
#'
#'  df <- load_chip_matrix(fname, samples_table)
#'
#' @export
load_chip_matrix <- function(fname, sample_table) {

  print('parsing matrix')
  cnts <- load_deeptoolsmatrix(fname)

  print('merging sample_table into matrix')
  cnts %<>% mutate(sample_name = plyr::mapvalues(sample_name,
                                                 levels(sample_name),
                                                 simplify_ChIP_name(levels(sample_name)))) %>%
    left_join(., sample_table) %>%
    mutate(siRNA = factor(siRNA, levels=c('siARS2', 'siFFL', 'siRRP40', 'siZ18')),
           run = factor(run, levels=c('CI62', 'CI63',
                                      'MA72_z18', 'MA73_z18',
                                      'MA72_rrp40', 'MA73_rrp40')))
}


#' Scale ChIP
#'
#' Applies manual scaling to a loaded deeptools chip matrix.
#'
#' @param matrix Loaded deeptools matrix with sample_names assigned in tidy format.
#' @param scales A dataframe containing scales to apply.
#' @param sample_table A dataframe mapping ChIP to Input and scales.
#' @param scale_to_peak Scale each sample to peak instead of "scales". Default
#'   is FALSE.
#'
#' @details Assigns Input to IP using sample_table
#'   and scales according to scales.
#'   LOGIC:
#'     scaled_input = Input * input_scale_factor
#'     BGSub_ChIP = ChIP - scaled_input
#'     scaled_BGSub_ChIP = chip_scale_factor * BGSub_ChIP
#'
#'    scales requires named columns:
#'      . sample_name (should matches sample_names in the loaded matrix.)
#'      . siRNA
#'      . ab (usually PolII or Input)
#'      . run
#'      . input_scale_factor
#'      . chip_scale_factor
#'
#'   sample_table requires named columns:
#'      . sample_name
#'      . siRNA (siFFL, siARS2, siZ18 or siRRP40)
#'      . ab (usually PolII or Input)
#'      . run (CI62, CI63, MA72_rrp40, MA72_z18, MA73_rrp40, MA73_z18)
#'      . used_input (a sample_name, blank for Inputs, ie siFFL_Input_CI623)
#'
#'   scale_to_peak: Uses max of BGSub_ChIP as chip_scale_factor instead of the
#'     value from "scales". Note: still requires input_scale_factor from scales!
#'
#'   scale_using_values: Uses mean values of another matrix (ie unscaled gene
#'   body, 1 or more bins to determine chip_scale_factor instead of the value
#'   from "scales". Note: still requires input_scale_factor from scales!
#'
#' @return A tidy tibble with columns id, siRNA, ab, run, rel_pos, scaled_BGSub_ChIP.
#'
#' @examples
#'
#'  fname <- system.file("extdata", "PolII_snRNA_TES_matrix.gz", package = "arsRtools")
#'  data(Refseq_body_scales, package='arsRtools')
#'  data(chip_samples_table, package='arsRtools')
#'
#'  df <- load_chip_matrix(fname, chip_samples_table)
#'  df_scaled <- scale_chip_matrix(df, Refseq_body_scales)
#'
#'  \dontrun{ df_scaled %>% group_by(run, siRNA, main_run, rel_pos) %>%
#'  summarize(mean_value = mean(scaled_BGSub_ChIP)) %>%ggplot(., aes(x=rel_pos,
#'  y=mean_value, color=siRNA)) + geom_line() + facet_grid(.~run) }
#'
#'  df_scaled2 <- scale_chip_matrix(df, Refseq_body_scales, scale_to_peak=TRUE)
#'
#'  \dontrun{ df2 %>% group_by(run, siRNA, main_run, rel_pos) %>%
#'  summarize(mean_value = mean(scaled_BGSub_ChIP)) %>%ggplot(., aes(x=rel_pos,
#'  y=mean_value, color=siRNA)) + geom_line() + facet_grid(.~run) }
#'
#' @export
scale_chip_matrix <- function(cnts,
                              scales,
                              scale_to_peak=FALSE,
                              scale_using_values=NULL) {

  print('finding ChIPs and Inputs')

  inputs <- dplyr::filter(cnts, ab=='Input') %>%
    dplyr::select(id, rel_pos, sample_name, value) %>%
    dplyr::rename(used_input = sample_name, Input = value)

  ips <- dplyr::filter(cnts, ab=='PolII') %>%
    dplyr::select(id, siRNA, ab, main_run, run, rel_pos, sample_name, used_input, value) %>%
    dplyr::rename(ChIP = value)

  print('assigning inputs to ChIP and applying scaling')


  d <- left_join(ips, inputs) %>%
    dplyr::select(id, siRNA, ab, main_run, run, rel_pos, ChIP, Input) %>%
    left_join(., scales) %>%
    mutate(scaled_input = Input*input_scale_factor,
           BGSub_ChIP = ChIP-scaled_input)

  if (scale_to_peak) {
    print('scaling to highest point in metagene profile')
    internal_scale <- d %>%
      dplyr::select(-chip_scale_factor) %>%
      group_by(siRNA, ab, main_run, run, rel_pos) %>%
      summarize(mean_BGSub_ChIP = mean(BGSub_ChIP)) %>%
      group_by(siRNA, ab, main_run, run) %>%
      summarize(chip_scale_factor = 1/max(mean_BGSub_ChIP))

    d %<>%
      dplyr::select(-chip_scale_factor) %>%
      left_join(., internal_scale)

  } else if ( !is.null(scale_using_values) ) {
    print('scaling to external provided scale')
    external_scale <- scale_using_values %>%
      group_by(id, siRNA, ab, main_run, run) %>%
      summarize(mean_value = mean(value)) %>%
      group_by(siRNA, ab, main_run, run) %>%
      summarize(chip_scale_factor = 1/mean(mean_value))

    d %<>% dplyr::select(-chip_scale_factor) %>%
      left_join(., external_scale)
  }

  d %>%
    mutate(scaled_BGSub_ChIP = BGSub_ChIP * chip_scale_factor,
           siRNA = factor(siRNA, levels=c('siARS2', 'siFFL', 'siRRP40', 'siZ18')),
           run = factor(run, levels=c('CI62', 'CI63',
                                      'MA72_z18', 'MA73_z18',
                                      'MA72_rrp40', 'MA73_rrp40')),
           main_run = factor(main_run, levels=c('CI siARS2',
                                                'MA siZ18',
                                                'MA siRRP40'))) %>%
    dplyr::select(id, main_run, run, siRNA, rel_pos, scaled_BGSub_ChIP)

}


#' Load and Manual Scale ChIP
#'
#' Loads a deeptools matrix of unscaled IP and inputs and applies a manual
#' scaling.
#'
#' @param matrix Filename for deeptools_matrix.
#' @param scales A dataframe containing scales to apply.
#' @param sample_table A dataframe mapping ChIP to Input and scales.
#' @param scale_to_peak Scale each sample to peak instead of "scales". Default
#'   is FALSE.
#'
#' @details Loads the deeptools matrix. Assigns Input to IP using sample_table
#'   and scales according to scales.
#'   LOGIC:
#'     scaled_input = Input * input_scale_factor
#'     BGSub_ChIP = ChIP - scaled_input
#'     scaled_BGSub_ChIP = chip_scale_factor * BGSub_ChIP
#'
#'    scales requires named columns:
#'      . sample_name (should matches sample_names in the loaded matrix.)
#'      . siRNA
#'      . ab (usually PolII or Input)
#'      . run
#'      . input_scale_factor
#'      . chip_scale_factor
#'
#'   sample_table requires named columns:
#'      . sample_name
#'      . siRNA (siFFL, siARS2, siZ18 or siRRP40)
#'      . ab (usually PolII or Input)
#'      . run (CI62, CI63, MA72_rrp40, MA72_z18, MA73_rrp40, MA73_z18)
#'      . used_input (a sample_name, blank for Inputs, ie siFFL_Input_CI623)
#'
#'   scale_to_peak: Uses max of BGSub_ChIP as chip_scale_factor instead of the
#'     value from "scales"
#'
#' @return A tidy tibble with columns id, siRNA, ab, run, rel_pos, scaled_BGSub_ChIP.
#'
#' @examples
#'
#'  fname <- system.file("extdata", "PolII_snRNA_TES_matrix.gz", package = "arsRtools")
#'  data(Refseq_body_scales, package='arsRtools')
#'  data(chip_samples_table, package='arsRtools')
#'
#'  df <- load_and_scale_chip_matrix(fname, Refseq_body_scales, chip_samples_table)
#'
#'  \dontrun{ df %>% group_by(run, siRNA, main_run, rel_pos) %>%
#'  summarize(mean_value = mean(scaled_BGSub_ChIP)) %>%ggplot(., aes(x=rel_pos,
#'  y=mean_value, color=siRNA)) + geom_line() + facet_grid(.~run) }
#'
#'  df2 <- load_and_scale_chip_matrix(fname, Refseq_body_scales, chip_samples_table, scale_to_peak=TRUE)
#'
#'  \dontrun{ df2 %>% group_by(run, siRNA, main_run, rel_pos) %>%
#'  summarize(mean_value = mean(scaled_BGSub_ChIP)) %>%ggplot(., aes(x=rel_pos,
#'  y=mean_value, color=siRNA)) + geom_line() + facet_grid(.~run) }
#'
#' @export
load_and_scale_chip_matrix <- function(fname,
                                       scales,
                                       sample_table,
                                       scale_to_peak=FALSE) {

  cnts <- load_chip_matrix(fname, sample_table)

  scale_chip_matrix(cnts, scales)

}
