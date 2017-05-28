#' Load RNAseq Matrices
#'
#' Load deeptools matrices for RNAseq.
#'
#' @param iasillo_matrix Filename for iasillo gzipped matrix.
#' @param meola_matrix Filename for Meola gzipped matrix.
#' @param sample_names List of sample_names to include.
#'
#' @details Probably only useful in ARS2 paper context.
#'
#' @return combined matrix with extra column "study"
#'
#' @examples
#'
#'  deeptools_file_iasillo <- '/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/RNAseq_deeptools_out/snRNA_MS31_TES_sensitivity_joined_sensitivity.gz'
#'  deeptools_file_meola <- '/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/RNAseq_deeptools_out_meola/snRNA_MS31_TES_sensitivity_joined_sensitivity.gz'
#'  bin_values <- load_RNAseq_matrices(deeptools_file_iasillo, deeptools_file_meola)
#'  data('RNAseq_value_heatmap_theme', package='arsRtools')
#'
#'  ggplot(filter(bin_values, value>0), aes(x=rel_pos, y=id, fill=log2(value))) +
#'    geom_raster(interpolate = FALSE) +
#'    facet_grid(.~siRNA) +
#'    RNAseq_value_heatmap_theme
#
#'
#' @export
load_RNAseq_matrices <- function(deeptools_file_iasillo,
                                deeptools_file_meola,
                                sample_names = c('EGFP', 'eGFP', 'Ars2', 'Z18', 'RRP40')) {

  iasillo <- load_deeptoolsmatrix(sub('_sensitivity.gz', '.gz', deeptools_file_iasillo),
                                  na.omit = FALSE,
                                  na.fill = 0) %>%
    mutate(study='Iasillo')

  meola <- load_deeptoolsmatrix(sub('_sensitivity.gz', '.gz', deeptools_file_meola),
                                na.omit = FALSE,
                                na.fill = 0) %>%
    mutate(study='Meola')

  bind_rows(iasillo, meola) %>%
    mutate(sample_name = gsub('-', '_', sample_name) %>%
             sub('_tot', '', .) %>%
             sub('01_', '', .)) %>%
    separate(sample_name, c('siRNA', 'rep'), sep='_', extra='warn' ) %>%
    filter(siRNA %in% sample_names)

}



#' RNAseq Sensitivity Matrix
#'
#' Converts tidy matrix to tidy sensitivity matrix for RNAseq.
#'
#' @param bin_values A tidy deeptools matrix containing value in single bins.
#'
#' @details Only useful in ARS2 paper context. matrix needs to contain a column study
#'
#' @return Sensitivity matrix
#'
#' @examples
#'
#'  deeptools_file_iasillo <- '/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/RNAseq_deeptools_out/snRNA_MS31_TES_sensitivity_joined_sensitivity.gz'
#'  deeptools_file_meola <- '/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/RNAseq_deeptools_out_meola/snRNA_MS31_TES_sensitivity_joined_sensitivity.gz'
#'  bin_values <- load_RNAseq_matrices(deeptools_file_iasillo, deeptools_file_meola)
#'  bin_sensitivities <- RNAseq_sensitivity_matrix(bin_values)
#'
#'  data('RNAseq_logFC_heatmap_theme', package='arsRtools')
#'
#'  p <- ggplot(filter(bin_sensitivities, value>0), aes(x=rel_pos, y=id, fill=log2(value))) +
#'    geom_raster(interpolate = FALSE) +
#'    facet_grid(.~siRNA) +
#'    RNAseq_logFC_heatmap_theme
#'
#'  #default plot
#'  p
#'
#'  # can still override some plot options
#'  p + theme(axis.text.y=element_blank())
#'
#' @export
RNAseq_sensitivity_matrix <- function(bin_values) {

    min_values <- group_by(bin_values, study) %>%
      filter(value > 0) %>%
      summarise(min_val = min(value))

    mean_bin_values <- left_join(bin_values, min_values) %>%
      mutate(value = value + min_val) %>%
      group_by(id, rel_pos, siRNA, study) %>%
      summarize(value = mean(value)) %>%
      ungroup

    egfp <- filter(mean_bin_values, grepl('GFP', siRNA)) %>%
      dplyr::select(-siRNA) %>%
      dplyr::rename(egfp_value = value)

    kd <- filter(mean_bin_values, !grepl('GFP', siRNA))

    left_join(kd, egfp) %>%
      mutate(value = value/egfp_value,
             siRNA = factor(siRNA, levels=c('Ars2',
                                            'Z18',
                                            'RRP40')))
}
