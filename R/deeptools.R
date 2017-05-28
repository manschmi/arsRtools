## Functions for parsing deeptools matrix object

#' load_deeptoolsmatrix
#'
#' Loads a deeptools matrix file.
#'
#' @param file the file to load.
#' @param na.omit Skip na values in the matrix.
#' @param na.fill Fills NAs with this value zero.
#' @param as.matrix Return matrix as matrix (default is FALSE, returns tidy data.frame)
#'
#' @details Loads a gzipped deeptools matrix file and collects this in a \code{tidy} data.frame. If as.matrix is true, the header information is added as attribute.
#'  The tidy dataframe contains the following columns
#'   . chr start end id name strand (first 6 columns of the bed file used for computeMatrix, content may differ depending on deeptools version used)
#'   . sample_name (name of the bigwig file used for that sample)
#'   . group if deeptools matrix contained groups (ie grouped bed files)
#'   . rel_pos relative position to reference-point ((ambigous for scale-regions mode ..))
#'   . value
#'
#' @return A tidy tibble or matrix.
#'
#' @examples
#'
#'  fname <- "/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/pA_seq/deeptools_out/Evgenia_PASseq_all_snRNA_joined.gz"
#'  df <- load_deeptoolsmatrix(fname)
#'  head(df)
#'  \dontrun{df %>% group_by(sample_name, rel_pos) %>% summarize(value=sum(value)) %>% ggplot(., aes(x=rel_pos, y=value, color=sample_name)) + geom_line()}
#'  df_withNA <- load_deeptoolsmatrix(fname, na.omit=FALSE)
#'  head(df_withNA)
#'  df_withNAas0 <- load_deeptoolsmatrix(fname, na.omit=FALSE, na.fill=0)
#'  head(df_withNAas0)
#'
#' @export
load_deeptoolsmatrix <- function(fname, na.omit=TRUE, na.fill=0, as.matrix=FALSE){
  header <- parse_matrix_header(fname)
  values <- read.table(gzfile(fname), header=F, skip=1)

  bed_info <- values[,1:6]
  msamples <- unlist(header['sample_labels'])
  sample_bounds <- as.numeric(unlist(header['sample_boundaries']))
  bin <- as.numeric(header['bin size'])
  body <- as.numeric(header['body'])
  up <- as.numeric(header['upstream'])
  unscaled_5 <- as.numeric(header['unscaled 5 prime'])
  unscaled_3 <- as.numeric(header['unscaled 3 prime'])
  dn <- as.numeric(header['downstream'])

  positions <- seq(-up,body+unscaled_5+unscaled_3+dn-1,bin)
  if (header['ref point'] == "null") {
    xlabels <- c()
    if(up >= bin) xlabels <- c(paste('TSS',seq(-up,-1,bin)))
    xlabels <- c(xlabels,'TSS')
    if(unscaled_5 > bin) xlabels <- c(xlabels, paste0('TSS +',seq(bin,unscaled_5-1,bin)))
    if(body > bin) xlabels <- c(xlabels, paste('body',seq(bin,body-1,bin)))
    if(unscaled_3 >= bin) xlabels <- c(xlabels, paste('TES',seq(-unscaled_3,-1,bin)))
    if(dn > bin) {
      xlabels <- c(xlabels,'TES')
      xlabels <- c(xlabels, paste0('TES +',seq(bin,dn-1,bin)))
    }
    if(length(positions)==1) xlabels <- 'body'
    else if(length(positions)==2) xlabels <- c('TSS', 'mid')
  } else {
    xlabels <- positions
  }

  mgroups <- unlist(header['group_labels'])
  group_bounds <- as.numeric(unlist(header['group_boundaries']))

  if (as.matrix) {
    attr(values, 'header') <- header
    bednames <- c('chr', 'start', 'end', 'id', 'name', 'strand')

    name_vec <- c()
    for(i in seq_along(msamples)) {
      name <- msamples[i]
      left <- sample_bounds[i] + 7
      right <- sample_bounds[i+1] + 6
      name_vec <- c(name_vec, paste0(name, xlabels))
    }
    colnames(values) <- c(bednames, name_vec)
    return(values)
  }

  dfs <- list()
  for(i in seq_along(msamples)) {
    name <- msamples[i]
    left <- sample_bounds[i] + 7
    right <- sample_bounds[i+1] + 6
    s_values <- cbind(bed_info, values[,left:right])
    colnames(s_values) <- c('chr', 'start', 'end', 'id', 'name', 'strand', xlabels)

    s_values$sample_name <- name
    s_values$group <- mgroups[1]
    for(g in seq_along(mgroups)) {
      upper <- group_bounds[g] + 1
      lower <- group_bounds[g+1]
      s_values$group[upper:lower] <- mgroups[g]
    }
    dfs[[name]] <- gather(s_values, rel_pos, value, -sample_name, -group, -c(1:6))
  }

  d <- bind_rows(dfs)
  if(na.omit) {
    d <- filter(d, !is.na(value))
  } else if (!missing(na.fill)) {
    fill_val <- as.numeric(na.fill)
    d <- mutate(d, value = ifelse(is.na(value), fill_val, value))
  }

  d %<>% mutate(sample_name = factor(sample_name))

  if (header['ref point'] != "null") {
    d %<>% mutate(rel_pos = as.numeric(rel_pos))
  }

  as_data_frame(d)
}


#' parse_matrix_header
#'
#' Parses a deeptools matrix file header.
#' @param file the file to get the header from.
#'
#' @details Loads a gzipped deeptools matrix file and collects this in a \code{tidy} data.frame.
parse_matrix_header <- function(fname) {

  header <-  scan(gzfile(fname), nlines=1, skip=0, what='character', sep='#') %>% sub('@\\{', ',', .) %>% sub('\\}', '', .) %>% sub('"', '', .)

  chunks <- strsplit(header, ',')

  meta <- c()
  keep_pasting <- FALSE
  for (c in chunks[[1]][2:length(chunks[[1]])]) {
    if (keep_pasting){
      meta[length(meta)] <- paste(meta[length(meta)], c, sep=';')
      if (grepl('\\]',c)) {
        keep_pasting <- FALSE
      }
      next()
    }
    if (grepl('\\[',c) && !grepl('\\]',c)) {
      keep_pasting <- TRUE
    }

    meta <- c(meta, c)
  }

  tag_id <- sapply(meta, function(m) strsplit(m, ':')[[1]][1])
  tag_vals <- lapply(meta, function(m) strsplit(m, ':')[[1]][2])

  names(tag_vals) <- tag_id

  tag_vals <- lapply(tag_vals, function(m) sub('\\[', '', m))
  tag_vals <- lapply(tag_vals, function(m) sub('\\]', '', m))
  tag_vals <- lapply(tag_vals, function(m) strsplit(m, ';')[[1]])

  tag_vals
}



#' plot_heatmap
#'
#' Plots tidy dataframe from a deeptools matrix file to a ggplot2 object.
#'
#' @param df the df to load.
#' @param sort_y how to sort y axis of the heatmap ('region_len', 'sum', 'mean', 'max' )
#' @param facets a string for facetting of the plot.
#'
#' @details Plots a df from a deeptools matrix file.
#'
#' @return An ggplot2 plot object.
#'
#' @examples
#'
#'  fname <- "/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/pA_seq/deeptools_out/Evgenia_PASseq_all_snRNA_joined.gz"
#'  df <- load_deeptoolsmatrix(fname)
#'  head(df)
#'  plot_heatmap(df, sort_y='sum', facets='.~sample_name')
#'
#'  #can be further modified
#'  p <- plot_heatmap(df, sort_y='sum', facets='.~sample_name')
#'  p + ylab('')
#'  p + scale_fill_gradient2(low='white', high='darkred')
#'  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#'  p + scale_fill_gradient2(low='white', high='darkred', name ='log10(norm. reads)')
#'
#' @export
plot_heatmap <- function( df, sort_y='sum', facets=NULL ) {

  if ( missing(sort_y) ) {
    y_rank <- df %>% distinct(id) %>% mutate(sort_value = 1:nrow(.))
  } else if (sort_y == 'sum') {
    y_rank <- df %>% group_by(id) %>% summarize(sort_value = sum(value)) %>% select(id, sort_value)
  } else if (sort_y == 'mean') {
    y_rank <- df %>% group_by(id) %>% summarize(sort_value = mean(value)) %>% select(id, sort_value)
  } else if (sort_y == 'max') {
    y_rank <- df %>% group_by(id) %>% summarize(sort_value = max(value)) %>% select(id, sort_value)
  } else if (sort_y == 'region_length') {
    y_rank <- df %>% distinct(id, .keep_all=T) %>% mutate(sort_value = end-start)
  }

  y_rank %<>% mutate(yrank = rank(sort_value, ties.method='first'))

  p <- left_join(df,y_rank)  %>%
    ggplot(., aes(x=rel_pos, y=yrank, fill=log10(value))) +
    geom_raster(interpolate=FALSE) +
    theme_bw()

  if (facets) {
    p <- p +
      facet_grid(facets)
  }

  p
}



#' Sort Matrix
#'
#' Sorts a deeptools matrix file using a bed file.
#'
#' @param matrix_file filename for the matrix (gzipped).
#' @param bed_file filename for the bed file used for sorting
#'
#' @details Sorts the rows of the matrix without header.
#'
#' @return The sorted matrix.
#'
#' @examples
#'
#'  matrix_name <- '/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/RNAseq_deeptools_out/eRNA_tsspm5000_sorted_by_INT_binding_joined_sensitivity.gz'
#'  bed_name <- '/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Shiekkatar_INT_ChIP/eRNAs_sorted_by_INTS_binding.bed'
#'
#'  m_sorted <- sort_matrix(matrix_name, bed_name)
#'
#'  write.table(m_sorted, file='/Users/schmidm/Documents/other_people_to_and_from/ClaudiaI/Effie_RNAseq_bws/RNAseq_deeptools_out/eRNA_tsspm5000_sorted_by_INT_binding_joined_sensitivity_sorted', row.names = FALSE, col.names = FALSE, quote = F, dec = '.')
#'
#'
#' @export
sort_matrix <- function( matrix_file, bed_file ) {

  m <- load_deeptoolsmatrix(matrix_file, as.matrix = TRUE) %>%
    tbl_df

  b <- read.table(bed_file)[,1:6]
  colnames(b) <- c('chr', 'start', 'end', 'id', 'name', 'strand')
  b %<>% tbl_df

  semi_join(m, b)
}
