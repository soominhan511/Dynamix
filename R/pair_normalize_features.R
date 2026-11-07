# Script : pair_normalize_features.R
#   This script contains functions for normalizing GC content and mappability for each interaction (i,j)

#' @importFrom dplyr %>% mutate filter pull select case_when rename
#' @importFrom tidyr separate crossing
#' @importFrom rlang sym
#' @importFrom data.table fwrite as.data.table tstrsplit `:=`
#' @importFrom tibble as_tibble

# Add Hi-C Interaction Counts to Pair Features


normalize_pair_features <- function(pairs_df, col_i, col_j) {
  
  sym_i <- rlang::sym(col_i)
  sym_j <- rlang::sym(col_j)
  norm_col_name <- paste0(col_i, "_norm") 
  
  pairs_df <- pairs_df %>%
    mutate(
      is_na_exist = is.na(!!sym_i) | is.na(!!sym_j),
      is_zero_exist = !is_na_exist & (!!sym_i == 0 | !!sym_j == 0),
      is_valid_for_log = !is_na_exist & (!!sym_i != 0 & !!sym_j != 0)
    )
  
  valid_vals_i <- pairs_df %>% filter(is_valid_for_log) %>% pull(!!sym_i)
  valid_vals_j <- pairs_df %>% filter(is_valid_for_log) %>% pull(!!sym_j)
  
  valid_logs <- log(valid_vals_i * valid_vals_j)
  mean_val <- mean(valid_logs, na.rm = TRUE)
  sd_val <- sd(valid_logs, na.rm = TRUE)
  
  # normalize pair values 
  pairs_df <- pairs_df %>%
    mutate(
      !!norm_col_name := case_when(
        is_na_exist ~ NA_real_,     
        is_zero_exist ~ 0,          # zero if col_i=0 or col_j=0 
        is_valid_for_log ~ (log((!!sym_i) * (!!sym_j)) - mean_val) / sd_val,        # normalize pair values with both non-zeros
        TRUE ~ NA_real_
      )
    ) %>%
    
    select(-is_na_exist, -is_zero_exist, -is_valid_for_log)
  
  return(pairs_df)
}


#' @export
get_pair_features_df <- function(features_df, output_path, max_distance = 2000000, resolution = 10000) {
  
  dt <- as.data.table(features_df)
  dt[, c("chr", "start", "end") := tstrsplit(bins, "-", fixed = TRUE)]
  dt[, start := as.numeric(start)]
  dt[, end := as.numeric(end)]
  chr_name <- dt$chr[1]

  dt_i <- dt[, .(chr_i = chr, start_i = start, start_i_rep = start, end_i = end, gc_i = gc, map_i = map, len_i = len)]
  dt_j <- dt[, .(chr_j = chr, start_j = start, start_j_rep = start, end_j = end, gc_j = gc, map_j = map, len_j = len)]
  dt_i[, start_limit_i := start_i + max_distance]
  
  # pair features
  pairs_dt <- dt_j[dt_i, 
                   on = .(start_j >= start_i, start_j <= start_limit_i), 
                   nomatch = NULL, 
                   allow.cartesian = TRUE]
  pairs_df <- as_tibble(pairs_dt[, .(chr_i, chr_j, start_i_rep, end_i, gc_i, map_i, len_i, 
                                     start_j_rep, end_j, gc_j, map_j, len_j)])
  pairs_df <- pairs_df %>%
    rename(
      start_i = start_i_rep,
      start_j = start_j_rep
    )
  
  # normalize gc and map 
  pairs_df <- pairs_df %>%
    normalize_pair_features("gc_i", "gc_j") %>%   
    normalize_pair_features("map_i", "map_j")    

  result_df <- pairs_df %>%
    filter(len_i != 0 & len_j != 0) %>%
    mutate(
        chr = chr_name,
        idx1 = as.integer(start_i / resolution),
        idx2 = as.integer(start_j / resolution),
        D = abs(start_i - start_j)
    ) %>%
    select(
        chr, 
        start1 = start_i, 
        start2 = start_j, 
        end1 = end_i, 
        end2 = end_j, 
        gc = gc_i_norm, 
        map = map_i_norm,
        idx1,
        idx2,
        D
    )

  # save features to file 
  if (is.null(output_path)) {  output_path <- getwd()  }
  if (!dir.exists(output_path)) {
      dir.create(output_path, showWarnings = FALSE, recursive = TRUE, mode = "0777")
  }
  file_name <- file.path(output_path, paste0(chr_name, '_pair_features_tmp.csv'))
  data.table::fwrite(result_df, file_name, row.names = FALSE, quote = FALSE, sep = ",")

  return(result_df)
  
}

#' @export
add_counts_to_features <- function(feature_path, hic_path, chrom, resolution) {
    
    py_script_path <- system.file("python", "call_stripes.py", package = "Dynamix")
    
    # source python functions 
    reticulate::source_python(py_script_path)
    
    # save csv feature file added hic counts 
    add_hic_counts_to_features(
        feature_file = as.character(feature_path),
        hic_file     = as.character(hic_path),
        chrom        = as.character(chrom),
        resolution   = as.integer(resolution)
    )
    
}

