# Script : run_Dynamix.R
#   This script contains the function Call_Stripes, which identifies chromatin stripes 
#   from the signal probability matrix using a Bayesian NDLM.
#   The signal probability matrix is constructed by ZINB+NB mixture EM algorithm

#' @importFrom readr read_csv
#' @export
Call_Stripes <- function(hic_file=hic_file, wg_file=wg_file, output_path=output_path, chr=chr, resolution=1e4){

    if (!dir.exists(output_path)) {
        dir.create(output_path, recursive = TRUE)
    }

    if (is.null(wg_file) || !file.exists(wg_file)) {
        message("Please specify the path of mappability bigWig file using wg_file = 'your_file_path'.")
        stop("Example of wg_file: wgEncodeCrgMapabilityAlign50mer.bigWig (you can download it from the UCSC website).")
    }
    
    # extract bin-level GC ratio and Mappability
    features_df <- construct_features(
        gen          = "Hsapiens", 
        gen_ver      = "hg19",
        sig          = c("GATC", "GANTC"),   
        bin_type     = "Bins-uniform",       
        binsize      = resolution,           
        wg_file      = wg_file,                    
        feature_type = "RE-based",           
        chrs         = c(chr)                
    )

    # compute the normalized sum of log-transformed GC, Map values for bin pairs
    pair_df <- get_pair_features_df(features_df, output_path=output_path, resolution=resolution)
    feature_path_tmp <- file.path(output_path, paste0(chr, '_pair_features_tmp.csv'))

    # add hic interaction counts to bin pairs
    add_counts_to_features(feature_path_tmp, hic_file, chr, resolution)
    feature_path <- file.path(output_path, paste0(chr, '_pair_features.csv'))
    feature_df <- read_csv(feature_path, show_col_types=FALSE)

    # calculate signal probabilies per bin pair 
    # via Zero-Infated Negative Binomial + Negative Binomial mixture EM algorithm
    em_path <- file.path(output_path, paste0(chr, '_stats/'))
    print('Calculating signal probabilities for each interaction via ZINB-NB EM mixture ...')
    get_signal_p(feature_df, em_path)
    
    # call stripes from signal probability matrix
    # via Bayesian inference in Normal Dynamic Linear Model
    signal_p_path <- paste0(em_path , paste0(chr, '_pair_features_stats.csv'))
    call_stripes_by_ndlm(signal_p_path, chr, resolution=resolution)

}