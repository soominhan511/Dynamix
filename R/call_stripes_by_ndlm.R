# Script : call_stripes_by_ndlm.R
#   This script contains functions to identify stripe anchors and length ranges
#   from signal probability matrix using a Bayesian NDLM

#' @export
save_ndlm_plot <- function(ts_data, results_smoothed_mean, ci_smoothed, direction, output_path, chrom, resolution){
    
    # genomic position index for x-axis
    T = length(ts_data)
    chrom_coords <- (0:(T-1)) * resolution 
    x_axis_vals <- chrom_coords / 1e3   # kbp unit
    xlab_text   <- "Genomic Position (kbp)"

    if (direction==1){
        dir = '5_prime'
    }else if (direction==2){
        dir = '3_prime'
    }
    filename <- file.path(output_path, paste0(chrom,"_rowsum_smoothed_ndlm_", dir, ".png"))
    png(filename, width = 3600, height = 900, res = 300)

    plot(x_axis_vals, ts_data, main = paste0(chrom," Row Sums NDLM"), type = 'n', 
        xlab = xlab_text, ylab = "Row sums (Signal probability)", 
        cex.main = 1.5, cex.lab = 1.2)
    
    # Observed row sums
    points(x_axis_vals, ts_data, pch = 20, cex = 0.5, col = "black")

    # Smoothed row sums and 95% CI 
    lines(x_axis_vals, results_smoothed_mean, type = 'l', col = 'blue', lwd = 2)
    lines(x_axis_vals, ci_smoothed[,1], type = 'l', col = 'blue', lty = 2)
    lines(x_axis_vals, ci_smoothed[,2], type = 'l', col = 'blue', lty = 2)

    # legend
    legend('topleft', 
           legend = c("Observed row sums", "Smoothed row sums (NDLM)", "95% CI"),
           col = c("black", "blue", "blue"), 
           pch = c(20, NA, NA),   
           lty = c(NA, 1, 2),   
           lwd = c(NA, 2, 1),
           cex = 0.6,          
           pt.cex = 0.6,        
           inset = c(0.01, 0.01), 
           bty = "n")           

    dev.off()
}

#' @export
fit_NDLM_on_rowsum <- function(row_sum_vector, chrom, smoothing_factor=100, var_prior=1e4, resolution=1e4, output_path){
    
    row_sums <- as.numeric(row_sum_vector)
    T <- length(row_sums)
    ts_data=row_sums[1:T]
    data <- list(yt = ts_data)

    # ratio of system variance to observational variance (smaller w_v_ratio results in heavier smoothing)
    w_v_ratio <- 1/smoothing_factor

    # Bayesian inference in Normal Dynamic Linear Model
    # y_t = F_t * theta_t + nu_t,          nu_t ~ N(0, V)
    # theta_t = G_t * theta_{t-1} + omega_t, omega_t ~ N(0, W)
    # set up matrices
    FF <- as.matrix(1)
    GG <- as.matrix(1)
    VV <- as.matrix(10*var(ts_data)) 
    WW <- as.matrix(w_v_ratio*10*var(ts_data))
    m0 <- as.matrix(mean(ts_data))              # prior mean for theta_0
    C0 <- as.matrix(var_prior)                        # prior variance for theta_0

    # wrap up all matrices and initial values
    matrices <- set_up_dlm_matrices(FF,GG,VV,WW)
    initial_states <- set_up_initial_states(m0, C0)

    # filtering
    results_filtered <- forward_filter(data, matrices, 
                                    initial_states)
    
    # smoothing
    results_smoothed <- backward_smoothing(data, matrices, 
                                        results_filtered)

    return(list(mean = results_smoothed$fnt, var=results_smoothed$Qnt))

}

#' @export
calc_signal_p_from_pvalues <- function(p_value, log_file_name=NULL){

    # get logit p values
    p_value <- as.numeric(unlist(p_value))
    filter_idx <- which(p_value > 0 & p_value<1.0)
    p_filt <- p_value[filter_idx]
    logit_pval <- qlogis(p_filt)
    logit_p_all <- rep(NA, length(p_value))
    logit_p_all[filter_idx] <- logit_pval

    # fit gaussian em mix 2 with logit p values
    fit_lp <- fit_gaussian_mix_em(logit_pval, KK=2, log_file_name)

    # get signal probability from fit result
    signal_cluster_idx <- which.min(fit_lp$mu)
    noise_cluster_idx <- which.max(fit_lp$mu)
    signal_w <- fit_lp$w[signal_cluster_idx]
    signal_p <- rep(0,length(p_value))
    signal_p[p_value==0] <- 1
    signal_p[p_value==1.0] <- 0 
    signal_mean <- fit_lp$mu[signal_cluster_idx]
    signal_std <- fit_lp$sigma[signal_cluster_idx]
    signal_p_filt <- pnorm(logit_pval, mean = signal_mean, sd = signal_std, lower.tail = FALSE)
    signal_p[filter_idx] <- signal_p_filt
    
    return( signal_p )
}

#' @export
get_signal_p_filtered_idx <- function(rowsum, rowsum_smoothed_mean, rowsum_smoothed_var, signal_p_threshold=0.5){

    p_value <- 1-pnorm(rowsum, rowsum_smoothed_mean, sqrt(rowsum_smoothed_var))
    signal_p <- calc_signal_p_from_pvalues(p_value)
    filtered_idx <- which(signal_p > signal_p_threshold)
    print('total rowsums : ')
    print(length(p_value))
    print('high signal rowsums : ')
    print(length(filtered_idx))

    return (list(filtered_idx=filtered_idx, p_value=p_value))
}


#' @export
fit_NDLM_and_EM_on_rowsum <- function(rowsum_vec, direction, params, chrom, resolution=1e4, output_path){

    # get parametrs 
    smoothing_factor <- params[[1]]
    var_prior <- params[[2]]
    rs_signal_p_threshold <- params[[3]]

    # Fit the NDLM to the row sums along the forward genomic coordinates
    rowsum_vec <- as.numeric(rowsum_vec)
    rowsum_smoothed_forward <- fit_NDLM_on_rowsum(rowsum_vec, chrom, smoothing_factor=smoothing_factor, var_prior=var_prior, resolution=resolution, output_path=output_path)

    # Fit the NDLM to the row sums along the reverse genomic coordinates
    rowsum_smoothed_reverse <- fit_NDLM_on_rowsum(rev(rowsum_vec), chrom, smoothing_factor=smoothing_factor, var_prior=var_prior, resolution=resolution, output_path=output_path)
    
    # get smoothing mean and variance by forward and reverse smoothing result 
    rowsum_smoothed_mean <- (rowsum_smoothed_forward$mean + rev(rowsum_smoothed_reverse$mean)) / 2
    rowsum_smoothed_var <- (rowsum_smoothed_forward$var + rev(rowsum_smoothed_reverse$var)) / 4

    # save NDLM smoothing plot
    ci_smoothed <- get_credible_interval(rowsum_smoothed_mean, rowsum_smoothed_var)
    save_ndlm_plot(rowsum_vec, rowsum_smoothed_mean, ci_smoothed, direction=direction, output_path=output_path, chrom, resolution)

    # Fit a Gaussian mixture model via EM from calculated p-values
    # and identify genomic row indices with high signal probabilities.
    result_em <- get_signal_p_filtered_idx(rowsum_vec, rowsum_smoothed_mean, rowsum_smoothed_var, signal_p_threshold=rs_signal_p_threshold )

    return(list(mod=list(mean=rowsum_smoothed_mean, var=rowsum_smoothed_var), high_signal_idx=result_em$filtered_idx, p_value=result_em$p_value))

}

#' @export
call_stripes_by_ndlm <- function(signal_p_path, chrom, resolution=1e4, max_distance=2e6){

    py_script_path <- system.file("python", "call_stripes.py", package = "Dynamix")
    params_path <- system.file("params.txt", package = "Dynamix")
    
    # source python functions 
    reticulate::source_python(py_script_path)
    
    # call stripes
    call_stripes(
        signal_p_path = as.character(signal_p_path),
        params_path = as.character(params_path),
        chrom = as.character(chrom),
        resolution    = as.integer(resolution),
        max_distance  = as.integer(max_distance)
    )

}

