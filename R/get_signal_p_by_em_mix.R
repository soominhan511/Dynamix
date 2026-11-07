# Script : get_signal_p_by_em_mix.R
#   This script contains functions to calculate signal probabilities using a ZINB+NB mixture EM algorithm

#' @importFrom dplyr %>% mutate filter pull select case_when left_join
#' @importFrom tidyr separate crossing
#' @importFrom rlang sym
#' @importFrom data.table fwrite
#' @import parallel
#' @import MASS
#' @import splines

# Run ZINB + NB mix EM to get signal probability for each interactions 

#' @export
split_indices_evenly <- function(min_idx, max_idx, group_num) {
 
  full_vector <- min_idx:max_idx
  breakpoints <- round(seq(0, length(full_vector), length.out = group_num + 1))
  result_groups <- list()
  for (i in 1:group_num) {
    start_pos <- breakpoints[i] + 1
    end_pos   <- breakpoints[i + 1]
    result_groups[[i]] <- full_vector[start_pos:end_pos]
  }
  return(result_groups)

}

#' @export
calc_zinb_pval <- function(x, mu, theta, w_zero){
    p_value = rep(1.0, length(x))
    non_zero_idx = (x!=0)
    if (length(theta)==1){ theta <- rep(theta, length(x)) }
    if (length(w_zero)==1){ w_zero <- rep(w_zero, length(x)) }
    p_value[non_zero_idx] <- (1-w_zero[non_zero_idx])*(1-pnbinom((x[non_zero_idx]-1), mu=mu[non_zero_idx], size=theta[non_zero_idx]))
    return(p_value)
}

#' @export
glm_nb_stable <- function(dat, fit_weight, df_D=6, boundary=NULL){
    
    nb_fit <- tryCatch({
        if (is.null(boundary)){
            glm.nb(count ~ splines::ns(log(D+2), df = df_D) + gc + map, data = dat, weights = fit_weight)
        }else{
            glm.nb(count ~ splines::ns(log(D+2), df = df_D, Boundary.knots = boundary) + gc + map, data = dat, weights = fit_weight)
        }
    }, error = function(e) {
        NULL 
    })
    return(nb_fit)
}

#' @export
merge_sparse_grp <- function(grp_l, dat, resolution, threshold=0.75){

    grp_d_merged <- list()
    for (i in 1:length(grp_l)){
        grp_i <- grp_l[[i]]
        d_l <- grp_i*resolution
        sample_idx <- which(dat$D %in% d_l)
        low_count_ratio <- sum(dat$count[sample_idx]<10)/length(sample_idx)
        if (low_count_ratio > threshold){
            remaining_groups <- grp_l[i:length(grp_l)]
            grp_d_merged[[i]] <- unlist(remaining_groups)
            return(grp_d_merged)
        }else{
            grp_d_merged[[i]] <- grp_i
        }
    }
    return (grp_d_merged)
}

#' @export
save_convergence_plot <- function(Q_result, plot_file_name, str_main){
    
    png(filename = plot_file_name, width = 800, height = 600, res = 100)
    plot(Q_result, 
        type = "b",              
        pch = 20,                  
        col = "darkgreen", 
        main = str_main, 
        xlab = "EM Iteration (s)", 
        ylab = "Q-Function Value (QQn)")
    dev.off()

}

#' @export
fit_gaussian_mix_em <- function(x, KK, log_file_name=NULL) {

    n= length(x)
    
    ## Run the EM algorithm
    ## Initialize the parameters
    w   = rep(1,KK)/KK 
    mu    = rnorm(KK, mean(x), sd(x))   #Random cluster centers randomly spread over the support of the data
    sigma    = rep(sqrt(var(x)/KK), KK)
    
    s  = 0
    sw = FALSE
    QQ = -Inf
    QQ.out = NULL
    epsilon = 10^(-6)


    ##Checking convergence of the algorithm
    while(!sw){
        
        ## E step
        v = array(0, dim=c(n,KK))
        for (k in 1:KK){
            v[,k] = log(w[k]) + dnorm(x, mu[k], sigma[k], log=TRUE)    #Compute the log of the weights
        }
        max_v <- do.call(pmax, as.data.frame(v[, 1:KK]))
        exp_v_matrix <- exp(v[, 1:KK] - max_v)
        sum_exp_v <- rowSums(exp_v_matrix)
        v <- sweep(exp_v_matrix, 1, sum_exp_v, "/")
        
        ## M step
        # Weights
        w <- colMeans(v)
        mu = rep(0, KK)
        for(k in 1:KK){
            mu[k]    = sum(v[,k]*x)/sum(v[,k])
        }
        # Standard deviations
        sigma =rep(0,KK)
        for(k in 1:KK){
            sigma[k] = sum( v[,k]*(x - mu[k])^2 )
            sigma[k] = sqrt(sigma[k]/sum(v[,k]))
        }
        
        
        ##Check convergence
        QQn = 0
        for (k in 1:KK){
            QQn = QQn + sum( v[,k]*(log(w[k]) + dnorm(x, mu[k], sigma[k], log=TRUE)) )
        }
        
        #print(QQn)
        if(abs(QQn-QQ)/abs(QQn)<epsilon){
            sw=TRUE
        }
        QQ = QQn
        QQ.out = c(QQ.out, QQ)
        s = s + 1
        log_msg <- paste(c(s, QQn, round(c(mu, sigma, w), 4)), collapse = " ")
        if (!is.null(log_file_name)){ write(log_msg, file = log_file_name, append = TRUE) }
        
        }
    
    # save qq plot 
    if (!is.null(log_file_name)){
        plot_file_name <- file.path(dirname(log_file_name), "EM_Qplot_gauss_mix.png")
        str_main <- "Gauss-Mix-EM Q value"
        save_convergence_plot(QQ.out, plot_file_name, str_main)
    }
    
    return(list(mu=mu, sigma=sigma, v=v, w=w))
    

}

#' @export
down_sample_dat <- function(x_dat, resolution = 1e4, down_ratio = 0.01, seed = 1010) {
  
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    x_dat$D_idx <- floor(x_dat$D / resolution)
    d_idx_unique <- unique(x_dat$D_idx)
    sampled_indices_list <- list()
    
    for (d_val in d_idx_unique) {
        sub_idx <- which(x_dat$D_idx == d_val)
        if (length(sub_idx) == 0) next
        
        sub_df <- x_dat[sub_idx, ]
        zero_idx <- which(sub_df$count == 0)
        non_zero_idx <- which(sub_df$count != 0)
        
        # apply down_ratio
        n_zero_sample <- floor(length(zero_idx) * down_ratio)
        n_non_zero_sample <- floor(length(non_zero_idx) * down_ratio)
        
        if (n_zero_sample == 0 && n_non_zero_sample == 0) {
        next
        }

        selected_idx <- c()
        if (n_zero_sample > 0) {
        sampled_zero <- sample(zero_idx, n_zero_sample, replace = FALSE)
        selected_idx <- c(selected_idx, sub_idx[sampled_zero])
        }
        if (n_non_zero_sample > 0) {
        sampled_non_zero <- sample(non_zero_idx, n_non_zero_sample, replace = FALSE)
        selected_idx <- c(selected_idx, sub_idx[sampled_non_zero])
        }
        
        if (length(selected_idx) > 0) {
        sampled_indices_list[[length(sampled_indices_list) + 1]] <- selected_idx
        }
    }
    
    sampled_indices <- unlist(sampled_indices_list)

    x_sampled <- x_dat[sampled_indices, ]
    
    return(x_sampled)
}

#' @export
fit_zinb_nb_mix_em <- function(x_dat, chr, down_ratio=0.01, KK=2, max_iter=350, df_D=6, file_path) {

    
    s       <- 1        # iteration number
    sw      <- FALSE        # EM convergence flag 
    QQ      <- -Inf
    QQ.out  <- NULL
    epsilon <- 10^(-5)      # EM convergence threshold

    # data 
    x_dat_all <- x_dat 
    dat_idx   <- paste0(x_dat_all$start1, "_", x_dat_all$start2)
    boundary <- range(log(x_dat_all$D+2))
    X_matrix_all <- model.matrix(~ splines::ns(log(D+2), df = df_D, Boundary.knots=boundary) + gc + map, data = x_dat_all)
    
    # down sampled data
    x_dat <- down_sample_dat(x_dat, down_ratio=down_ratio)
    x <- x_dat$count
    zero_idx <- which(x==0)
    non_zero_idx <- which(x!=0)
    n= length(x)
    
    # set log file
    current_D <- sort(unique(x_dat$D))
    str_D <- paste0(min(current_D),'-',max(current_D))
    logs_path <- file.path(file_path, "logs")
    if (!dir.exists(logs_path)) {
        dir.create(logs_path, showWarnings = FALSE, recursive = TRUE, mode = "0777")
    }
    log_file_name <- file.path(logs_path, paste0("log_D_", str_D, ".txt"))


    # monitoring
    zero_ratio <- sum(x==0)/length(x)
    log_msg_1 <- paste(c("zero ratio : ", zero_ratio))
    log_msg_2 <- capture.output(print(table(x)))
    log_msg_3 <- capture.output(length(x))
    write(log_msg_1, file = log_file_name, append = TRUE)
    write(log_msg_2, file = log_file_name, append = TRUE)
    write(log_msg_3, file = log_file_name, append = TRUE)
    
    # intialize weights of ZINB(zero-inflated negative binomial) + NB(negative binomial) mixture group
    # w : (ZINB weight, NB weight)
    # w_zero : structural zero weight inside ZINB group 
    w <- rep(1/KK, KK)
    w_zero <- 0.5 

    # initialize v : likelihood ratio of ZINB and NB 
    v = array(1/KK, dim=c(n,KK))
    v[zero_idx,1:KK] <- 0
    v[zero_idx,1] <- 1        # zero count only in ZINB
    
    
    # initialize v_zin : likelihood ratio of structural zero and non-structural zero inside ZINB
    v_zin = array(0, dim=c(n,2))
    v_zin[zero_idx,1] <- w_zero
    v_zin[zero_idx,2] <- 1-w_zero
    v_zin[non_zero_idx,1] <- 0      # non zero count only belongs to non-structural zero group
    v_zin[non_zero_idx,2] <- 1
    
    
    # build design matrix
    # predict counts from linear function of three coefficients
    # coefficient1 : natural spline of log( distance+2 ), 
    # coefficient2 : normalized log( gc ratio )
    # coefficient3 : normalized log( map ratio )
    X_matrix <- model.matrix(~ splines::ns(log(D+2), df = df_D, Boundary.knots=boundary) + gc + map, data = x_dat)
    
    while(!sw){   

        ## E step
        q = array(0, dim=c(n,KK))
        nb_fit <- list()
        for (k in 1:KK){
            
            if (k==1){      # ZINB group 

                # nb fit weights
                fit_w <- v[,k]*v_zin[,2]

                # nb fit 
                nb_fit_k <- glm_nb_stable(x_dat,fit_weight=fit_w, df_D=df_D, boundary=boundary)
                if (is.null(nb_fit_k)) { return(NULL) }
                mu_full <- exp(X_matrix %*% coef(nb_fit_k))

                # calculate zinb likelihood 
                log_p_noise <- log(1-w_zero) + dnbinom(x, mu=mu_full, size=nb_fit_k$theta, log=TRUE)
                q[zero_idx,k] = log(w[k]) + log(w_zero + exp(log_p_noise[zero_idx]))
                q[non_zero_idx,k] = log(w[k]) + log_p_noise[non_zero_idx]

            }else{      # NB group

                # nb fit weights 
                fit_w <- v[,k]

                # nb fit
                nb_fit_k <- glm_nb_stable(x_dat,fit_weight=fit_w, df_D=df_D, boundary=boundary)
                if (is.null(nb_fit_k)) { return(NULL) }
                mu_full <- exp(X_matrix %*% coef(nb_fit_k))

                # calculate nb likelihood
                q[,k] = log(w[k]) + dnbinom(x, mu=mu_full, size=nb_fit_k$theta, log=TRUE)

            }
            
            nb_fit[[k]] <- nb_fit_k
            
        }
        
        ##Check convergence
        if (s > 1){
            QQn <- 0 
            if (KK==1){
                QQn_zero <- sum(v_zin[zero_idx,1]*log(w_zero)) + sum(v_zin[zero_idx,2]*log_p_noise[zero_idx])
                QQn_non_zero <- sum(v_zin[non_zero_idx,2]*log_p_noise[non_zero_idx])
                QQn = QQn_zero + QQn_non_zero
            }else{
                for (k in 1:KK){
                    QQn = QQn + sum( v[,k]*q[,k] ) 
                }
            }
            
            
            if (s > 2) { 
                # check if QQ is increasing and iterations are sufficient to prevent convergence to a suboptimal state
                if ((abs(QQn - QQ) / abs(QQ) < epsilon) && ((QQn > QQ && s >= 130) || (QQn < QQ && s >= 200))) {      
                    sw <- TRUE   
                    QQ.out <- c(QQ.out, QQ)   
                    log_msg_converged <- paste(c("EM converged", "| w:", round(w, 4), "| w_zero:", round(w_zero, 4)), collapse = " ")
                    print(paste0(str_D, ' -> ', log_msg_converged))
                    write(log_msg_converged, file = log_file_name, append = TRUE)
                    break
                } else if (s == max_iter) {
                    sw <- TRUE
                    QQ.out <- c(QQ.out, QQ)   
                    log_msg_max_reached <- paste(c("Max iter reached", "| w:", round(w, 4), "| w_zero:", round(w_zero, 4)), collapse = " ")
                    print(paste0(str_D, ' -> ', log_msg_max_reached))
                    write(log_msg_max_reached, file = log_file_name, append = TRUE)
                    break
                }            
            }
            QQ     <- QQn   
            QQ.out <- c(QQ.out, QQ)   
            
        }
        
        # update v to likelihood ratio of ZINB and NB with numerical stability
        v = array(0, dim=c(n,KK))
        v <- q 
        max_v <- do.call(pmax, as.data.frame(v[, 1:KK]))
        exp_v_matrix <- exp(v[, 1:KK] - max_v)
        sum_exp_v <- rowSums(as.matrix(exp_v_matrix))
        v <- sweep(as.matrix(exp_v_matrix), 1, sum_exp_v, "/")
        v[zero_idx, 1:KK] <- 0
        v[zero_idx, 1] <- 1         # for zero count, likelihood ratio of ZINB to NB should be 1:0

        # update v_zin to likelihood ratio of structural-zero and non-structural-zero with numerical stability
        v_zin <- cbind(log(w_zero), log_p_noise)
        max_tmp <- do.call(pmax, as.data.frame(v_zin))
        exp_diff_tmp <- exp(v_zin - max_tmp)
        sum_tmp <- rowSums(as.matrix(exp_diff_tmp))
        v_zin <- sweep(as.matrix(exp_diff_tmp), 1, sum_tmp, "/")
        v_zin[non_zero_idx,1] <- 0 
        v_zin[non_zero_idx,2] <- 1         # for non-zero count, likelihood ratio of structural-zero to non-structural-zero should be 0:1
        

        # update w, w_zero
        w <- colMeans(v)
        w_zero <- sum(v_zin[,1]*v[,1])/sum(v[,1])
        
        # write log for monitoring
        qq_print <- if(s == 1) "Initial" else round(QQn, 2)
        log_msg <- paste(c("Iteration:", s, "| Q:", qq_print, "| w:", round(w, 3), "| w_zero:", round(w_zero, 3)), collapse = " ")
        #tbl <- table(x_dat$count[v[, 2] > 1/KK])
        #captured_tbl <- capture.output(print(tbl))
        
        #write(captured_tbl, file = log_file_name, append = TRUE)
        write(log_msg, file = log_file_name, append = TRUE)
        
        s <- s + 1
        
        
        }

    
    # save QQ plot
    plot_file_name <- file.path(logs_path, paste0("EM_Qplot_D_", str_D, ".png"))
    str_main <- paste0("ZINB-NB-EM Q value (Distance: ", str_D, ")")
    save_convergence_plot(QQ.out, plot_file_name, str_main)
    
    # write log for converged result
    for (k in 1:KK){
        write(paste0("\n--- Component [", k, "] Summary ---"), file = log_file_name, append = TRUE)
        model_summary_txt <- capture.output(print(summary(nb_fit[[k]])))
        write(model_summary_txt, file = log_file_name, append = TRUE)
    }    

    # calculate p-value with ZINB
    zinb_mod <- nb_fit[[1]]
    theta_zinb <- zinb_mod$theta
    mu_zinb <- exp(X_matrix_all %*% coef(zinb_mod))
    p_value <- calc_zinb_pval(x_dat_all$count, mu_zinb, theta_zinb, w_zero)

    return(list(mod=nb_fit, w=w, w_zero = w_zero, p_value=p_value, mu=mu_zinb, names=dat_idx, D_grp=str_D))
}

#' @export
get_nb_model_info <- function(mod){
    
    mod_summary <- summary(mod)
    coef_table <- mod_summary$coefficients 
    
    # 1. get coefficients
    theta <- mod$theta
    all_coefs <- coef(mod)
    ns_coef_names <- grep("ns\\(", names(all_coefs), value = TRUE)
    
    coef_df <- data.frame(
        intercept = if("(Intercept)" %in% names(all_coefs)) all_coefs["(Intercept)"] else NA_real_,
        ns1 = if(length(ns_coef_names) >= 1) all_coefs[ns_coef_names[1]] else NA_real_,
        ns2 = if(length(ns_coef_names) >= 2) all_coefs[ns_coef_names[2]] else NA_real_,
        ns3 = if(length(ns_coef_names) >= 3) all_coefs[ns_coef_names[3]] else NA_real_,
        ns4 = if(length(ns_coef_names) >= 4) all_coefs[ns_coef_names[4]] else NA_real_,
        ns5 = if(length(ns_coef_names) >= 5) all_coefs[ns_coef_names[5]] else NA_real_,
        ns6 = if(length(ns_coef_names) >= 6) all_coefs[ns_coef_names[6]] else NA_real_,
        gc  = if("gc" %in% names(all_coefs)) all_coefs["gc"] else NA_real_,
        map = if("map" %in% names(all_coefs)) all_coefs["map"] else NA_real_,
        theta = theta, 
        stringsAsFactors = FALSE
    )
    
    # 2. get p_values of coefficients
    p_vals <- if (!is.null(coef_table) && ncol(coef_table) >= 4) coef_table[, 4] else rep(NA_real_, nrow(coef_table))
    if (is.null(names(p_vals))) names(p_vals) <- rownames(coef_table)
    
    ns_pval_names <- grep("ns\\(", names(p_vals), value = TRUE)
    
    pval_df <- data.frame(
        intercept = if("(Intercept)" %in% names(p_vals)) p_vals["(Intercept)"] else NA_real_,
        ns1 = if(length(ns_pval_names) >= 1) p_vals[ns_pval_names[1]] else NA_real_,
        ns2 = if(length(ns_pval_names) >= 2) p_vals[ns_pval_names[2]] else NA_real_,
        ns3 = if(length(ns_pval_names) >= 3) p_vals[ns_pval_names[3]] else NA_real_,
        ns4 = if(length(ns_pval_names) >= 4) p_vals[ns_pval_names[4]] else NA_real_,
        ns5 = if(length(ns_pval_names) >= 5) p_vals[ns_pval_names[5]] else NA_real_,
        ns6 = if(length(ns_pval_names) >= 6) p_vals[ns_pval_names[6]] else NA_real_,
        gc  = if("gc" %in% names(p_vals)) p_vals["gc"] else NA_real_,
        map = if("map" %in% names(p_vals)) p_vals["map"] else NA_real_,
        stringsAsFactors = FALSE
    )

    return(list(coef = coef_df, pval = pval_df))
}

#' @export
save_result_model_gmix <- function(mod, file_path){
    
    signal_idx <- which.min(mod$mu)
    noise_idx <- which.max(mod$mu)

    # mean of logit(p_value) gauss-mix
    signal_mean <- mod$mu[signal_idx]
    noise_mean <- mod$mu[noise_idx]

    # standard deviation of logit(p_value) gauss-mix
    signal_std <- mod$sigma[signal_idx]
    noise_std <- mod$sigma[noise_idx]

    # weight of logit(p_value) gauss-mix
    noise_weight <- mod$w[noise_idx]
    signal_weight <- mod$w[signal_idx] 

    # data frame 
    row_df <- data.frame(
        signal_mean = signal_mean,
        noise_mean = noise_mean,
        signal_std = signal_std,
        noise_std = noise_std,
        noise_weight = noise_weight,
        signal_weight = signal_weight,
        stringsAsFactors = FALSE
    )
    rownames(row_df) <- "logit_pvalue_guass_mix"

    # save result to csv 
    output_file_name <- file.path(file_path, "logit_pvalue_gauss_mix_EM_coefficients.csv")
    data.table::fwrite(row_df, output_file_name, row.names = TRUE, quote = FALSE, sep = ",")
    cat("Result of gaussian mix EM of logit(p_value) saved in : ", output_file_name, "\n")

    return(row_df)

}

#' @export
save_result_features <- function(feature_df, nb_mix_result, signal_p_result, file_path){
    
    feature_df$feature_names <- paste0(feature_df$start1, '_', feature_df$start2)
    chr_name <- feature_df$chr[1]

    all_names <- unlist(lapply(nb_mix_result, function(res) res$names))
    all_p_val <- unlist(lapply(nb_mix_result, function(res) res$p_value))
    all_mu    <- unlist(lapply(nb_mix_result, function(res) res$mu))
    
    metrics_df <- data.frame(
        feature_names = as.character(all_names),
        p_value = as.numeric(all_p_val),
        mu = as.numeric(all_mu),
        stringsAsFactors = FALSE
    )
    
    signal_df <- data.frame(
        feature_names = as.character(signal_p_result$names),
        signal_p = as.numeric(signal_p_result$signal_p),
        stringsAsFactors = FALSE
    )
    
    result_df <- feature_df %>%
        left_join(metrics_df, by = "feature_names") %>%
        left_join(signal_df, by = "feature_names") %>%
        select(-feature_names)
    
    # save csv 
    output_file_name <- file.path(file_path, paste0(chr_name,"_pair_features_stats.csv"))
    data.table::fwrite(result_df, output_file_name, row.names = FALSE, quote = FALSE, sep = ",")
    cat("signal probability and p value saved in :", output_file_name, "\n")
    feature_file_name <- file.path(dirname(file_path), paste0(chr_name, "_pair_features.csv"))
    file.remove(feature_file_name)
    
    return(result_df)
}


#' @export
save_result_model <- function(nb_mix_result, file_path){
    
    row_list <- list()
    
    for (res in nb_mix_result) {
        D_grp <- res$D_grp
        zinb_mod <- res$mod[[1]]
        nb_mod <- res$mod[[2]]
        
        noise_weight <- res$w[1]
        signal_weight <- res$w[2]
        zinb_struct_zero_weight <- res$w_zero
        
        noise_info <- get_nb_model_info(zinb_mod)
        signal_info <- get_nb_model_info(nb_mod)
        
        noise_coef <- as.data.frame(noise_info$coef)
        if (ncol(noise_coef) > 0) colnames(noise_coef) <- paste0("noise_", colnames(noise_coef))
        
        noise_pval <- as.data.frame(noise_info$pval)
        if (ncol(noise_pval) > 0) colnames(noise_pval) <- paste0("p_value_noise_", colnames(noise_pval))
        
        signal_coef <- as.data.frame(signal_info$coef)
        if (ncol(signal_coef) > 0) colnames(signal_coef) <- paste0("signal_", colnames(signal_coef))
        
        signal_pval <- as.data.frame(signal_info$pval)
        if (ncol(signal_pval) > 0) colnames(signal_pval) <- paste0("p_value_signal_", colnames(signal_pval))
        
        weights_df <- data.frame(
            noise_weight = noise_weight,
            signal_weight = signal_weight,
            zinb_struct_zero_weight = zinb_struct_zero_weight,
            stringsAsFactors = FALSE
        )
        
        row_df <- cbind(
            data.frame(D_grp = D_grp, stringsAsFactors = FALSE),
            noise_coef, noise_pval, signal_coef, signal_pval, weights_df
        )
        
        row_list[[length(row_list) + 1]] <- row_df
    }
    
    result_df <- dplyr::bind_rows(row_list)
    if (!dir.exists(file_path)) {
        dir.create(file_path, showWarnings = FALSE, recursive = TRUE, mode = "0777")
    }
    
    # save model result to csv 
    output_file_name <- file.path(file_path, "ZINB_NB_mix_EM_coefficients.csv")
    data.table::fwrite(result_df, output_file_name, row.names = FALSE, quote = FALSE, sep = ",")
    cat("ZINB NB mix EM coefficients saved in : ", output_file_name, "\n")
    
    return(result_df)
}




#' @export
calc_signal_p_from_pvalues_nb_mix <- function(nb_mix_result, log_file_name){

    names <- unlist(lapply(nb_mix_result, function(res) res$names))
    p_value <- unlist(lapply(nb_mix_result, function(res) res$p_value))

    # get logit p values
    filter_idx <- which(p_value > 0 & p_value<1.0)
    p_filt <- p_value[filter_idx]
    logit_pval <- qlogis(p_filt)
    logit_p_all <- rep(NA, length(p_value))
    logit_p_all[filter_idx] <- logit_pval

    # fit gaussian em mix 2 with logit p values
    fit_lp <- fit_gaussian_mix_em(logit_pval, KK=2, log_file_name)

    # save model result 
    file_path <- dirname(dirname(log_file_name))
    save_result_model_gmix(fit_lp, file_path)

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
    
    return(list(names=names, signal_p=signal_p))

}

#' @export
get_signal_p <- function(feature_df, file_path, n_cores=1, down_ratio=0.01, D_grp_size=200e4, resolution=1e4, max_distance=2e6, max_iter=350){
    
    # split by distance groups
    zinb_nb_mix_result<- list()
    min_D_idx <- 0 
    max_D_idx <- as.integer(max_distance/resolution)
    D_grp_num <- as.integer(D_grp_size/resolution)
    split_D_grp_num <- as.integer((max_D_idx-min_D_idx+1)/D_grp_num)
    grp_D_tmp <- split_indices_evenly(min_D_idx, max_D_idx, split_D_grp_num)
    grp_D <- grp_D_tmp

    # get down_sample_ratio
    params_path <- system.file("params.txt", package = "Dynamix")
    params_df <- read.table(params_path, header = FALSE, sep = "", stringsAsFactors = FALSE, comment.char = "#")
    down_ratio <- as.numeric(params_df[params_df$V1 == "down_sample_ratio", 2])
    
    # run parallel loop 
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl,library("MASS"))
    parallel::clusterEvalQ(cl,library("splines"))
    parallel::clusterExport(cl, varlist = c("feature_df", "resolution", "fit_zinb_nb_mix_em", "file_path", "glm_nb_stable", "grp_D", "down_ratio",
                             "calc_zinb_pval", "fit_gaussian_mix_em", "calc_signal_p_from_pvalues_nb_mix","max_iter", "max_distance", "save_convergence_plot", "save_result_features", "save_result_model", "down_sample_dat"), envir = environment())

    # run ZINB+NB mix EM parallel
    zinb_nb_mix_result <- parallel::parLapply(cl, 1:length(grp_D), function(d_idx) {

        # get samples in D group
        D_idx_grp <- grp_D[[d_idx]]
        D_grp <- resolution * D_idx_grp
        sampled_idx <- which(feature_df$D %in% D_grp)
        sampled_dat <- feature_df[sampled_idx, ]
        current_D <- sort(unique(D_grp))
        str_D <- paste0(min(current_D),'-',max(current_D))
        log_file_name <- file.path(file_path, paste0("error_log_", str_D, ".txt"))
        chr_name <- feature_df$chr[1]

        # run EM 
        result <- tryCatch({
            fit_zinb_nb_mix_em(sampled_dat, chr_name, down_ratio=down_ratio, KK=2, max_iter=max_iter, df_D=6, file_path)
        }, error = function(e) {
            error_time <- Sys.time()
            error_msg <- paste0(
                "\n===================================================\n",
                "[ERROR ALERT] ", error_time, "\n",
                "Distance (D): ", str_D, " (d_idx: ", str_D, ") -> error occured.\n",
                "error message : ", e$message, "\n",
                "===================================================\n"
            )
            write(error_msg, file = log_file_name, append = TRUE)
            return(NULL)
        })
        
        return(result)
    })
        
    if (any(sapply(zinb_nb_mix_result, function(x) any(is.na(x)))) || any(sapply(zinb_nb_mix_result, is.null))) {
        return(NULL)
    }

    # get signal_p with p-values by fitting guassian-mix2 EM 
    log_file_name_gauss_EM <- file.path(file_path, "logs", "logit_pvalue_gauss_mix_EM_logs.txt")
    signal_p_result <- calc_signal_p_from_pvalues_nb_mix(zinb_nb_mix_result, log_file_name_gauss_EM)

    
    # save signal_p and p_values
    save_result_features(feature_df, zinb_nb_mix_result, signal_p_result, file_path)
    
    # save model results
    save_result_model(zinb_nb_mix_result, file_path)


}
