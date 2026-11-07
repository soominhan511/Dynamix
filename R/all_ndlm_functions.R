# Script : all_ndlm_functions.R
#   This script contains functions to perform Bayesian NDLM filtering and smoothing
#
# Reference : 
#   Prado, R., Ferreira, M. A. R. & West, M. (2021)
#   Time Series: Modeling, Computation, and Inference. (Chapman and Hall/CRC, Boca Raton).


#' @export
set_up_dlm_matrices <- function(FF, GG, VV, WW){
  return(list(FF=FF, GG=GG, VV=VV, WW=WW))
}

set_up_initial_states <- function(m0, C0){
  return(list(m0=m0, C0=C0))
}

#' @export
forward_filter <- function(data, matrices, initial_states){
  ## retrieve dataset
  yt <- data$yt
  T <- length(yt)
  
  ## retrieve a set of quadruples 
  # FF, GG, VV, WW are scalar
  FF <- matrices$FF  
  GG <- matrices$GG
  VV <- matrices$VV
  WW <- matrices$WW
  
  ## retrieve initial states
  m0 <- initial_states$m0
  C0 <- initial_states$C0
  
  ## create placeholder for results
  d <- dim(GG)[1]
  at <- matrix(NA, nrow=T, ncol=d)
  Rt <- array(NA, dim=c(d, d, T))
  ft <- numeric(T)
  Qt <- numeric(T)
  mt <- matrix(NA, nrow=T, ncol=d)
  Ct <- array(NA, dim=c(d, d, T))
  et <- numeric(T)
  
  
  for(i in 1:T){
    # moments of priors at t
    if(i == 1){
      at[i, ] <- GG %*% t(m0)
      Rt[, , i] <- GG %*% C0 %*% t(GG) + WW
      Rt[,,i] <- 0.5*Rt[,,i]+0.5*t(Rt[,,i]) 
    }else{
      at[i, ] <- GG %*% t(mt[i-1, , drop=FALSE])
      Rt[, , i] <- GG %*% Ct[, , i-1] %*% t(GG) + WW
      Rt[,,i] <- 0.5*Rt[,,i]+0.5*t(Rt[,,i]) 
    }
    
    # moments of one-step forecast:
    ft[i] <- t(FF) %*% (at[i, ]) 
    Qt[i] <- t(FF) %*% Rt[, , i] %*% FF + VV
    
    # moments of posterior at t:
    At <- Rt[, , i] %*% FF / Qt[i]
    et[i] <- yt[i] - ft[i]
    mt[i, ] <- at[i, ] + t(At) * et[i]
    Ct[, , i] <- Rt[, , i] - Qt[i] * At %*% t(At)
    Ct[,,i] <- 0.5*Ct[,,i] + 0.5*t(Ct[,,i]) 
  }
  cat("Forward filtering (NDLM) is completed!\n") # indicator of completion
  return(list(mt = mt, Ct = Ct, at = at, Rt = Rt, 
              ft = ft, Qt = Qt))
}

#' @export
forecast_function <- function(posterior_states, k, matrices){
  
  ## retrieve matrices
  FF <- matrices$FF
  GG <- matrices$GG
  WW <- matrices$WW
  VV <- matrices$VV
  mt <- posterior_states$mt
  Ct <- posterior_states$Ct
  
  ## set up matrices
  T <- dim(mt)[1] # time points
  d <- dim(mt)[2] # dimension of state parameter vector
  
  ## placeholder for results
  at <- matrix(NA, nrow = k, ncol = d)
  Rt <- array(NA, dim=c(d, d, k))
  ft <- numeric(k)
  Qt <- numeric(k)
  
  
  for(i in 1:k){
    ## moments of state distribution
    if(i == 1){
      at[i, ] <- GG %*% t(mt[T, , drop=FALSE])
      Rt[, , i] <- GG %*% Ct[, , T] %*% t(GG) + WW
      Rt[,,i] <- 0.5*Rt[,,i]+0.5*t(Rt[,,i]) 
    }else{
      at[i, ] <- GG %*% t(at[i-1, , drop=FALSE])
      Rt[, , i] <- GG %*% Rt[, , i-1] %*% t(GG) + WW
      Rt[,,i] <- 0.5*Rt[,,i]+0.5*t(Rt[,,i]) 
    }
    
    ## moments of forecast distribution
    ft[i] <- t(FF) %*% t(at[i, , drop=FALSE])
    Qt[i] <- t(FF) %*% Rt[, , i] %*% FF + VV
  }
  cat("Forecasting is completed!\n") # indicator of completion
  return(list(at=at, Rt=Rt, ft=ft, Qt=Qt))
}

#' @export
get_credible_interval <- function(mu, sigma2, 
                          quantile = c(0.025, 0.975)){
  z_quantile <- qnorm(quantile)
  bound <- matrix(0, nrow=length(mu), ncol=2)
  bound[, 1] <- mu + z_quantile[1]*sqrt(as.numeric(sigma2)) # lower bound
  bound[, 2] <- mu + z_quantile[2]*sqrt(as.numeric(sigma2)) # upper bound
  return(bound)
}

#' @export
backward_smoothing <- function(data, matrices, 
                               posterior_states){
  ## retrieve data 
  yt <- data$yt
  T <- length(yt) 
  
  ## retrieve matrices
  FF <- matrices$FF
  GG <- matrices$GG
  VV <- matrices$VV
  
  ## retrieve matrices
  mt <- posterior_states$mt
  Ct <- posterior_states$Ct
  at <- posterior_states$at
  Rt <- posterior_states$Rt
  
  ## create placeholder for posterior moments 
  mnt <- matrix(NA, nrow = dim(mt)[1], ncol = dim(mt)[2])
  Cnt <- array(NA, dim = dim(Ct))
  fnt <- numeric(T)
  Qnt <- numeric(T)
  for(i in T:1){
    # moments for the distributions of the state vector given D_T
    if(i == T){
      mnt[i, ] <- mt[i, ]
      Cnt[, , i] <- Ct[, , i]
      Cnt[, , i] <- 0.5*Cnt[, , i] + 0.5*t(Cnt[, , i]) 
    }else{
      inv_Rtp1<-solve(Rt[,,i+1])
      Bt <- Ct[, , i] %*% t(GG) %*% inv_Rtp1
      mnt[i, ] <- mt[i, ] + Bt %*% (mnt[i+1, ] - at[i+1, ])
      Cnt[, , i] <- Ct[, , i] + Bt %*% (Cnt[, , i + 1] - Rt[, , i+1]) %*% t(Bt)
      Cnt[,,i] <- 0.5*Cnt[,,i] + 0.5*t(Cnt[,,i]) 
    }
    # moments for the smoothed distribution of the mean response of the series
    fnt[i] <- t(FF) %*% t(mnt[i, , drop=FALSE])
    Qnt[i] <- t(FF) %*% t(Cnt[, , i]) %*% FF + VV
  }
  cat("Backward smoothing (NDLM) is completed!\n")
  return(list(mnt = mnt, Cnt = Cnt, fnt=fnt, Qnt=Qnt))
}

