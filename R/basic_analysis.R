#' Calculate the necessary.
#'
#' @param Y_mat an nxT matrix of outcome.
#' @param W_mat an nxT matrix of endogenous treatment.
#' @param Z an nx1 vector of an aggregate instrument.
#' @param psi aggregate covariates.
#' @param X_unit unit covariate.
#' @param D an nx1 vector of exposure.
#' @param T_0 the size of the learning period.
#' @param seed a seed to set.
#'
#' @return A list with weights, coefficients, standard errors, delta and pi and aggregated W, Y (observed and predicted).
#'
#' @export
#' 




basic_analysis <- function(Y_mat, W_mat, Z, X_unit, psi, D, T_0, seed = 1234) {
  
  set.seed(seed)
  
  
  n <- dim(Y_mat)[1]
  T <- dim(Y_mat)[2]
  
  ## basic analysis
  
  omega_or <-  weights_function(Y_mat,W_mat,Z, D_unit = D,
                                X_unit = X_unit, psi = psi, lambda = 100000)
  
  Z_dem <- lm(Z~psi)$residuals #substitute expectation 
  tau <- as.numeric(t(omega_or)%*%Y_mat%*%Z_dem/(t(omega_or)%*%W_mat%*%Z_dem))
  
  pi <- as.numeric(t(omega_or)%*%W_mat%*%Z_dem/(var_biased(Z_dem)*n*(T)))
  
  delta <- as.numeric(t(omega_or)%*%Y_mat%*%Z_dem/(var_biased(Z_dem)*n*(T)))
  
  ## our estimator
  
  Y_pre <- Y_mat[,1:T_0]
  W_pre <- W_mat[,1:T_0]
  Z_pre <- Z[1:T_0]
  psi_pre <- psi[1:T_0,]
  
  omega_rob <- weights_function(Y_pre, W_pre,Z_pre, D_unit = D,
                                X_unit = X_unit, psi = psi_pre)
  
  
  Y_post <- Y_mat[,(T_0+1):T]
  W_post <- W_mat[,(T_0+1):T]
  Z_post <- Z[(T_0+1):T]
  psi_post <- psi[(T_0+1):T,]
  
  Z_dem_pre <- lm(Z_pre~psi_pre)$residuals
  Z_dem_post <- lm(Z_post~psi_post)$residuals
  
  tau_rob <- as.numeric(t(omega_rob)%*%Y_post%*%Z_dem_post/(t(omega_rob)%*%W_post%*%Z_dem_post))
  
  pi_rob <- as.numeric(t(omega_rob)%*%W_post%*%Z_dem_post/(var_biased(Z_dem_post)*n*(T-T_0)))
  delta_rob <- as.numeric(t(omega_rob)%*%Y_post%*%Z_dem_post/(var_biased(Z_dem_post)*n*(T-T_0)))
  
  ## alt estimator
  
  tau_alt <- as.numeric(t(omega_or)%*%Y_post%*%Z_dem_post / (t(omega_or)%*%W_post%*%Z_dem_post))
  
  pi_alt <- as.numeric(t(omega_or)%*%W_post%*%Z_dem_post / (var_biased(Z_dem_post)*n*(T-T_0)))
  delta_alt <- as.numeric(t(omega_or)%*%Y_post%*%Z_dem_post / (var_biased(Z_dem_post)*n*(T-T_0)))
  
  
  
  ## variance estimation
  
  
  ## Cluster at yearly level
  
  #res_agg_rob <- (t(omega_rob)%*%Y_post - tau_rob*t(omega_rob)%*%W_post) / n
  #se_rob <- sd(res_agg_rob*Z_dem_post)/(abs(pi_rob)*var_biased(Z_dem_post)) / sqrt(T-T_0)
  
  #res_agg_alt <- (t(omega_or)%*%Y_post - tau*t(omega_or)%*%W_post)/n
  #se_alt <- sd(res_agg_alt*Z_dem_post)/(abs(pi_alt)*var_biased(Z_dem_post)) / sqrt(T-T_0)
  
  #res_agg <- (t(omega_or)%*%Y_mat - tau*t(omega_or)%*%W_mat) / n
  #se_or <- sd(res_agg*Z_dem) / (abs(pi)*var_biased(Z_dem)) / sqrt(T)
  
  ## computation based on simulations
  
  Z_fit <- auto.arima(Z_dem)
  S <- 200
  results_sd <- matrix(0, ncol = 3, nrow = S)
  
  for (j in 1:S) {
    art_Z_j <- as.numeric(simulate(Z_fit, T))
    art_Z_post_j <- art_Z_j[(T_0+1):T]
    art_Z_dem_j <- art_Z_j - mean(art_Z_j)                          
    art_Z_post_dem_j <- art_Z_post_j - mean(art_Z_post_j)
    rob_j <- mean(res_agg_rob*art_Z_post_dem_j) / (abs(pi_rob)*var_biased(art_Z_post_dem_j))
    or_alt_j <- mean(res_agg_alt*art_Z_post_dem_j) / (abs(pi_rob)*var_biased(art_Z_post_dem_j))
    or_j <-  mean(res_agg*art_Z_dem_j) / (abs(pi)*var_biased(art_Z_dem_j))
    results_sd[j,] <- c(rob_j,or_j,or_alt_j)	
  }
  
  
  se_rob <- sd(results_sd[,1])
  se_or <- sd(results_sd[,2])
  se_or_alt <- sd(results_sd[,3])
  
  
  ## inputs for graphs 
  
  rf_coefs_or <- Y_mat%*%Z_dem / var_biased(Z_dem)
  fs_coefs_or <- W_mat%*%Z_dem / var_biased(Z_dem)
  
  rf_coefs <- Y_post%*%Z_dem_post / var_biased(Z_dem_post)
  fs_coefs <- W_post%*%Z_dem_post / var_biased(Z_dem_post)
  
  ###
  W_agg <-  (t(omega_or) %*% W_mat) / n
  W_fit <- lm(t(W_agg)~Z+psi)$fitted.values
  
  Y_agg <- (t(omega_or) %*% Y_mat)/n
  Y_fit <- lm(t(Y_agg)~Z+psi)$fitted.values
  
  W_agg_post_rob <- t((t(omega_rob)%*%W_post)/n)
  W_fit_post_rob <- lm(W_agg_post_rob~Z_post+psi_post)$fitted.values
  
  ###Robust weigts
  Y_agg_post_rob <- t((t(omega_rob)%*%Y_post)/n)
  Y_fit_post_rob <- lm(Y_agg_post_rob~Z_post+psi_post)$fitted.values
  
  W_agg_pre_rob <- t((t(omega_rob)%*%W_pre)/n)
  W_fit_pre_rob <- lm(W_agg_pre_rob~Z_pre+psi_pre)$fitted.values
  
  Y_agg_pre_rob <- t((t(omega_rob)%*%Y_pre)/n)
  Y_fit_pre_rob <- lm(Y_agg_pre_rob~Z_pre+psi_pre)$fitted.values
  
  ## Reporting results
  
  weights_matrix <- cbind(omega_or,omega_rob)
  colnames(weights_matrix) <- c('original_weights','robust_weights')
  
  coefs_matrix <- rbind(c(pi,delta,tau),c(pi_rob,delta_rob,tau_rob),c(pi_alt,delta_alt, tau_alt))
  colnames(coefs_matrix) <- c('pi','delta','tau')
  rownames(coefs_matrix) <- c('original','robust','alternative')
  
  #se_vec <- c(se_or,se_rob,se_alt,se_or_alt,se_rob_alt,se_or_alt_alt)
  se_vec <- c(se_rob_alt,se_or_alt_alt)
  names(se_vec) <- c('se_or_cor','se_rob_cor','se_alt_cor')
  
  ind_coefs <- cbind(rf_coefs_or,fs_coefs_or,rf_coefs,fs_coefs)
  colnames(ind_coefs) <- c('rf_or','fs_or','rf_post','fs_post')
  
  agg_fits <- cbind(t(W_agg),W_fit,t(Y_agg),Y_fit,c(W_agg_pre_rob,W_agg_post_rob),
                    c(W_fit_pre_rob,W_fit_post_rob),c(Y_agg_pre_rob,Y_agg_post_rob),
                    c(Y_fit_pre_rob,Y_fit_post_rob))
  colnames(agg_fits) <- c('W_agg','W_fit','Y_agg','Y_fit','W_agg_rob','W_fit_rob',
                          'Y_agg_rob','Y_fit_rob')
  
  results <- list(weights_matrix,coefs_matrix,se_vec,ind_coefs,agg_fits)
  
  return(results)
  
}
