#' Perform simulations.
#'
#' @param Y_mat_or an nxT matrix of outcome.
#' @param W_mat_or an nxT matrix of endogenous treatment.
#' @param Z n nx1 vector of an aggregate instrument.
#' @param share_t the share of the learning period.
#' @param share_rank the share of the rank of generalized fixed effects matrix in total rank.
#' @param rho_agg the share of Z_t during H_t generation.
#' @param rho_theta_w the exposure of W (endogenous treatment) to unobserved confounder. 
#' @param rho_theta_y the exposure of Y (outcome) to unobserved confounder. 
#' @param B number of replications for each design of simulation.
#' @param S if 'test' == TRUE, the number of simulations for standard error estimation.
#' @param test If true, compute the tests for the coverage rates.
#' @param K number of domain splits. It is used in density plot.
#' @param deg degrees of freedom for natural cubic splines. It is used in density plot.
#' @param height the height of the plot.
#' @param width the width of the plot.
#' @param folder_plot The folder where the plot is saved. NULL stands for the working directory.
#' @param save_sim If true, the function saves the simulation.
#' @param folder_sim If 'save_sim' is true, the folder where the simulation is saved. NULL stands for the working directory.
#' @param sim_name If 'save_sim' is true, the folder where the simulation is saved.
#' @param seed seed to set.
#'
#' @export



simulation <- function(Y_mat_or, W_mat_or, Z, share_t = 1/3, share_rank = 1/3, 
                       rho_agg = 0.5, rho_theta_w = 0.2, rho_theta_y = 0.3, B = 1000,
                       S = 300, test = FALSE, K = 300, deg = 4, height = 9*0.75, 
                       width = 16*0.75, folder_plot = NULL, save_sim = FALSE, 
                       folder_sim = NULL, sim_name = 'simulation_result',  seed = 1234){
  
  set.seed(seed)
  
  n <- dim(W_mat_or)[1]
  T <- dim(W_mat_or)[2]
  
  ### Parameters
  T_0 <- floor(share_t * T)
  rank <- floor(share_rank * T)
  pi_unit <- (1/T)*W_mat_or[,1:T]%*%(Z[1:T]-mean(Z[1:T])) / var(Z[1:T])
  Z_dem <- Z - mean(Z)
  
  
  ## Subtracting effects
  tau <- as.numeric(t((pi_unit-mean(pi_unit)))%*%Y_mat_or%*%Z_dem /
                      t((pi_unit-mean(pi_unit)))%*%W_mat_or%*%Z_dem)
  
  M_z <- diag(T) - Z_dem%*%solve(t(Z_dem)%*%Z_dem)%*%t(Z_dem)
  Y_mat <- Y_mat_or %*% M_z 
  W_mat <- W_mat_or %*% M_z
  
  
  ### Preparations
  results_Y <- data_preparation(Y_mat, rank)
  L_mat_Y <- results_Y[[1]]
  F_mat_Y <- results_Y[[2]]
  E_mat_Y <- results_Y[[3]]
  
  
  results_W <-data_preparation(W_mat, rank)
  L_mat_W <- results_W[[1]]
  F_mat_W <- results_W[[2]]
  E_mat_W <- results_W[[3]]
  
  Z_ts <- ts(Z, start = 1, end = T, frequency = 1) 
  Z_fit <- auto.arima(Z_ts)
  
  size_L_y <- norm(L_mat_Y,'f')^2 / (n*T)
  size_F_y <- norm(F_mat_Y,'f')^2 / (n*T)
  size_E_y <- norm(E_mat_Y,'f')^2 / (n*T)
  
  
  size_L_w <- norm(L_mat_W,'f')^2 / (n*T)
  size_F_w <- norm(F_mat_W,'f')^2 / (n*T)
  size_E_w <- norm(E_mat_W,'f')^2 / (n*T)
  
  rho_cross <- sum(diag(E_mat_Y%*%t(E_mat_W))) / (norm(E_mat_W,'f')*norm(E_mat_Y,'f'))
  
  scale_fs <- size_E_w 
  scale_rf <- size_E_y
  cov_mat_fs <- scale_fs*diag(1,T)
  cov_mat_rf <- scale_rf*diag(1,T)
  
  
  scale_str <- sqrt(size_L_w)  / (sd(Z)*sd(pi_unit)) 
  bias_scale_w <- 1
  bias_scale_y <- 3
  theta_w <- scale_str*bias_scale_w*(rho_theta_w*pi_unit  + sqrt(1-(rho_theta_w)^2)*rnorm(n)*sd(pi_unit))
  theta_y <- scale_str*bias_scale_y*(rho_theta_y*pi_unit  + sqrt(1-(rho_theta_y)^2)*rnorm(n)*sd(pi_unit)) 

  
  
  sigma_rf <- sqrt(size_E_y)
  sigma_fs <- sqrt(size_E_w)
  
  
  ### Simulations

  
  results_sim_4 <- do.call(rbind,lapply(1:B, function(b) { 
    estimate <- basic_sim(F_mat_W, L_mat_W ,F_mat_Y, L_mat_Y, cov_mat_fs, cov_mat_rf, scale_str*pi_unit, 
                          theta_w, theta_y, tau, rho_agg, rho_cross, T_0, Z_fit, test = test, S = S)
  })) 
  
  results_sim_3 <- do.call(rbind,lapply(1:B, function(b) { 
    estimate <- basic_sim(F_mat_W, 0, F_mat_Y, 0, cov_mat_fs, cov_mat_rf, scale_str*pi_unit,
                          theta_w, theta_y, tau, rho_agg, rho_cross, T_0, Z_fit, test = test, S = S)
  }))
  
  
  
  theta_w <- matrix(0,ncol = 1, nrow = n)
  theta_y <- theta_w
  
  results_sim_2 <- do.call(rbind,lapply(1:B, function(b) { 
    estimate <- basic_sim(F_mat_W, L_mat_W, F_mat_Y, L_mat_Y, cov_mat_fs, cov_mat_rf, scale_str*pi_unit, 
                          theta_w, theta_y, tau, rho_agg, rho_cross,T_0, Z_fit, test = test, S = S)
  }))
  
  results_sim_1 <- do.call(rbind,lapply(1:B, function(b) { 
    estimate <- basic_sim(F_mat_W, 0, F_mat_Y, 0, cov_mat_fs, cov_mat_rf, scale_str*pi_unit, 
                          theta_w, theta_y, tau, rho_agg, rho_cross, T_0, Z_fit, test = test, S = S)
  }))
  
  

  ### RMSE and bias calculations
  rmse_res_4 <- round(sqrt(colMeans(results_sim_4[,c(1:6)]^2)),3)
  bias_res_4 <- round(colMeans(results_sim_4[,c(1:6)]),3)
  table_res_4 <- cbind(rmse_res_4,bias_res_4)
  colnames(table_res_4) <- c('RMSE','Bias')
  
  
  rmse_res_3 <- round(sqrt(colMeans(results_sim_3[,c(1:6)]^2)),3)
  bias_res_3 <- round(colMeans(results_sim_3[,c(1:6)]),3)
  table_res_3 <- cbind(rmse_res_3,bias_res_3)
  colnames(table_res_3) <- c('RMSE','Bias')
  
  
  rmse_res_2 <- round(sqrt(colMeans(results_sim_2[,c(1:6)]^2)),3)
  bias_res_2 <- round(colMeans(results_sim_2[,c(1:6)]),3)
  table_res_2 <- cbind(rmse_res_2,bias_res_2)
  colnames(table_res_2) <- c('RMSE','Bias')
  
  
  rmse_res_1 <- round(sqrt(colMeans(results_sim_1[,c(1:6)]^2)),3)
  bias_res_1 <- round(colMeans(results_sim_1[,c(1:6)]),3)
  table_res_1 <- cbind(rmse_res_1, bias_res_1)
  colnames(table_res_1) <- c('RMSE','Bias')
  
  
  #### Table 1
  table_full_tau_0 <- cbind(table_res_1,table_res_2,table_res_3,table_res_4)
  print(xtable(table_full_tau_0, digits = 2))
  
  
  ### Coverage rates
  if (test == TRUE){
    
    row_rob <- c(mean(abs(results_sim_1[,7]) < qnorm(0.975)),
                 mean(abs(results_sim_2[,7]) < qnorm(0.975)),
                 mean(abs(results_sim_3[,7]) < qnorm(0.975)),
                 mean(abs(results_sim_4[,7]) < qnorm(0.975)))
    
    row_or <- c(mean(abs(results_sim_1[,8])< qnorm(0.975)),
                mean(abs(results_sim_2[,8]) < qnorm(0.975)),
                mean(abs(results_sim_3[,8]) < qnorm(0.975)),
                mean(abs(results_sim_4[,8]) < qnorm(0.975)))
    
    #### Table 2
    print(xtable(rbind(row_rob, row_or),digit = 2))
  }
  
  
  
  ### Density plot
  dens_our_des_4 <- my_density_function(results_sim_4[,5], K = K, deg = deg)
  dens_tsls_des_4 <- my_density_function(results_sim_4[,6], K = K,deg = deg)
  
  dens_our_des_2 <- my_density_function(results_sim_2[,5], K = K,deg = deg)
  dens_tsls_des_2 <- my_density_function(results_sim_2[,6], K = K,deg = deg)
  
  
  if (is.null(folder_plot)){
    pdf('fig_dens_full_orig.pdf', width = width, height = height)
  } else{
    dir <- paste(folder_plot, "fig_dens_full_orig.pdf", sep = "/")
    pdf(dir_2, width = width, height = height)
  }
  
  par(mfrow=c(1,2)) 
  
  plot(dens_our_des_2[,c(1,3)],lwd = 2, xlim = c(-1.5,1.5), type = 'l',lty = 1,xlab = 'estimate',
       ylab = 'Density',main = 'No unobserved shocks')
  lines(dens_tsls_des_2 [,c(1,3)],lwd = 2,lty = 2)
  abline(v = 0,lwd = 1,lty =2)
  legend('topright',lty = c(1,2),legend = c('Robust','TSLS'))
  
  plot(dens_our_des_4[,c(1,3)],lwd = 2, xlim = c(-1.5, 1.5), type = 'l',lty = 1,xlab = 'estimate',
       ylab = 'Density',main = 'Unobserved shocks')
  lines(dens_tsls_des_4 [,c(1,3)],lwd = 2,lty = 2)
  abline(v = 0,lwd = 1,lty =2)
  legend('topright',lty = c(1,2),legend = c('Robust','TSLS'))
  
  
  dev.off()
  
  
  if (save_sim) {
    
    if (is.null(folder_sim)){
      
    #sim_results <- list(results_sim_1, results_sim_2, results_sim_3, results_sim_4)
    save(results_sim_1, results_sim_2, results_sim_3, results_sim_4,
         file = 'sim_name')
    } else {
      dir <- paste(folder_sim, sim_name, sep = '/')
      save(results_sim_1, results_sim_2, results_sim_3, results_sim_4,
           file = dir)
    }
    
  }
  
  
}