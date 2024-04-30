#' Present the simulation results from simulation file.
#'
#' @param sim_file String giving the simulation file path (including the file name).
#' @param plot_file String giving the plot file path (including the file name).
#' @param K number of domain splits. It is used in density plot.
#' @param deg degrees of freedom for natural cubic splines. It is used in density plot.
#' @param height the height of the plot.
#' @param width the width of the plot.
#'
#' @export



sim_from_file <- function(sim_file= 'simulation_result', plot_file = 'density_plot',
                          K = 300, deg = 4, height = 9*0.75, width = 16*0.75) {
  
  load(sim_file)
  
  
  ### RMSE and bias calculations
  rmse_res_4 <- round(sqrt(colMeans(results_sim_4[,c(1:6)]^2)), 3)
  bias_res_4 <- round(colMeans(results_sim_4[,c(1:6)]), 3)
  table_res_4 <- cbind(rmse_res_4,bias_res_4)
  colnames(table_res_4) <- c('RMSE','Bias')
  
  
  rmse_res_3 <- round(sqrt(colMeans(results_sim_3[,c(1:6)]^2)), 3)
  bias_res_3 <- round(colMeans(results_sim_3[,c(1:6)]), 3)
  table_res_3 <- cbind(rmse_res_3,bias_res_3)
  colnames(table_res_3) <- c('RMSE','Bias')
  
  
  rmse_res_2 <- round(sqrt(colMeans(results_sim_2[,c(1:6)]^2)), 3)
  bias_res_2 <- round(colMeans(results_sim_2[,c(1:6)]), 3)
  table_res_2 <- cbind(rmse_res_2,bias_res_2)
  colnames(table_res_2) <- c('RMSE','Bias')
  
  
  rmse_res_1 <- round(sqrt(colMeans(results_sim_1[,c(1:6)]^2)), 3)
  bias_res_1 <- round(colMeans(results_sim_1[,c(1:6)]), 3)
  table_res_1 <- cbind(rmse_res_1, bias_res_1)
  colnames(table_res_1) <- c('RMSE','Bias')
  
  #### Table 1
  table_full_tau_0 <- cbind(table_res_1,table_res_2,table_res_3,table_res_4)
  print(xtable(table_full_tau_0, digits = 2))
  
  
  ### Coverage rates
  if (dim(results_sim_1)[2] == 8){
    
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
  dens_our_des_4 <- density_function(results_sim_4[,5], K = K, deg = deg)
  dens_tsls_des_4 <- density_function(results_sim_4[,6], K = K,deg = deg)
  
  dens_our_des_2 <- density_function(results_sim_2[,5], K = K,deg = deg)
  dens_tsls_des_2 <- density_function(results_sim_2[,6], K = K,deg = deg)
  
  
  if (nchar(plot_file) >= 4) {
    if (substring(plot_file, nchar(plot_file) - 3) != ".pdf") {
      plot_file <- paste(plot_file, '.pdf', sep = '')
    }
  } else{
    plot_file <- paste(plot_file, '.pdf', sep = '')
  }
  
  pdf(plot_file, width = width, height = height)
  par(mfrow=c(1,2)) 
  
  plot(dens_our_des_2[,c(1,3)],lwd = 2, xlim = c(-1.5,1.5), type = 'l',lty = 1,xlab = 'estimate',
       ylab = 'Density',main = 'No unobserved shocks')
  lines(dens_tsls_des_2 [,c(1,3)],lwd = 2,lty = 2)
  abline(v = 0,lwd = 1,lty =2)
  legend('topright',lty = c(1,2),legend = c('Robust','TSLS'))
  
  plot(dens_our_des_4[,c(1,3)],lwd = 2, xlim = c(-1.5, 1.5), type = 'l',lty = 1,xlab = 'estimate',
       ylab = 'Density',main = 'Unobserved shocks')
  lines(dens_tsls_des_4 [,c(1,3)],lwd = 2,lty = 2)
  abline(v = 0,lwd = 1,lty = 2)
  legend('topright',lty = c(1,2),legend = c('Robust','TSLS'))
  
  dev.off()
  
}  
  