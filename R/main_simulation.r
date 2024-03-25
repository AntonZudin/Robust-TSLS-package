`library(mvtnorm)
library('forecast')
library(xtable)
library(parallel)



rm(list = ls())
seed = 1234
set.seed(seed)


#setwd('/Users/darkhangelsky/Library/CloudStorage/Dropbox-Personal/Research/aggregate_iv/code/simulation_paper/Nakamura')
setwd("C:/Users/Serge/OneDrive - экономический факультет МГУ им. М.В.Ломоносова/Митя/Robust-TSLS")
#source('implementation/synth_weights.r')
source('package/density_function.R')
source('package/weights_function.R')
source('package/data preparation.R')
source('implementation/basic_sim.r')
source('implementation/test_function.r')



## Loading raw data

load("data/nak_out_data.RData")

Y_mat_or <- data_nak[[1]]#[-c(1,26,29),]
W_mat_or <- data_nak[[2]]#[-c(1,26,29),]
Z <- data_nak[[3]]



n <- dim(W_mat_or)[1]
T <- dim(W_mat_or)[2]

## Parameters 

share_t <- 1/3
share_rank <- 1/3
T_0 <- floor(T*share_t)
rank <- floor(T*share_rank)
rho_agg <- 0.5
rho_theta_w <- 0.2
rho_theta_y <- 0.3
pi_unit <- (1/T)*W_mat_or[,1:T]%*%(Z[1:T]-mean(Z[1:T]))/var(Z[1:T]) #calculate D_i
Z_dem <- Z - mean(Z)



tau <- as.numeric(t((pi_unit-mean(pi_unit)))%*%Y_mat_or%*%Z_dem/
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

Z_ts <- ts(Z, start=1, end=T, frequency=1) 
Z_fit <- auto.arima(Z_ts)

size_L_y <- norm(L_mat_Y,'f')^2 / (n*T)
size_F_y <- norm(F_mat_Y,'f')^2 / (n*T)
size_E_y <- norm(E_mat_Y,'f')^2 / (n*T)


size_L_w <- norm(L_mat_W,'f')^2 / (n*T)
size_F_w <- norm(F_mat_W,'f')^2 / (n*T)
size_E_w <- norm(E_mat_W,'f')^2 / (n*T)

rho_cross <- sum(diag(E_mat_Y%*%t(E_mat_W)))/(norm(E_mat_W,'f')*norm(E_mat_Y,'f'))

scale_fs <- size_E_w 
scale_rf <- size_E_y
cov_mat_fs <- scale_fs*diag(1,T)
cov_mat_rf <- scale_rf*diag(1,T)


scale_str <- sqrt(size_L_w)/(sd(Z)*sd(pi_unit)) 
bias_scale_w <- 1
bias_scale_y <- 3
theta_w <- scale_str*bias_scale_w*(rho_theta_w*pi_unit  + sqrt(1-(rho_theta_w)^2)*rnorm(n)*sd(pi_unit))
theta_y <- scale_str*bias_scale_y*(rho_theta_y*pi_unit  + sqrt(1-(rho_theta_y)^2)*rnorm(n)*sd(pi_unit)) 


sigma_rf <- sqrt(size_E_y)
sigma_fs <- sqrt(size_E_w)


### Simulation 1: results
B <- 1000

start_time <- Sys.time()

#simulation 5 for Table 2
n_sim_5 <- 100
T_sim_5 <- 80
T_0_sim_5 <- floor(T_sim_5/3)


results_sim_5 <- do.call(rbind,lapply(1:B, function(b){ 
  estimate <-test_function(n_sim_5, T_sim_5, T_0_sim_5, pi_unit, theta_w, theta_y, sigma_rf, sigma_fs,
                           tau, rho_agg, Z_fit, S = 200)
}))



results_sim_4 <- do.call(rbind,lapply(1:B, function(b){ 
		estimate <- basic_sim(F_mat_W, L_mat_W ,F_mat_Y, L_mat_Y, cov_mat_fs, cov_mat_rf, scale_str*pi_unit, 
		                      theta_w,theta_y, tau, rho_agg, rho_cross, T_0, Z_fit, test = TRUE)
})) 




rmse_res <- round(sqrt(colMeans(results_sim_4[,c(1:6)]^2)),3)
bias_res <- round(colMeans(results_sim_4[,c(1:6)]),3)

table_res_4 <- cbind(rmse_res,bias_res)
colnames(table_res_4) <- c('RMSE','Bias')



results_sim_3 <- do.call(rbind,lapply(1:B, function(b){ 
		estimate <- basic_sim(F_mat_W, 0, F_mat_Y, 0, cov_mat_fs, cov_mat_rf, scale_str*pi_unit,
		                      theta_w, theta_y, tau, rho_agg, rho_cross, T_0, Z_fit, test = TRUE)
}))


rmse_res <- round(sqrt(colMeans(results_sim_3[,c(1:6)]^2)),3)
bias_res <- round(colMeans(results_sim_3[,c(1:6)]),3)

table_res_3 <- cbind(rmse_res,bias_res)
colnames(table_res_3) <- c('RMSE','Bias')



theta_w <- matrix(0,ncol = 1, nrow = n)
theta_y <-  theta_w

results_sim_2 <- do.call(rbind,lapply(1:B, function(b){ 
		estimate <- basic_sim(F_mat_W, L_mat_W, F_mat_Y, L_mat_Y, cov_mat_fs, cov_mat_rf, scale_str*pi_unit, 
		                      theta_w, theta_y, tau, rho_agg, rho_cross,T_0, Z_fit, test = TRUE)
}))


rmse_res <- round(sqrt(colMeans(results_sim_2[,c(1:6)]^2)),3)
bias_res <- round(colMeans(results_sim_2[,c(1:6)]),3)

table_res_2 <- cbind(rmse_res,bias_res)
colnames(table_res_2) <- c('RMSE','Bias')



theta_w <- matrix(0,ncol = 1, nrow = n)
theta_y <-  theta_w

results_sim_1 <- do.call(rbind,lapply(1:B, function(b){ 
		estimate <- basic_sim(F_mat_W, 0, F_mat_Y, 0, cov_mat_fs, cov_mat_rf, scale_str*pi_unit, 
		                      theta_w,theta_y, tau, rho_agg, rho_cross, T_0, Z_fit, test = TRUE)
}))


rmse_res <- round(sqrt(colMeans(results_sim_1[,c(1:6)]^2)),3)
bias_res <- round(colMeans(results_sim_1[,c(1:6)]),3)

table_res_1 <- cbind(rmse_res,bias_res)
colnames(table_res_1) <- c('RMSE','Bias')





table_full_tau_0 <- cbind(table_res_1,table_res_2,table_res_3,table_res_4)

end_time <- Sys.time()
end_time - start_time

xtable(table_full_tau_0, digits = 2)


########### Tests


row_rob <- c(mean(abs(results_sim_1[,7]) < qnorm(0.975)),
mean(abs(results_sim_2[,7]) < qnorm(0.975)),
mean(abs(results_sim_3[,7]) < qnorm(0.975)),
mean(abs(results_sim_4[,7]) < qnorm(0.975)),
mean(abs(results_sim_5[,3]) < qnorm(0.975)))

row_or <- c(mean(abs(results_sim_1[,8])< qnorm(0.975)),
mean(abs(results_sim_2[,8]) < qnorm(0.975)),
mean(abs(results_sim_3[,8]) < qnorm(0.975)),
mean(abs(results_sim_4[,8]) < qnorm(0.975)),
mean(abs(results_sim_5[,4]) < qnorm(0.975)))

xtable(rbind(row_rob, row_or),digit = 2)

#mean(abs(results_sim_2[,7]) > qnorm(0.975))
#mean(abs(results_sim_3[,7]) > qnorm(0.975))
#mean(abs(results_sim_4[,7]) > qnorm(0.975))



######### Plots
dens_our_des_4 <- my_density_function(results_sim_4[,5], K = 300,deg = 4)
dens_tsls_des_4 <- my_density_function(results_sim_4[,6], K = 300,deg = 4)

dens_our_des_2 <- my_density_function(results_sim_2[,5], K = 300,deg = 4)
dens_tsls_des_2 <- my_density_function(results_sim_2[,6], K = 300,deg = 4)


pdf(file = "./plots/fig_dens_full_orig.pdf",   width = 16*0.75,
   	 height = 9*0.75) 
par(mfrow=c(1,2)) 

plot(dens_our_des_2[,c(1,3)],lwd = 2, xlim = c(-1.5,1.5), type = 'l',lty =1,xlab = 'estimate',
ylab = 'Density',main = 'No unobserved shocks')
lines(dens_tsls_des_2 [,c(1,3)],lwd = 2,lty =2)
abline(v = 0,lwd = 1,lty =2)
legend('topright',lty = c(1,2),legend = c('Robust','TSLS'))

plot(dens_our_des_4[,c(1,3)],lwd = 2, xlim = c(-1.5,1.5),type = 'l',lty =1,xlab = 'estimate',
ylab = 'Density',main = 'Unobserved shocks')
lines(dens_tsls_des_4 [,c(1,3)],lwd = 2,lty =2)
abline(v = 0,lwd = 1,lty =2)
legend('topright',lty = c(1,2),legend = c('Robust','TSLS'))


dev.off()
