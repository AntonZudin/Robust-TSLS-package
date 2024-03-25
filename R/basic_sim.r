basic_sim <- function(F_mat_W, L_mat_W, F_mat_Y, L_mat_Y, cov_mat_fs,cov_mat_rf,pi_unit,theta_w,theta_y,
				tau,rho_agg, rho_cross,T_0,Z_fit, test = FALSE, S = 300){
				
	
	### Data generation 	
				
				
	n <- dim(F_mat_Y)[1]
	T <- dim(F_mat_Y)[2]

	art_Z <- as.numeric(simulate(Z_fit, T))
	art_H <- sqrt((1-rho_agg^2))*as.numeric(simulate(Z_fit,T)) +rho_agg*art_Z
	
	noise_fs <- rmvnorm(n,sigma = cov_mat_fs)
	noise_rf <- rmvnorm(n,sigma = cov_mat_rf)*sqrt((1-rho_cross^2)) + noise_fs*rho_cross 
	
	baseline_W <- F_mat_W + L_mat_W + noise_fs
	baseline_Y <- F_mat_Y + L_mat_Y + noise_rf
	
	agg_shocks_W <- pi_unit%*%t(art_Z) + theta_w%*%t(art_H)
	agg_shocks_Y <- tau*agg_shocks_W +	theta_y%*%t(art_H)
	
	W_b <- baseline_W + agg_shocks_W 
	Y_b <- baseline_Y + agg_shocks_Y
	
	### Weights

	D_b <- (1/T_0) * W_b[,1:T_0] %*% (art_Z[1:T_0] - mean(art_Z[1:T_0])) / var(art_Z[1:T_0])
	
	omega_rob_b <- weights_function(Y_b[,1:T_0],W_b[,1:T_0],art_Z[1:T_0],D_unit =D_b, X_unit =matrix(1,nrow = n, ncol = 1),
	                                psi=matrix(1,nrow = T_0, ncol = 1),lambda = 'basic')
	omega_or_b <- weights_function(Y_b[,1:T_0],W_b[,1:T_0],art_Z[1:T_0],D_unit =D_b, X_unit =matrix(1,nrow = n, ncol = 1), 
	                                psi=matrix(1,nrow = T_0, ncol = 1),lambda = 100000)
	
	### Aggregates
	
	agg_Y_rob_b <- as.numeric(t(omega_rob_b) %*% Y_b[,(T_0+1):T]) / n
	agg_W_rob_b <- as.numeric(t(omega_rob_b) %*% W_b[,(T_0+1):T]) / n

	agg_Y_or_b <- as.numeric(t(omega_or_b) %*% Y_b[,(T_0+1):T]) / n
	agg_W_or_b <- as.numeric(t(omega_or_b) %*% W_b[,(T_0+1):T]) / n
	
	### Coefficients
	
	
	pi_or <- lm(agg_W_or_b~art_Z[(T_0+1):T])$coefficients[2]
	delta_or <- lm(agg_Y_or_b~art_Z[(T_0+1):T])$coefficients[2]
	tau_or <- delta_or/pi_or


	
	pi_rob <- lm(agg_W_rob_b~art_Z[(T_0+1):T])$coefficients[2]
	delta_rob <- lm(agg_Y_rob_b~art_Z[(T_0+1):T])$coefficients[2]
	tau_rob <- delta_rob/pi_rob


	### Results
	
	pi_rob_true <- mean(omega_rob_b*pi_unit)
	delta_rob_true <- tau*mean(omega_rob_b*pi_unit)
	pi_or_true <- mean(omega_or_b*pi_unit)
	delta_or_true <- tau*mean(omega_or_b*pi_unit)
	
	
	
	dif_or_pi <- pi_or - pi_or_true
	dif_or_delta <- delta_or - delta_or_true
	dif_rob_pi <- pi_rob - pi_rob_true
	dif_rob_delta <- delta_rob - delta_rob_true
	dif_tau_rob <- tau_rob - tau
	dif_tau_or <- tau_or - tau
	
	
	results <- as.numeric(c(dif_rob_pi,dif_or_pi,
							dif_rob_delta, dif_or_delta,
							dif_tau_rob, dif_tau_or
	))
	
	
	if (test == TRUE){
			
		art_Z_ts <- ts(art_Z, start=1, end=T, frequency=1) 
		art_Z_fit <- auto.arima(art_Z_ts)
		art_Z_dem <- art_Z[(T_0+1):T]-mean(art_Z[(T_0+1):T])
	
		results_sd <- matrix(0, ncol = 2, nrow = S)
	
		for(j in 1:S){
	
			art_Z_j <- as.numeric(simulate(art_Z_fit, T-T_0))
			art_Z_dem_j <- art_Z_j-mean(art_Z_j)
			rob_j <- mean((agg_Y_rob_b - tau_rob*agg_W_rob_b)*art_Z_dem_j)/(pi_rob*var(art_Z_dem_j))
			or_j <-  mean((agg_Y_or_b - tau_or*agg_W_or_b)*art_Z_dem_j)/(pi_or*var(art_Z_dem_j))
			results_sd[j,] <- c(rob_j,or_j)	
		}
	
	
	se_rob <- sd(results_sd[,1])
	test_rob <- (tau_rob -tau)/se_rob
	
	
	se_or <- sd(results_sd[,2])
	test_or <- (tau_or -tau)/se_or
	
	results <- c(results, test_rob, test_or)
	
	}
	
	
	
	return(results)
}


