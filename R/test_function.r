test_function <- function(n,T,T_0, pi_unit,theta_w,theta_y, sigma_rf, sigma_fs,
				tau,rho_agg,Z_fit, S = 200){
				
	
	### Data generation


    index_b <- sample(length(pi_unit),n,replace = TRUE)
    pi_unit_b <- pi_unit[index_b]
    theta_y_b <- theta_y[index_b]
    theta_w_b <- theta_w[index_b]
    add_shock <- pi_unit_b # 0.3*pi_unit_b + sqrt(1-0.3^2)*rnorm(n)*sd(pi_unit_b)

	art_Z <- as.numeric(simulate(Z_fit, T))
  art_Z_add <- c(rep(0, T_0),as.numeric(simulate(Z_fit, T-T_0)))
	art_H <- sqrt((1-rho_agg^2))*as.numeric(simulate(Z_fit,T)) +rho_agg*art_Z
	
	noise_fs <- sigma_fs*rmvnorm(n,sigma = diag(T))
	noise_rf <- sigma_rf*rmvnorm(n,sigma = diag(T))
	
	baseline_W <- noise_fs
	baseline_Y <- noise_rf #+0.5*outer(add_shock,art_Z_add)
	
	agg_shocks_W <- pi_unit_b%*%t(art_Z) + theta_w_b%*%t(art_H)
	agg_shocks_Y <- tau*agg_shocks_W +	theta_y_b%*%t(art_H)
	
	W_b <- baseline_W + agg_shocks_W 
	Y_b <- baseline_Y + agg_shocks_Y
	
	### Weights

	D_b <- (1/T_0)*W_b[,1:T_0]%*%(art_Z[1:T_0]-mean(art_Z[1:T_0]))/var(art_Z[1:T_0])
	
	omega_rob_b <- weights_function(Y_b[,1:T_0],W_b[,1:T_0],art_Z[1:T_0],D_unit =D_b, X_unit =matrix(1,nrow = n, ncol = 1), psi=matrix(1,nrow = T_0, ncol = 1),lambda = 'basic')
	omega_or_b <- weights_function(Y_b[,1:T_0],W_b[,1:T_0],art_Z[1:T_0],D_unit =D_b, X_unit =matrix(1,nrow = n, ncol = 1), psi=matrix(1,nrow = T_0, ncol = 1),lambda = 100000)
	
	### Aggregates
	
	agg_Y_rob_b <- as.numeric(t(omega_rob_b)%*%Y_b[,(T_0+1):T])/n
	agg_W_rob_b <- as.numeric(t(omega_rob_b)%*%W_b[,(T_0+1):T])/n

	agg_Y_or_b <- as.numeric(t(omega_or_b)%*%Y_b[,(T_0+1):T])/n
	agg_W_or_b <- as.numeric(t(omega_or_b)%*%W_b[,(T_0+1):T])/n
	
	### Coefficients
	
	
	pi_or <- lm(agg_W_or_b~art_Z[(T_0+1):T])$coefficients[2]
	delta_or <- lm(agg_Y_or_b~art_Z[(T_0+1):T])$coefficients[2]
	tau_or <- delta_or/pi_or


	
	pi_rob <- lm(agg_W_rob_b~art_Z[(T_0+1):T])$coefficients[2]
	delta_rob <- lm(agg_Y_rob_b~art_Z[(T_0+1):T])$coefficients[2]
	tau_rob <- delta_rob/pi_rob


	### Results

			
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
	test_rob <- (tau_rob - tau) / se_rob
	
	
	se_or <- sd(results_sd[,2])
	test_or <- (tau_or -tau)/se_or
	
	results <- c(tau_rob-tau, tau_or-tau, test_rob, test_or)
	
	return(results)
}