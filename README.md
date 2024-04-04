# Robust-TSLS-package


### Installation
The current version of the package can be installed from source using devtools. 

 ```R  
 devtools::install_github("AntonZudin/Robust-TSLS-package")
```

### Estimation

```R
library(RobustTSLS)
data('nak_out_data')

Y_mat_or <- data_nak[[1]][-1,]
W_mat_or <- data_nak[[2]][-1,]
Z <- data_nak[[3]]
start_year <- 1968
state_names <- data_nak[[4]][-1]
T <- dim(W_mat_or)[2]
T_0 <- 10

pi_unit <- W_mat_or%*%(Z-mean(Z)) / var_biased(Z)/T #calculate D_i

index_sub <- pi_unit >= quantile(pi_unit, 0.04) & pi_unit <= quantile(pi_unit, 1)#drop inappropriate states

robust_estimates <- robust_estimate(Y_mat_or, W_mat_or, Z, NULL, NULL, index_sub, T_0,
                                    state_names, time_column = TRUE) 

plot_3(robust_estimates)
plot_1(robust_estimates)
plot_2(robust_estimates, T_0, 1968)
```

### Simulation

``` R
library(RobustTSLS)
data('nak_out_data')

Y_mat_or <- data_nak[[1]]
W_mat_or <- data_nak[[2]]
Z <- data_nak[[3]]
start_year <- 1968
state_names <- data_nak[[4]]
T <- dim(W_mat_or)[2]
T_0 <- 10

simulation(Y_mat_or, W_mat_or, Z, share_t = 1/3, share_rank = 1/3, 
           rho_agg = 0.5, rho_theta_w = 0.2, rho_theta_y = 0.3, 
           B = 1000, S = 300, test = FALSE, K = 300, deg = 4, 
           height = 9*0.75, width = 16*0.75, plot_folder = NULL, 
           save_sim = TRUE, sim_folder = NULL, sim_name = 'simulation_result', 
           seed = 1234)

```

### Simulation resultes from file

```R
sim_from_file(file_folder = NULL, file_name = 'simulation_result',
              plot_folder = NULL, K = 300, deg = 4, 
              height = 9 * 0.75, width = 16 * 0.75)
```
