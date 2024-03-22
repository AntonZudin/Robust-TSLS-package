# Robust-TSLS-package


### Installation
The current version of the package can be installed from source using devtools. 

 ```R  
 devtools::install_github("AntonZudin/Robust-TSLS-package")
```

### Example

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

robust_estimates <- robust_estimate(Y_mat_or, W_mat_or, Z, index_sub, T_0,
                                    state_names, time_column = TRUE) 

plot_3(robust_estimates)
plot_1(robust_estimates)
plot_2(robust_estimates, T_0, 1968)
```
