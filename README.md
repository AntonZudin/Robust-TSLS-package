# Robust-TSLS-package


### Installation
The current version of the package can be installed from source using devtools. 

 ```R  
 devtools::install_github("AntonZudin/Robust-TSLS-package")
```

The [online vignettes](https://antonzudin.github.io/Robust-TSLS-package/) contains paper rusults and functions documentation.

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
estimate <- robust_estimate(Y_mat_or, W_mat_or, Z, NULL, NULL, index_sub, T_0,
                                    state_names, time_column = TRUE)

tau.hat <- estimate$tau
se <- estimate$se
sprintf('robust estimator: %.2f, s.e.: %.2f,  95%% confidence interval: [%.2f, %.2f]',
          tau.hat, se, tau.hat - 1.96*se, tau.hat + 1.96*se)
```

#### References
Dmitry Arkhangelsky and Vasily Korovkin.

<b>On Policy Evaluation with Aggregate Time-Series Shocks</b>, 2024.
[<a href="https://arxiv.org/abs/1905.13660">arxiv</a>]
