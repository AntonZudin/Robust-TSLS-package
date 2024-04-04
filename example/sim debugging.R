library(mvtnorm)
library('forecast')
library(xtable)
library(plotrix)



rm(list = ls())
set.seed(1234)


#setwd('/Users/darkhangelsky/Library/CloudStorage/Dropbox-Personal/Research/aggregate_iv/code/simulation_paper/Nakamura')
setwd("C:/Users/Serge/OneDrive - экономический факультет МГУ им. М.В.Ломоносова/Митя/Robust-TSLS-package")
source('R/simulation.R')
source('R/basic_sim.R')
source('R/weights_function.R')
source('R/data preparation.R')
source('R/density_function.R')
source('R/sim_from_file.R')


load('data/nak_out_data.RData')

Y_mat_or <- data_nak[[1]]
W_mat_or <- data_nak[[2]]
Z <- data_nak[[3]]
start_year <- 1968
state_names <- data_nak[[4]]
T <- dim(W_mat_or)[2]
T_0 <- 10

simulation(Y_mat_or, W_mat_or, Z, share_t = 1/3, share_rank = 1/3, 
           rho_agg = 0.5, rho_theta_w = 0.2, rho_theta_y = 0.3, B = 1000,
           S = 300, test = FALSE, K = 300, deg = 4, height = 9*0.75, width = 16*0.75, 
           save_sim = TRUE, plot_folder = NULL, sim_folder = NULL, 
           sim_name = 'simulation_result', seed = 1234)

sim_from_file(file_folder = NULL, file_name = 'simulation_result', plot_folder = NULL,
              K = 300, deg = 4, height = 9*0.75, width = 16*0.75)



