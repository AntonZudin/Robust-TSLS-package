#' Compute the robust estimator for tau.
#' @param Y_mat_or an nxT matrix of outcome.
#' @param W_mat_or an nxT matrix of endogenous treatment.
#' @param Z an nx1 vector of an aggregate instrument.
#' @param unit_covariates an optional nxC matrix of unit covariate(s), where C is the number of covariates. NULL stands for no covariates.
#' @param time_covariates an optional TxC matrix of time covariate(s), where C is the number of covariates. NULL stands for no covariates. 
#' @param index_sub an nx1 vector with binary coordinates. If the coordinate is TRUE, we include this state in the estimation. NULL stands for selecting all the units.
#' @param T_0 the size of the learning period.
#' @param unit_names an nx1 vector with unit names. It is used only in plot_1 function. 
#' @param time_column If true, the units are in the rows and the time periods are in the columns in Y_mat_or and W_mat_or.
#' @param D an nx1 vector of exposure.
#' @param seed seed to set.
#' @return A robust estimator (S3 class) with result, attached as an attribute. It is an output of basic_analysis function. 
#' @export
#'



robust_estimate <- function(Y_mat_or, W_mat_or, Z, 
                            unit_covariates = NULL, time_covariates = NULL,
                            index_sub = NULL, T_0 = NULL,
                            unit_names = NULL, time_column = TRUE, 
                            seed = 1234){
  
  if (time_column == TRUE) {
    Y_mat_or <- as.matrix(Y_mat_or)
    W_mat_or <- as.matrix(W_mat_or)
    Z <- as.vector(Z)
  } else {
    Y_mat_or <- t(as.matrix(Y_mat_or))
    W_mat_or <- t(as.matrix(W_mat_or))
    Z <- as.vector(Z)
  }  
  
  
  if (!is.null(index_sub)){
    Y_mat_or <- Y_mat_or[index_sub, ]
    W_mat_or <- W_mat_or[index_sub, ]
  } else {
    index_sub <- rep(TRUE, length.out = nrow(Y_mat_or))
  }
  
  n <- dim(Y_mat_or)[1]
  T <- dim(Y_mat_or)[2]
  
  if (is.null(T_0)){
    T_0 <- floor(T / 3)
  }
  
  
  pi_unit <- ((W_mat_or %*% (Z-mean(Z))) / T) / var_biased(Z) #calculate D_i  
  
  
  
  basic_res <- basic_analysis(Y_mat = Y_mat_or, W_mat = W_mat_or, Z = Z, 
                              unit_covariates = unit_covariates, time_covariates = time_covariates,
                              add_const = TRUE,
                              D = pi_unit, T_0 = T_0, seed)
  
  
  
  
  robust_estimate <- list(result = basic_res,
                          index_sub = index_sub, 
                          unit_names = unit_names,
                          tau = basic_res[[2]][2, 3],
                          se = basic_res[[3]][['se_rob_cor']])
  
  class(robust_estimate) <- "robust_estimate"
  
  return(robust_estimate)
  
}