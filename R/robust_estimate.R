#' Compute the robust estimator for tau.
#' @param Y_mat_or an nxT matrix of outcome.
#' @param W_mat_or an nxT matrix of endogenous treatment.
#' @param Z an nx1 vector of an aggregate instrument.
#' @param index_sub an nx1 vector with binary coordinates. If one of the coordinates is TRUE, we include this state in the estimation.
#' @param T_0, the size of the learning period.
#' @param state_names an nx1 vector with state names. It is used in plot_1 function. 
#' @param time_column If true, the units are in the rows and the time periods are in the columns in Y_mat_or and W_mat_or.
#' @param D an nx1 vector of exposure.
#' @param seed, seed to set.
#' @return A robust estimator with result, attached as an attribute. It is an output of basic_analysis function. 
#' @export
#'



robust_estimate <- function(Y_mat_or, W_mat_or, Z, index_sub = NULL, T_0 = NULL, 
                            state_names = NULL, time_column = TRUE, seed = 1234){
  
  if (time_column == TRUE) {
    Y_mat_or <- as.matrix(Y_mat_or)
    W_mat_or <- as.matrix(W_mat_or)
    Z <- as.vector(Z)
  } else{
    Y_mat_or <- t(as.matrix(Y_mat_or))
    W_mat_or <- t(as.matrix(W_mat_or))
    Z <- as.vector(Z)
  }  
  
  
  
  if (!is.null(index_sub)){
    Y_mat_or <- Y_mat_or[index_sub,]
    W_mat_or <- W_mat_or[index_sub,]
  }
  
  n <- dim(Y_mat_or)[1]
  T <- dim(Y_mat_or)[2]
  
  if (is.null(T_0)){
    T_0 <- floor(T / 3)
  }
  
  
  pi_unit <- ((W_mat_or%*%(Z-mean(Z))) / T) / var_biased(Z) #calculate D_i  
  
  
  
  basic_res <- basic_analysis(Y_mat = Y_mat_or, W_mat = W_mat_or, Z = Z, 
                              X_unit = matrix(1,ncol = 1, nrow = n), psi = matrix(1,ncol = 1, nrow = T),
                              D = pi_unit, T_0 = T_0)
  
  
  
  robust_estimate <- list(result=basic_res, index_sub=index_sub, state_names=state_names)
  
  class(robust_estimate) <- "robust_estimate"
  
  
  return(robust_estimate)
  
}