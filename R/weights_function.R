#' Calculate the original and robust weights  
#' @param Y_mat an nxT matrix of outcome.
#' @param W_mat an nxT matrix of endogenous treatment.
#' @param Z_agg an nx1 vector of an aggregate instrument.
#' @param D_unit an nx1 vector of exposure.
#' @param X_unit unit covariate.
#' @param psi aggregate covariates.
#' @param lambda regularization hyperparameter.
#' 
#' @return A list with original and robust weight vectors.
#' 
#' @export
#'



weights_function <- function(Y_mat, W_mat, Z_agg, D_unit, unit_covariates = NULL, time_covariates = NULL, lambda = 'basic') {
  
  T <- dim(Y_mat)[2]
  n <- dim(Y_mat)[1]
  dim_x <- dim(unit_covariates)[2] + 1
  
  
  unit_covariates = cbind(matrix(1, ncol = 1, nrow = n), matrix(unit_covariates, nrow = n))
  time_covariates = cbind(matrix(1, ncol = 1, nrow = T), matrix(time_covariates, nrow = T))
  
  
  Y_dm <- Y_mat - outer(rep(1,n),colMeans(Y_mat)) -  outer(rowMeans(Y_mat),rep(1,T)) + mean(Y_mat)
  W_dm <- W_mat - outer(rep(1,n),colMeans(W_mat)) -  outer(rowMeans(W_mat),rep(1,T)) + mean(W_mat)
  
  Z_agg <- matrix(Z_agg, nrow = T)
  Z_full <- cbind(psi,Z_agg)
  M_z <- diag(T) - Z_full%*%solve(t(Z_full) %*% Z_full) %*% t(Z_full)
  Y_z <- Y_dm %*% M_z
  W_z <- W_dm %*% M_z
  Y_norm <- Y_z/norm(Y_z-mean(Y_z),'f')
  W_norm <- W_z/norm(W_z-mean(W_z),'f')
  
  if (lambda == 'basic') {lambda <- max(norm(Y_norm,'2')^2,norm(W_norm,'2')^2)}
  
  X <- cbind(X_unit,D_unit)
  A <- cbind(Y_norm,W_norm)
  
  D <- solve(A%*%t(A) + lambda*diag(n)/T)
  weights_un <- D%*%X%*%solve(t(X)%*%D%*%X)%*%c(rep(0,dim_x),1)
  weights_norm <- weights_un/mean(weights_un*D_unit)
  return(weights_norm)
  
}



