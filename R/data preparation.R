#' Data preparation for simulations.
#'
#' @param data_mat a matrix.
#' @param rank the rank of the low-rank matrix approximation with SVD.
#' @return A list with generalized fixed effects (low-rank matrix), two-way fixed effects matrix, residuals matrix, unit fixed effects and factor loadings. 
#' @export
#'


data_preparation <- function(data_mat, rank){
  
  n <- dim(data_mat)[1]
  T <- dim(data_mat)[2]
  
  svd_data_mat <- svd(data_mat)
  factor_unit <- as.matrix(svd_data_mat$u[,1:rank] * sqrt(n))
  factor_time <- as.matrix(svd_data_mat$v[,1:rank] * sqrt(T))
  
  magnitude <- svd_data_mat$d[1:rank] / sqrt(n * T)
  L_mat_orig <- factor_unit %*% diag(magnitude, nrow = rank, ncol = rank) %*%
    t(factor_time) #matrix approximation
  
  error_mat <- data_mat - L_mat_orig
  F_mat <- outer(rowMeans(L_mat_orig), rep(1,T)) +
    outer(rep(1,n), colMeans(L_mat_orig)) - mean(L_mat_orig)
  L_mat <- L_mat_orig - F_mat
  #unit_fe <- rowMeans(F_mat)
  
  return(list(L_mat, F_mat, error_mat))
}