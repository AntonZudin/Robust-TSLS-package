#' Biased covarience.
#' 
#' @description
#' The denominator is n.
#'   
#' @param x a numeric vector, matrix or data frame.
#' @param y `NULL` (default) or a vector, matrix or data frame with compatible dimensions to x. The default is equivalent to y = x (but more efficient).  
#' @export
#'

cov_biased <- function(x, y){
  return((length(x) - 1) / length(x) * cov(x, y)) 
  
}


#' Biased varience.
#' 
#' @description
#' The denominator is n.
#'   
#' @param x a numeric vector, matrix or data frame.
#' @export
#'


var_biased <- function(x){
  return((length(x) - 1) / length(x) * var(x)) 
}