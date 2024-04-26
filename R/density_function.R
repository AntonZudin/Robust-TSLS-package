#' Computations for density plots.
#' @param x a vector of observed.
#' @param K number of domain splits.
#' @param deg degrees of freedom for natural cubic splines.
#' @export
#'



density_function <- function(x, K, deg){
  
  x_min <- min(x)
  x_max <- max(x)
  range_x <- x_max - x_min
  low_x <- x_min - 0.2*range_x
  up_x <- x_max + 0.2*range_x
  range_full <-up_x- low_x
  splits <- seq(low_x,up_x,length.out = K+1)
  mesh_size <- splits[2] - splits[1]
  centers <- (splits[-1]+splits[-(K+1)])/2
  counts <- as.vector(table(cut(x,splits,include.lowest = TRUE)))
  scale <- sum(counts)*mesh_size
  
  data_matrix <- splines::ns(centers, df = deg)
  pois_reg_res <- glm(counts~data_matrix, family = 'poisson')
  freq_pois <- exp(pois_reg_res$linear.predictors)
  dens_pois <- freq_pois / scale
  
  
  return(cbind(centers,freq_pois,dens_pois))
}