#' Draw 2 time-series plots of aggregated data for original and robust weights.
#' @description
#' Solid lines represent aggregated Y_t and W_t. Dashed ones represent OLS predictions of the aggregate data with the instrument.
#' Blue and red dashed lines in the robust weights plot stand for learning and estimation periods respectively.
#' 
#' @param robust_estimate The object to use for building the plot.
#' @param T_0 The number of periods used for constructing the weights.
#' @param start_year The starting year in the plot's x-axis.
#' @param folder The folder where the plots are saved. NULL stands for the working directory.
#' @param height The height of the plot.
#' @param width The width of the plot.
#' @export
#'


plot_2 <- function(robust_estimate, T_0, start_year=1, 
                   folder = NULL, height = 9, width = 16) {
  
  basic_res <- robust_estimate$result  
  
  years <- start_year:(start_year+T-1)
  
  
  scale_w_agg <- (max(basic_res[[5]][,1])-min(basic_res[[5]][,1]))/2
  scale_y_agg <- (max(basic_res[[5]][,3])-min(basic_res[[5]][,3]))/2
  
  
  W_agg <- basic_res[[5]][,1]/scale_w_agg
  W_fit <- basic_res[[5]][,2]/scale_w_agg
  Y_agg <- basic_res[[5]][,3]/scale_y_agg
  Y_fit <- basic_res[[5]][,4]/scale_y_agg
  W_agg_rob <- basic_res[[5]][,5]/scale_w_agg
  W_fit_rob <- basic_res[[5]][,6]/scale_w_agg
  Y_agg_rob <- basic_res[[5]][,7]/scale_y_agg
  Y_fit_rob <- basic_res[[5]][,8]/scale_y_agg
  
  
  
  r_sq_w_rob <- 1 - var_biased(W_agg_rob[(T_0+1):T] - W_fit_rob[(T_0+1):T])/var_biased(W_agg_rob[(T_0+1):T])
  r_sq_y_rob <-1 - var_biased(Y_agg_rob[(T_0+1):T] - Y_fit_rob[(T_0+1):T])/var_biased(Y_agg_rob[(T_0+1):T])
  
  r_sq_w_or <- 1 - var_biased(W_agg[(T_0+1):T] - W_fit[(T_0+1):T]) / var_biased(W_agg[(T_0+1):T])
  r_sq_y_or <- 1 - var_biased(Y_agg[(T_0+1):T] - Y_fit[(T_0+1):T]) / var_biased(Y_agg[(T_0+1):T])
  
  
  if (is.null(folder)){
    pdf('nuk_original.pdf', width = width, height = height)
  } else{
    dir_1 <- paste(folder, "nuk_original.pdf", sep = "/")
    pdf(dir_1, width = width, height = height)
  }
  
  print(length(years))
  print(dim(W_agg))
  print(dim(Y_agg))
  
  par(mfrow=c(2,1)) 
  plot(years,W_agg, ylim = c(min(W_agg),max(W_agg)+0.1), type = 'b',xlab ="",ylab = 'Aggregate W',main = 'First Stage',lty = 1, pch = 19,frame = FALSE)
  lines(years,W_fit,col = 'red',lty = 2, type = 'b', pch = 18)
  abline(h = 0, lty = 2, col = 'black')
  legend('bottomright', legend = c('Aggregate data, original weights', 'OLS fit, full sample'), col = c('black','red'), pch = c(19,18))
  
  
  plot(years,Y_agg,ylim = c(min(Y_agg),max(Y_agg)), type = 'b',xlab ="", ylab = 'Aggregate Y',main = 'Reduced Form',lty = 1, pch = 19,frame = FALSE)
  lines(years,Y_fit,col = 'red',lty = 2, type = 'b', pch = 18)
  abline(h = 0, lty = 2, col = 'black')
  legend('topleft', legend = c('Aggregate data, original weights', 'OLS fit, full sample'), col = c('black','red'), pch = c(19,18))
  
  dev.off()
  
  
  
  if (is.null(folder)){
    pdf('nuk_new.pdf', width = width, height = height)
  } else{
    dir_2 <- paste(folder, "nuk_new.pdf", sep = "/")
    pdf(dir_2, width = width, height = height)
  }
  
  par(mfrow=c(2,1)) 
  
  plot(years,W_agg_rob, ylim = c(min(W_agg),max(W_agg)+0.1), type = 'b',xlab = '',ylab = 'Aggregate W',main = 'First Stage',lty = 1, pch = 19,frame = FALSE)
  lines(years[(T_0+1):T],W_fit_rob[(T_0+1):T],col = 'red',lty = 2, type = 'b', pch = 18)
  lines(years[1:T_0],W_fit_rob[1:T_0],col = 'blue',lty = 2, type = 'b', pch = 17)
  abline(h = 0, lty = 2, col = 'black')
  abline(v = start_year+T_0, lty = 2, col = 'grey')
  legend('bottomright', legend = c('Aggregate data, robust weights', 'OLS fit, learning periods', 'OLS fit, estimation periods'), col = c('black','blue','red'), pch = c(19,17,18))
  
  
  plot(years,Y_agg_rob, ylim = c(min(Y_agg),max(Y_agg)+0.1), type = 'b',xlab = '',ylab = 'Aggregate Y',main = 'Reduced Form',lty = 1, pch = 19,frame = FALSE)
  lines(years[(T_0+1):T],Y_fit_rob[(T_0+1):T],col = 'red',lty = 2, type = 'b', pch = 18)
  lines(years[1:T_0],Y_fit_rob[1:T_0],col = 'blue',lty = 2, type = 'b', pch = 17)
  abline(h = 0, lty = 2, col = 'black')
  abline(v = start_year+T_0, lty = 2, col = 'grey')
  legend('topleft', legend = c('Aggregate data, robust weights', 'OLS fit, learning periods', 'OLS fit, estimation periods'), col = c('black','blue','red'), pch = c(19,17,18))
  
  
  dev.off()	
  
}
