#' Draw a comparative scatterplot for original weights and robust weights.
#' 
#' @param robust_estimate The object to use for building the plot.
#' @param folder The folder where the plots are saved. NULL stands for the working directory.
#' @param height The height of the plot.
#' @param width The width of the plot.
#' @export
#'

plot_3 <- function(robust_estimate, folder = NULL, width = 9, height = 9) {
  
  basic_res <- robust_estimates$result
  index_sub <- robust_estimates$index_sub
  state_names <- robust_estimates$state_names
  
  or_weights <- basic_res[[1]][,1]/sd(basic_res[[1]][,1])
  rob_weights <- basic_res[[1]][,2]/sd(basic_res[[1]][,2])
  

  
  if (is.null(folder)){
    pdf('nak_weights.pdf', width = width, height = height)
  } else{
    dir <- paste(folder, "nak_weights.pdf", sep = "/")
    pdf(dir, width = width, height = height)
  }
  
  plot(or_weights,rob_weights, pch = 3, cex = 0.5, col = 'white', lwd = 1, xlab = 'Original weights', ylab = 'Robust weights')
  text(or_weights,rob_weights, labels=state_names[index_sub], cex=0.7, font=1)
  abline(a = 0, b = 1, lty = 2, lwd = 0.5, col = 'grey')
  dev.off()

}  