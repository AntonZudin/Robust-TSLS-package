#' Draw a comparative scatterplot for original weights and robust weights.
#' 
#' @param robust_estimate The object to use for building the plot.
#' @param save_pdf  If True, the scatterplot is saved and not displayed.
#' @param file If save_pdf is true, string giving the file path (including the file name).
#' @param height The height of the plot.
#' @param width The width of the plot.
#' @export
#'

weights_scatterplot <- function(robust_estimate, save_pdf = FALSE,
                                file = 'nak_weights', height = 9, width = 9) {
  
  basic_res <- robust_estimate@result
  index_sub <- robust_estimate@index_sub
  unit_names <- robust_estimate@unit_names
  
  or_weights <- basic_res[[1]][,1]/sd(basic_res[[1]][,1])
  rob_weights <- basic_res[[1]][,2]/sd(basic_res[[1]][,2])
  

  
  if (save_pdf){
    
    if (nchar(file) >= 4){
      if (substring(file, nchar(file) - 3) == ".pdf") {
        dir <- file
      } else {
        dir <- cat(file, '.pdf', sep = '')
      }
    } else{
        dir <- cat(file, '.pdf', sep = '')
    }
    
    pdf(dir, width = width, height = height) 
    plot(or_weights,rob_weights, pch = 3, cex = 0.5, col = 'white', lwd = 1, xlab = 'Original weights', ylab = 'Robust weights')
    text(or_weights,rob_weights, labels = unit_namess[index_sub], cex = 0.7, font = 1)
    abline(a = 0, b = 1, lty = 2, lwd = 0.5, col = 'grey')
    dev.off()
  } else {
    
    plot(or_weights,rob_weights, pch = 3, cex = 0.5, col = 'white', lwd = 1, xlab = 'Original weights', ylab = 'Robust weights')
    text(or_weights,rob_weights, labels = unit_names[index_sub], cex = 0.7, font = 1)
    abline(a = 0, b = 1, lty = 2, lwd = 0.5, col = 'grey')
  }
}  