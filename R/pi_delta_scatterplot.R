#' Draw 2 scatterplots of reduced-form and first-stage coefficients for original and robust weights.
#' @param robust_estimate The object to use for building the plot.
#' @param folder The folder where the plots are saved. NULL stands for the working directory.
#' @param height The height of the plot.
#' @param width The width of the plot.
#' @export
#'



pi_delta_scatterplot <- function(robust_estimate, folder = NULL, 
                                 save_pdf = FALSE, height = 9, width = 9) {
  
  basic_res <- robust_estimate$result
  
  
  tau_or <- basic_res[[2]][1,3]
  #se_or <- basic_res[[3]][1]
  se_or_cor <- basic_res[[3]][1]
  tau_rob <- basic_res[[2]][2,3]
  #se_rob <- basic_res[[3]][2]
  se_rob_cor <- basic_res[[3]][2]
  #se_alt_cor <- basic_res[[3]][3]
  #se_alt_cor/se_rob_cor-1
  
  
  
  ## Weights
  or_weights <- basic_res[[1]][,1]/sd(basic_res[[1]][,1])
  rob_weights <- basic_res[[1]][,2]/sd(basic_res[[1]][,2])
  
  
  
  scale_fs <- (max(basic_res[[4]][,2]) - min(basic_res[[4]][,2]))/2
  
  rf_coeffs_full <- basic_res[[4]][,1]/scale_fs
  fs_coeffs_full <- basic_res[[4]][,2]/scale_fs
  
  
  rf_coeffs_post <- basic_res[[4]][,3]/scale_fs
  fs_coeffs_post <- basic_res[[4]][,4]/scale_fs
  
  ### tau_or for post data
  tau_int <- mean(rf_coeffs_post*or_weights)/mean(fs_coeffs_post*or_weights)
  
  
  
  omega_or_norm <- basic_res[[1]][,1]/mean(abs(basic_res[[1]][,1]))
  omega_rob_norm <- basic_res[[1]][,2]/mean(abs(basic_res[[1]][,1]))
  
  
  index_pos_or <- omega_or_norm>=0
  index_pos_rob <- omega_rob_norm>=0
  
  range_rf <- c(min(rf_coeffs_full,rf_coeffs_post),max(rf_coeffs_full,rf_coeffs_post ))
  range_fs <- c(min(fs_coeffs_full,fs_coeffs_post),max(fs_coeffs_full,fs_coeffs_post ))
  
  
  ### Computing coordinates for blue triangles
  x_or_p <- mean(abs(omega_or_norm[index_pos_or])*fs_coeffs_full[index_pos_or]) / mean(abs(omega_or_norm[index_pos_or]))
  x_or_n <- mean(abs(omega_or_norm[!index_pos_or])*fs_coeffs_full[!index_pos_or]) / mean(abs(omega_or_norm[!index_pos_or]))
  y_or_p <- mean(abs(omega_or_norm[index_pos_or])*rf_coeffs_full[index_pos_or]) / mean(abs(omega_or_norm[index_pos_or]))
  y_or_n <- mean(abs(omega_or_norm[!index_pos_or])*rf_coeffs_full[!index_pos_or]) / mean(abs(omega_or_norm[!index_pos_or]))
  
  x_p <- mean(abs(omega_rob_norm[index_pos_rob])*fs_coeffs_post[index_pos_rob]) / mean(abs(omega_rob_norm[index_pos_rob]))
  x_n <- mean(abs(omega_rob_norm[!index_pos_rob])*fs_coeffs_post[!index_pos_rob]) / mean(abs(omega_rob_norm[!index_pos_rob]))
  y_p <- mean(abs(omega_rob_norm[index_pos_rob])*rf_coeffs_post[index_pos_rob]) / mean(abs(omega_rob_norm[index_pos_rob]))
  y_n <- mean(abs(omega_rob_norm[!index_pos_rob])*rf_coeffs_post[!index_pos_rob]) / mean(abs(omega_rob_norm[!index_pos_rob]))
  
  ## Plot 1
  if (save_pdf){
    if (is.null(folder)){
      pdf('nuk_points_1.pdf', width = width, height = height)
    } else{
      dir_1 <- paste(folder, "nuk_points_1.pdf", sep = "/")
      pdf(dir_1, width = width, height = height)
    }
    
    plot(x=fs_coeffs_full, y=rf_coeffs_full, ylim = range_rf,xlim = range_fs, pch = 1, frame = FALSE,
         bg=NULL , col=1 + as.numeric(omega_or_norm>0), cex = abs(omega_or_norm),
         xlab = 'First Stage Coefficients', ylab = 'Reduced Form Coefficients')
    abline(a = weighted.mean(rf_coeffs_full - tau_or*fs_coeffs_full, w = abs(omega_or_norm)), b =tau_or, lty = 2, col = 'grey')
    points( cbind(c(x_or_p,x_or_n),c(y_or_p,y_or_n)),cex = 2,pch = 2, col = 'blue' )
    legend( x= x_p, y = y_n-0.3, legend = substitute(paste(hat(tau)[TSLS], ' = ', tau_or, ', ', hat(se)(hat(tau)[TSLS]), ' = ', se_or_cor),
                                                     list(tau_or = round(tau_or,2), se_or_cor = round(se_or_cor,2))), cex=1)
    dev.off()
  } else {
    plot(x=fs_coeffs_full, y=rf_coeffs_full, ylim = range_rf,xlim = range_fs, pch = 1, frame = FALSE,
         bg=NULL , col=1 + as.numeric(omega_or_norm>0), cex = abs(omega_or_norm),
         xlab = 'First Stage Coefficients', ylab = 'Reduced Form Coefficients')
    abline(a = weighted.mean(rf_coeffs_full - tau_or*fs_coeffs_full, w = abs(omega_or_norm)), b =tau_or, lty = 2, col = 'grey')
    points( cbind(c(x_or_p,x_or_n),c(y_or_p,y_or_n)),cex = 2,pch = 2, col = 'blue' )
    legend( x= x_p, y = y_n-0.3, legend = substitute(paste(hat(tau)[TSLS], ' = ', tau_or, ', ', hat(se)(hat(tau)[TSLS]), ' = ', se_or_cor),
                                                     list(tau_or = round(tau_or,2), se_or_cor = round(se_or_cor,2))), cex=1)
  }
  
  ## Plot 2
  if (save_pdf){
    if (is.null(folder)){
      pdf('nuk_points_2.pdf', width = width, height = height)
    } else{
      dir_2 <- paste(folder, "nuk_points_2.pdf", sep = "/")
      pdf(dir_2, width = width, height = height)
    }
    
    plot(x=fs_coeffs_post, y=rf_coeffs_post, pch = 1,ylim = range_rf,xlim = range_fs, frame = FALSE,
         bg=NULL , col=1 + as.numeric(omega_rob_norm>0), cex = abs(omega_rob_norm),
         xlab = 'First Stage Coefficients', ylab = 'Reduced Form Coefficients')
    abline(a = weighted.mean(rf_coeffs_post - tau_rob*fs_coeffs_post,w = abs(omega_rob_norm)), b =tau_rob,lty = 2, col = 'grey')
    points( cbind(c(x_p,x_n),c(y_p,y_n)),cex = 2,pch = 2, col = 'blue' )
    legend( x= x_p+0.3, y = y_n-0.3, legend = substitute(paste(hat(tau)[rob], ' = ', tau_rob, ', ', hat(se)(hat(tau)[rob]), ' = ', se_rob_cor),
                                                         list(tau_rob = round(tau_rob,2), se_rob_cor = round(se_rob_cor,2))), cex=1)
  
    dev.off()
  } else{
    plot(x=fs_coeffs_post, y=rf_coeffs_post, pch = 1,ylim = range_rf,xlim = range_fs, frame = FALSE,
         bg=NULL , col=1 + as.numeric(omega_rob_norm>0), cex = abs(omega_rob_norm),
         xlab = 'First Stage Coefficients', ylab = 'Reduced Form Coefficients')
    abline(a = weighted.mean(rf_coeffs_post - tau_rob*fs_coeffs_post,w = abs(omega_rob_norm)), b =tau_rob,lty = 2, col = 'grey')
    points( cbind(c(x_p,x_n),c(y_p,y_n)),cex = 2,pch = 2, col = 'blue' )
    legend( x= x_p+0.3, y = y_n-0.3, legend = substitute(paste(hat(tau)[rob], ' = ', tau_rob, ', ', hat(se)(hat(tau)[rob]), ' = ', se_rob_cor),
                                                         list(tau_rob = round(tau_rob,2), se_rob_cor = round(se_rob_cor,2))), cex=1)
  }
  

}