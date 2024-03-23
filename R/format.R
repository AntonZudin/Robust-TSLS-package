#' Format a robust estimator object.
#' @param estimate The object to format.
#' @method format robust_estimate
#' @export


format.robust_estimate <-  function(estimate) {
  se <-  estimate$result[[3]][['se_rob_cor']]
  tau <-  estimate$result[[2]][2, 3]
  sprintf('robust estimator: %.2f, s.e.: %.2f,  95%% confidence interval: [%.2f, %.2f]',
          tau, se, tau - 1.96*se, tau + 1.96*se)
  
}