#' Print a robust_estimate object.
#' @param x The object to print.
#' @method print robust_estimate.
#' @export
#'
print.robust_estimate <- function(estimate) { 
  cat(format(estimate), "\n") 
  }


  