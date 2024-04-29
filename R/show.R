#' Show a robust_estimate object.
#' @param object The object to print.
#' @method print robust_estimate.
#' @export
#'


setMethod("show",
          "robust_estimate",
          function(object) {
            se <-  object@result[[3]][['se_rob_cor']]
            tau <-  object@result[[2]][2, 3]
            output <- sprintf('robust estimator: %.2f, s.e.: %.2f,  95%% confidence interval: [%.2f, %.2f]',
                    tau, se, tau - 1.96*se, tau + 1.96*se)
            cat(output, "\n")
          }
)

  