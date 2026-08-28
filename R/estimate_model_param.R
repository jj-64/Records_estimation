#' Estimate drift/power parameters in LDM and YNM
#'
#' Estimates the drift parameter theta in the Linear Drift Model (LDM)
#' and the gamma/power parameter in Yang-Nevzrovo model (YNM)
#' using different methods: "moments", or "MLE" using indicator series
#'
#' @param X Numeric vector of the LDM process.
#' @param method Character string, one of "moments" or "mle_indicator".
#' @param bias Logical, if biased estimator when moments estimation (default = TRUE)
#' @param model Character string, one of "LDM", "YNM"
#' @param variance logical to return variance or nor (default = TRUE)
#' @param obs_type string, one of "records" or "all" (Default = "records")
#' @param ... additional arguments specific for each model
#' @details choosing the "mle_indicator" should provide additional arguments of
#' min, max, step for the parameter grid search
#' approximate : Logical for variance.
#' scale in case of LDM Gumbel
#' For more info, see \code{\link{estimate_YNM_mle_indicator}} and \code{\link{estimate_LDM_mle_indicator}}
#' @return Numeric estimate of theta and variance (if TRUE)
#' @export
#' @examples
#' X = c(0.428,1.311,2.023,2.882,2.096,-0.197,1.339,  1.748,1.418, 0.711, 1.999,
#' 3.598, 3.308, 3.942,2.025,3.282,4.043, 0.492, 4.639, 1.408, 3.525, 5.398,
#'3.719, 3.741, 4.729)
#'
#' estimate_model_param(X, method="moments", bias = TRUE, model = "LDM", scale=1)
#' #  $param
#' # [1] 0.4462871
#'
#' #  $variance
#' # [1] 0.5625
#'
#'  estimate_model_param(X, method="moments", bias = FALSE, model = "LDM", scale=1)
#' #  $param
#' #  [1] 0.3601306
# '
#' #  $variance
#' #  [1] 0.3116952
#'
#' estimate_model_param(X, method="mle_indicator", model = "LDM", scale = 1,
#' min= 0.01, max=2, step = 0.001, approximate = FALSE)
#' # $param
#' # [1] 0.328
#'
#' # $variance
#' # [1] 0.02128002
#'
#'
#' estimate_model_param(X, method="moments", model = "YNM")
#' # $param
#' # 1.5625
#'
#' # $variance
#' # [1] 1.373291
#'
#' estimate_model_param(X, method="moments", bias = FALSE, model = "YNM")
#' # $param
#' # [1] 1.388625
#'
#' # $variance
#' # [1] 0.7497066
#'
#' estimate_model_param(X, method="mle_indicator", model = "YNM", min= 1,
#'  max=5, step = 0.001, approximate = FALSE)
#' # $param
#' # [1] 1.388
#'
#' # $variance
#' # [1] 0.0409829
estimate_model_param <- function(X, method = c("moments","mle_indicator"),
                                 bias = TRUE ,model=c("LDM","YNM"),
                                 variance = TRUE, obs_type = "records", ...) {
  method <- match.arg(method)
  model <- match.arg(model)
  args <- list(...)

  if(obs_type == "records"){
  if (model == "LDM"){
    if (method == "moments" && bias == TRUE) {
      est <- estimate_LDM_moments(X, variance = variance, scale=args$scale)

    } else if (method == "moments" && bias == FALSE) {
      est = estimate_LDM_moments_unbias(X, variance = variance, scale=args$scale)

    } else if (method == "mle_indicator") {
      est = estimate_LDM_mle_indicator(X, variance = variance, min = args$min, max=args$max, step=args$step, scale = args$scale)

    }
  } else if(model == "YNM"){
    if (method == "moments" && bias == TRUE) {
      est <- estimate_YNM_moments(X, variance = variance)

    } else if (method == "moments" && bias == FALSE) {
      est = estimate_YNM_moments_unbias(X, variance = variance)

    } else if (method == "mle_indicator") {
      est = estimate_YNM_mle_indicator(X, variance = variance, approximate=args$approximate, min = args$min, max=args$max, step=args$step)
    }
  }
  }
  return(est)
}

