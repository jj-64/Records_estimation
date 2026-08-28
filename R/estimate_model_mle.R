#' Maximize LogLikelihood function to estimate model parameters
#' @param data numeric vector, observed series or list
#' @param logLik_fun function, loglikelihood function stored in the registry with embeded model,
#' obs_type, data and dist
#' @param obs_type character "all" or "records" depending on the data used
#' @param lower_bounds named numeric vector, lower bound of the model parameters
#' @param upper_bounds named numeric vector, upper bound of the model parameters
#' @param start_values named numeric vector, start value for the optim model
#' @examples
#' \dontrun{
#'   ## Loglikelihood call
#' logLik_fun <- loglik_registry[["YNM"]][["all"]][["frechet"]]
#'   ## Optimize
#' estimate_model_mle(logLik_fun, obs_type = "all", data = c(1,1.5,1,3),
#'                    lower_bounds=c(gamma = 1.001, shape = 0.001, scale = 0.001),
#'                    upper_bounds = c(gamma = 5,shape = 20,scale = 20),
#'                    start_values = c(gamma = 1.001, shape = 0.001, scale = 0.001)
#'                      )
#'
#' #  $par
#' #  gamma    shape    scale
#' #  1.742214 3.302968 0.808554
#' #
#' #  $objective
#' #  [1] 2.948769
#' #
#' #  $convergence
#' #  [1] 0
#' #
#' #  $iterations
#' #  [1] 25
#' #
#' #  $evaluations
#' #  function gradient
#' #  29       84
#' #
#' #  $message
#' #  [1] "relative convergence (4)"
#' #
#'
#'
#'   ## Loglikelihood call
#' logLik_fun_rec <- loglik_registry[["YNM"]][["records"]][["frechet"]]
#'   ## Records Data
#' data_rec = list(rec_values = c(1,5,7,10), rec_times = c(1,2,5,6), time = 10)
#' data_rec = data.frame(rec_values = c(1,5,7,10), rec_times = c(1,2,5,6), time = 10)
#'   ## Optimize
#' estimate_model_mle(logLik_fun_rec, obs_type = "records", data = data_rec,
#'                     lower_bounds=c(gamma = 1.001, shape = 0.001, scale = 0.001),
#'                     upper_bounds = c(gamma = 5,shape = 20,scale = 20),
#'                     start_values = c(gamma = 1.001, shape = 0.001, scale = 0.001)
#'                      )
#' # $par
#' #  gamma    shape    scale
#' #  1.314292 1.871232 1.312383
#' #
#' #  $objective
#' #  [1] 12.45219
#' #
#' #  $convergence
#' #  [1] 0
#' #
#' #  $iterations
#' #  [1] 21
#' #
#' #  $evaluations
#' #  function gradient
#' #  22       73
#' #
#' #  $message
#' #  [1] "relative convergence (4)"
#' }
#' @export
estimate_model_mle = function(data,
                             logLik_fun,
                             obs_type = "all",
                             lower_bounds=c(gamma = 1.001, shape = 0.001, scale = 0.001),
                             upper_bounds = c(gamma = 5,shape = 20,scale = 20),  # Example: Bound params[1] by min(R/L)
                             start_values = c(gamma = 1.001, shape = 0.001, scale = 0.001)
                             ){
  if(obs_type == "all"){
    if(!is.numeric(data)) stop("data should be a numerical vector")
    if(length(rec_values(data))<=1) {return("Impossible to solve since we have only one trivial record")}  ## only one record, ignore

  } else{
    if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")
    if(length(data$rec_values)<=1) {return("Impossible to solve since we have only one trivial record")}  ## only one record, ignore
  }

  success <- FALSE

  ##Reverse to minimize
  Likelihood = function(params){ - logLik_fun(data=data, params) }

  MLE_C= nlminb (start_values,Likelihood,  lower=lower_bounds, upper = upper_bounds)

  # Check if the result is successful (i.e., no error occurred)
  if (inherits(MLE_C$par, "try-error")){
    MLE_C <- maxLik(  logLik_fun, start=start_values)
  } else if(inherits(MLE_C$par, "try-error")){
    MLE_C=   optim(par = start_values,           # Starting values for params
                   fn = Likelihood,        # The likelihood function
                   method = "L-BFGS-B",
                   lower = lower_bounds,         # Lower bounds
                   upper = upper_bounds)
  }
  return(MLE_C)
}
