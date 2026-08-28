##1) Recommended (proper) method  - observed information / Hessian-------------

#' Robust inversion of Hessian matrix
#' @param M matrix object
#' @param ridge very small error number, default = 1e-8
#' @param max_attempts integer, maximum number of attempts, default = 5
invert_posdef_matrix <- function(M, ridge = 1e-8, max_attempts = 5) {
  if (!is.matrix(M)) stop("M must be a matrix")
  n <- nrow(M)
  if (any(!is.finite(M))) stop("Non-finite entries in matrix")
  attempt <- 0
  cur_ridge <- 0
  while (attempt < max_attempts) {
    attempt <- attempt + 1
    try({
      if (cur_ridge > 0) {
        M2 <- M + diag(rep(cur_ridge, n))
      } else {
        M2 <- M
      }
      inv <- solve(M2)
      return(inv)
    }, silent = TRUE)
    cur_ridge <- if (cur_ridge == 0) ridge else cur_ridge * 10
  }
  stop("Failed to invert observed information matrix (not positive-definite).")
}

#' Compute numeric Hessian of log-likelihood (logLik_fun returns scalar log-lik)
#' @details
#' This is a numerical computation of the variance of MLE estimated paraemters.
#' For computations assuming independence, see \code{\link{compute_vcov_marginal}}
#'
#' @param logLik_fun function, loglikelihood function stored in the registry with embeded model,
#' obs_type, data and dist
#' @param params list or named vector, parameters of the likelihood function
#' @param method character "numDeriv" or "fd" for the hessian matrix derivation
#' @param eps numeric, epsilon very small number
#' @return variance-coavriance matrix
#' @examples
#' \dontrun{
#'  data = list(rec_values = c(1,2,4,6),
#'              rec_times =  c(1,3,6,7),
#'                time = 7)
#'  params = list(theta= 0.3, scale=2, shape=2)
#' logLik_wrapper <- function(params) logLik_records("LDM", "records", "frechet",data=data, params)
#' vcov_num <- compute_vcov_loglik(logLik_fun=logLik_wrapper, params = params, method = "numDeriv")
#' vcov_df <- compute_vcov_loglik(logLik_fun=logLik_wrapper, params = params, method = "fd")
#' vcov_num
#'              theta         scale         shape
#' # theta  0.0002413388 -0.0018339109 -9.914805e-04
#' # scale -0.0018339109 -0.0005549926  2.128553e-02
#' # shape -0.0009914805  0.0212855317-5.246428e-05
#' # attr(,"H")
#' #          [,1]       [,2]       [,3]
#' # [1,] -14516.5394 -679.30707 -1268.4226
#' # [2,]   -679.3071  -31.90424  -106.3396
#' # [3,]  -1268.4226 -106.33961  -112.0569
#'
#' vcov_df
#'               theta         scale         shape
#' # theta  0.0002338340 -0.0018340020 -0.0009199346
#' # scale -0.0018340020 -0.0002522914  0.0210677989
#' # shape -0.0009199346  0.0210677989 -0.0005627663
#' # attr(,"H")
#' #      [,1]  [,2]  [,3]
#' # [1,] -14520 -668 -1272
#' # [2,]   -668  -32  -106
#' # [3,]  -1272 -106  -112
#' }
#' @export
compute_vcov_loglik <- function(logLik_fun, params, method = c("numDeriv", "fd"),
                                     eps = sqrt(.Machine$double.eps)) {
  method <- match.arg(method)
  params = unlist(params)
  if (method == "numDeriv" && requireNamespace("numDeriv", quietly = TRUE)) {
    H <- numDeriv::hessian(func = logLik_fun, x = params)
  } else {
    # simple central difference Hessian (fallback)
    n <- length(params)
    H <- matrix(0, n, n)
    h <- eps * pmax(1, abs(params))
    for (i in seq_len(n)) {
      e_i <- rep(0, n); e_i[i] <- 1
      for (j in i:n) {
        e_j <- rep(0, n); e_j[j] <- 1
        fpp <- logLik_fun(params + h[i]*e_i + h[j]*e_j)
        fpm <- logLik_fun(params + h[i]*e_i - h[j]*e_j)
        fmp <- logLik_fun(params - h[i]*e_i + h[j]*e_j)
        fmm <- logLik_fun(params - h[i]*e_i - h[j]*e_j)
        H[i,j] <- (fpp - fpm - fmp + fmm) / (4*h[i]*h[j])
        H[j,i] <- H[i,j]
      }
    }
  }

  # observed information is -H; invert it
  obs_info <- -H
  colnames(obs_info) = names(params)
  rownames(obs_info) = names(params)
  # regularized inversion helper below
  vcov <- invert_posdef_matrix(obs_info)
  attr(vcov, "H") <- H
  return(vcov)
}


## 2 Build diagonal vcov from named per-param variance functions assuming independence-----
#' Compute var-cov matrix assuming independence
#' @details
#' This is a numerical computation of the variance of MLE estimated paraemters.
#' For numerical computations, see \code{\link{compute_vcov_loglik}}
#'
#' @param var_fun list of functions, each function stored in the registry for a defined model,
#' obs_type, dist and param_name
#' @param data list or dataframe, named with "rec_values", "rec_times" and "time"
#' @param params list or named vector, parameters of the likelihood function
#' @param param_names default takes the names of input vector params
#' @return variance-coavriance matrix
#' @examples
#' \dontrun{
#'  data = list(rec_values = c(1,2,4,6),
#'              rec_times =  c(1,3,6,7),
#'                time = 7)
#'  params = list(theta= 0.3, scale=2, shape=2)
#' var_funs <- list(
#'   theta = var_logLik_records("LDM","records","frechet","theta"),
#'   scale = var_logLik_records("LDM","records","frechet","scale"),
#'   shape = var_logLik_records("LDM","records","frechet","shape")
#' )
#' vcov_marg <- compute_vcov_marginal(var_funs, data, params)
#' #          theta        scale       shape
#' # theta 6.888694e-05 0.0000000000 0.000000000
#' # scale 0.000000e+00 0.0007125534 0.000000000
#' # shape 0.000000e+00 0.0000000000 0.002413605
#' }
#' @export
compute_vcov_marginal <- function(var_fun, data, params, param_names = names(params)) {
  p <- length(param_names)
  vcov <- matrix(0, p, p)
  rownames(vcov) <- colnames(vcov) <- param_names
  for (i in seq_along(param_names)) {
    nm <- param_names[i]
    if (!is.null(var_fun[[nm]])) {
      vcov[i,i] <- var_fun[[nm]](data = data, params = params)
    } else {
      warning("No variance function for ", nm, " leaving NA")
      vcov[i,i] <- NA_real_
    }
  }
  return(vcov)
}

## 3) Combine SEs for functions of parameters (Delta method) SE for a scalar function g(par)
# Example
# g_scale <- function(par) 1/as.numeric( par["scale"])  # par order theta,A,a
# se_scale <- se_delta(g_fun = g_scale, par = params, vcov = vcov_marg)
se_delta <- function(g_fun, par, vcov, method = c("numDeriv", "fd")) {
  method <- match.arg(method)
  par = unlist(par)
  if (method == "numDeriv" && requireNamespace("numDeriv", quietly = TRUE)) {
    grad <- numDeriv::grad(func = g_fun, x = par)
  } else {
    # tiny forward difference gradient
    eps <- sqrt(.Machine$double.eps)
    n <- length(par)
    grad <- numeric(n)
    h <- eps * pmax(1, abs(par))
    for (i in seq_len(n)) {
      e <- rep(0, n); e[i] <- 1
      grad[i] <- (g_fun(par + h[i]*e) - g_fun(par - h[i]*e)) / (2*h[i])
    }
  }
  var_g <- as.numeric(t(grad) %*% vcov %*% grad)
  sqrt(max(0, var_g))
}

# # compute vcov (numeric) ----------------------
#' master helper: try Hessian-based vcov, fallback to diag of provided var_fun
#'
#' @param logLik_fun function, loglikelihood function stored in the registry with embeded model,
#' obs_type, data and dist
#' @param var_fun list of functions, each function stored in the registry for a defined model,
#' obs_type, dist and param_name
#' @param data list or dataframe, named with "rec_values", "rec_times" and "time"
#' @param params list or named vector, parameters of the likelihood function
#' @param quiet default is FALSE, for quiet computation
#' @examples
#' \dontrun{
#'  data = list(rec_values = c(1,2,4,6),
#'              rec_times =  c(1,3,6,7),
#'                time = 7)
#'  params = list(theta= 0.3, scale=2, shape=2)
#' logLik_wrapper <- function(params) logLik_records("LDM", "records", "frechet",
#' data=data, params)
#' var_funs <- list(
#'   theta = var_logLik_records("LDM","records","frechet","theta"),
#'   scale = var_logLik_records("LDM","records","frechet","scale"),
#'   shape = var_logLik_records("LDM","records","frechet","shape")
#' )
#' compute_vcov_mle(logLik_wrapper, params, var_funs, data)
#' }
#' @export
compute_vcov_mle <- function(logLik_fun = NULL, params, var_fun = NULL,
                              data, quiet = FALSE) {
  # 1) Try numeric Hessian from logLik_fun
  if (!is.null(logLik_fun)) {
    try({
      vcov <- compute_vcov_loglik(logLik_fun, params)
      if (!quiet) message("Using Hessian-based vcov.")
      return(vcov)
    }, silent = TRUE)
    if (!quiet) message("Hessian-based vcov failed or not available; trying marginal variances.")
  }

  # 2) Fallback: use var_fun (diagonal)
    param_names <- names(params)
    vcov <- compute_vcov_marginal(var_fun, data, params = params, param_names = param_names)
    if (!quiet) message("Using diagonal vcov from marginal variance functions (assumes zero covariances).")
    return(vcov)


  stop("No method available to construct vcov (provide logLik_fun or var_fun + data).")
}
