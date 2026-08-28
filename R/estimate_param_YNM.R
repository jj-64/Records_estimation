
#########------------- Biased based on Moments -------------- ##########
#' Estimation of Gamma Parameter Based on Number of Records (with optional Variance)
#'
#' Estimates the \eqn{\gamma} parameter of the Yang-Nevzorov Model (YNM)
#' based on on Moments and number of records, and optionally computes its analytical variance.
#'
#'
#' @details
#' The estimator of \eqn{\gamma} is defined as:
#' \deqn{
#' \hat{\gamma} = \frac{1}{1 - \frac{N_T}{T}}
#' }
#' where:
#' \itemize{
#'   \item \eqn{N_T} is the number of records in the observed sequence \eqn{X_t},
#'   \item \eqn{T} is the total number of observations.
#' }
#'
#' This estimator is known to be biased for small sample sizes.
#' A bias-corrected version is available in \code{\link{estimate_YNM_moments_unbias}}.
#'
#' The analytical variance of the estimator is:
#' \deqn{
#' \text{Var}(\hat{\gamma}) = \gamma^2 (\gamma - 1)
#' }
#'
#' @param X Numeric vector representing the time series or dataset.
#' @param variance Logical; if TRUE (default), compute the analytical variance.
#'
#' @return A list containing:
#' \item{gamma}{The estimated \eqn{\gamma} parameter.}
#' \item{variance}{The variance of the estimator (or NA if \code{variance = FALSE}).}
#'
#' @examples
#' X=YNM_series(T=50, gamma=1.2, dist="weibull", shape=2, scale=1)
#' X
#' # [1] 0.3984360 0.9623867 1.5793358 0.5810157 1.1990188 1.1023248 0.7779932
#' # [8] 1.7265598 1.6089237 1.8683241 1.9077715 1.1938959 2.1730424 1.8272045
#' # [15] 1.7397198 1.5702099 2.1379330 1.8103224 1.9107890 2.0386229 1.9813020
#' # [22] 2.6492839 2.3424281 2.4793168 2.4569918 2.0060300 2.2629934 2.3609622
#' # [29] 2.1392445 2.1252441 3.2127646 2.2451416 2.6506264 2.7209208 2.5084019
#' # [36] 2.8734849 3.1062846 2.5480120 2.6141440 2.5286307 2.5044793 2.8896440
#' # [43] 2.7151352 2.8182222 3.3822301 3.7379543 3.1460838 3.1213262 3.1302067
#' # [50] 3.1171138
#'
#' estimate_YNM_moments(X)
#' # $gamma
#' # [1] 1.282051
#'
#' # $variance
#' # [1] 0.4635951
#'
#' @export
estimate_YNM_moments <- function(X, variance = TRUE) {
  if (!is.numeric(X) || length(X) < 4)
    stop("X must be a numeric vector of length > 4.")

  # --- Compute number of records and sample size ---
  N <- rec_count(X)
  T <- length(X)

  # --- Gamma estimator ---
  gamma_hat <- 1 / (1 - N / T)

  # --- Analytical variance (optional) ---
  var_hat <- if (isTRUE(variance)) gamma_hat^2 * (gamma_hat - 1) else NA_real_

  # --- Return structured output ---
  return(list(
    param = gamma_hat,
    variance = var_hat
  ))
}

#########------------- Unbiased based on Moments ------------- ##########
#' Unbiased Estimation of Gamma Parameter based on number of records (with optional variance)
#'
#' Computes **bias-corrected** estimator of the (\eqn{\gamma}) parameter using a
#' bias correction, correcting
#' \code{\link{estimate_YNM_moments}}, and
#' optionally computes its analytical variance.
#'
#' @details
#' The unbiased estimator is given by:
#' \deqn{\hat{\gamma} = \hat{\gamma_{biased}} - \text{Bias}(\hat{\gamma_{biased}})}
#' The bias term is approximated as:
#' \deqn{
#' \text{Bias} =
#' \frac{h_T - P}{(1 - P)^2} +
#' \frac{1}{(1 - P)^3}
#' \left[
#'   \frac{h_T}{T}
#'   - \frac{1}{T^2}\sum_{t=1}^T \{ P_t^2 + (h_T - P)^2 \}
#' \right]
#' }
#' where:
#' \itemize{
#'   \item \eqn{P_t = P_t(\gamma)} - record rate at time finite \eqn{t}
#'   \code{\link{rec_rate_YNM}},
#'   \item \eqn{h_T = \frac{1}{T}\sum_{t=1}^T P_t},
#'   \item \eqn{P = \lim_{t\to\infty} P_t(\gamma)} - record rate at time infinite
#'   \eqn{t} \code{\link{rec_rate_YNM}}.
#' }
#'
#' ### Variance Approximation
#' The analytical variance is obtained via the **delta method**, as:
#' \deqn{
#' \text{Var}(\hat{\gamma}_u)
#'   \approx \lambda_\gamma
#'   \left( \frac{d H_T(\gamma)}{d \gamma} \right)^2,
#' }
#' where \eqn{\lambda_\gamma} is the variance of the biased estimator,
#' and the derivative term is expanded as:
#'
#' \deqn{
#' \frac{d H_T(\gamma)}{d \gamma}
#' = 1
#' - \frac{2\, dP(\gamma)/d\gamma}{(1-P(\gamma))^3}
#'   \left[
#'     \frac{1}{T}\sum_{t=1}^{T} P_t(\gamma) - P(\gamma)
#'   \right]
#' + \frac{1}{(1-P(\gamma))^2}
#'   \left[
#'     \frac{1}{T}\sum_{t=1}^{T} \frac{dP_t(\gamma)}{d\gamma} - \frac{dP(\gamma)}{d\gamma}
#'   \right]
#' }
#' \deqn{
#' + \frac{3\, dP(\gamma)/d\gamma}{(1-P(\gamma))^4}
#'   \left[
#'     \frac{1}{T^2}\sum_{t=1}^{T} P_t(\gamma)(1-P_t(\gamma))
#'     + \left(
#'       \frac{1}{T}\sum_{t=1}^{T} P_t(\gamma) - P(\gamma)
#'     \right)^2
#'   \right]
#' + \frac{1}{(1-P(\gamma))^3}
#'   \left[
#'     \frac{1}{T^2}\sum_{t=1}^{T} \frac{dP_t(\gamma)}{d\gamma}
#'     - \frac{1}{T^2}\sum_{t=1}^{T} \frac{d^2P_t(\gamma)}{d\gamma^2}
#'   \right]
#' + 2
#' \left(
#'   \frac{1}{T}\sum_{t=1}^{T} P_t(\gamma) - P(\gamma)
#' \right)
#' \left(
#'   \frac{1}{T}\sum_{t=1}^{T} \frac{dP_t(\gamma)}{d\gamma} - \frac{dP(\gamma)}{d\gamma}
#' \right).
#' }
#'
#'
#' @param X Numeric vector representing the observed sequence.
#' @param variance Logical; if TRUE (default), compute analytical variance.
#' @return A list with:
#' \item{gamma}{Bias-corrected (unbiased) gamma estimate.}
#' \item{variance}{Analytical variance of the unbiased estimator (if computed).}
#'
#' @export
#' @examples
#' set.seed(123)
#' X=YNM_series(T=50, gamma=1.2, dist="weibull", shape=2, scale=1)
#' X
#' # [1] 0.3984360 0.9623867 1.5793358 0.5810157 1.1990188 1.1023248 0.7779932
#' # [8] 1.7265598 1.6089237 1.8683241 1.9077715 1.1938959 2.1730424 1.8272045
#' # [15] 1.7397198 1.5702099 2.1379330 1.8103224 1.9107890 2.0386229 1.9813020
#' # [22] 2.6492839 2.3424281 2.4793168 2.4569918 2.0060300 2.2629934 2.3609622
#' # [29] 2.1392445 2.1252441 3.2127646 2.2451416 2.6506264 2.7209208 2.5084019
#' # [36] 2.8734849 3.1062846 2.5480120 2.6141440 2.5286307 2.5044793 2.8896440
#' # [43] 2.7151352 2.8182222 3.3822301 3.7379543 3.1460838 3.1213262 3.1302067
#' # [50] 3.1171138
#'
#' estimate_YNM_moments_unbias(X)
#' # $gamma
#' # [1] 1.212979
#' #
#' # $variance
#' # [1] 0.2997593
estimate_YNM_moments_unbias = function(X, variance = TRUE){ ## compute the second estimator*

  if (!is.numeric(X) || length(X) < 4) stop("X must be a numeric vector of length >= 4.")

  T <- length(X)

  # --- 0. Obtain baseline NT estimate (assumes this function exists and returns list(param, variance)) ---
  # Use variance = FALSE to avoid unnecessary computation inside that function if we will compute variance here
  if (!exists("estimate_YNM_moments", mode = "function")) {
    stop("Required helper function 'estimate_YNM_moments' not found in the environment.")
  }
  Estimated <- estimate_YNM_moments(X = X, variance = TRUE)
  if (!is.list(Estimated) || is.null(Estimated$param)) {
    stop("estimate_YNM_moments must return a list with at least element $gamma.")
  }

  gamma_biased = Estimated$param

  hT = mean(rec_rate_YNM(gamma_biased,1:T)) ## sum /T

  P = rec_rate_YNM(gamma_biased, t=Inf)

  ### Bias estimator function
    term1 = (hT-P)/(1-P)^2

    term2= 1/(1-P)^3

    term3 = rec_rate_YNM(gamma_biased,t=seq_along(T))^2 + (hT-P)^2

    term4=hT/T - (sum(term3))/T^2

    bias = term1+term2*term4


  gamma = gamma_biased - bias

  # --- 2. Optional: compute variance of the unbiased estimator ---
  var_out <- NA_real_
  if (isTRUE(variance)) {
    ## Record rate
    P = rec_rate_YNM(gamma, t=Inf)
    P_t = rec_rate_YNM(gamma, t=1:T)  ## vectorized

    ## First order derivative
    dP_t_f = function(gamma, t) (gamma^(t - 2) * (gamma^t - t * gamma + t - 1)) / (gamma^t - 1)^2
    dP_t = dP_t_f(gamma=gamma, t= 1:T )  ##vectorized
    dP =  1/gamma^2

    ## means
    mean_P_t = mean(P_t)
    mean_dP_t = mean(dP_t)
    term = mean_P_t - P

    fisher = 1
    fisher = fisher - (2*dP *( term)) /(1-P)^3
    fisher = fisher + (dP *( mean_dP_t - dP)) /(1-P)^2
    fisher = fisher + (3*dP *(sum(P_t*(1-P_t))/T^2) * (term)^2 ) /(1-P)^4
    fisher = fisher + (mean_dP_t/T - sum(dP_t^2)/T^2 + 2*(term) * (mean_dP_t - dP) ) /(1-P)^3

    # lambda_gamma: variance of the biased NT-estimator (obtained from estimate_LDM_moments if available)
        lambda_gamma <- NA_real_
    if (!is.null(Estimated$variance) && is.numeric(Estimated$variance)) {
      lambda_gamma <- Estimated$variance
    } else {
      warning("Could not locate variance for the base NT estimator. Returning NA for final variance.")
      lambda_gamma <- NA_real_
    }

    if (!is.na(lambda_gamma)) {
      var_out <- as.numeric(lambda_gamma * (fisher)^2)
    } else {
      var_out <- NA_real_
    }
  }

  return(list(param = gamma, variance = var_out))
}


#########------------- MLE indicator ------------------------- ##########
#' Maximum Likelihood Estimation of Gamma Based on Indicator Series
#'
#' Estimates the gamma (\eqn{\gamma}) parameter using the maximum likelihood
#' method based on the record indicator series, and optionally computes the
#' analytical or approximate variance.
#'
#' @details
#' The log-likelihood function is given by:
#' \deqn{
#' \log \ell(v = 1/\gamma) =
#' -N \log(1 - v)
#' + (T - N)\log(v)
#' - \log(1 - v^T)
#' - \sum_{t=2}^{T} \delta_t \log(1 - v^{t-1})
#' }
#' where:
#' \itemize{
#'   \item \eqn{v = 1 / \gamma},
#'   \item \eqn{N} is the number of records,
#'   \item \eqn{T} is the sequence length,
#'   \item \eqn{\delta_t} are the record indicator variables.
#' }
#'
#' The estimate \eqn{\hat{\gamma}} is found by setting the derivative
#' of the log-likelihood with respect to \eqn{v} to zero and solving numerically.
#'
#'
#' The **exact Fisher information** for \eqn{\gamma} is:
#' \deqn{
#' I_T(\gamma) =
#' \frac{1}{\gamma^2 (\gamma - 1)^2} \sum_{t=1}^{T} P_t(\gamma)
#' + \frac{1}{\gamma^2} (T - \sum_{t=1}^{T} P_t(\gamma))
#' - \frac{T (1 + \gamma^T (T - 1))}{\gamma^2 (\gamma^T - 1)^2}
#' - \sum_{t=2}^{T} \frac{(t-1)(1 + (t-2)\gamma^{t-1}) P_t(\gamma)}{\gamma^2 (\gamma^{t-1} -1)^2}
#' }
#'
#' The asymptotic approximation for large \eqn{T} is:
#' \deqn{
#' \text{Var}(\hat{\gamma}) \approx \frac{1 - v}{T v^3}, \quad v = 1/\gamma
#' }
#'
#' @param X Numeric vector representing the sequence.
#' @param min Lower bound for the search range of \eqn{\gamma}. Default is 1.
#' @param max Upper bound for the search range of \eqn{\gamma}. Default is 5.
#' @param step Step size for the search grid. Default is 0.001.
#' @param variance Logical; if TRUE (default), compute the analytical or approximate variance.
#' @param approximate Logical; if TRUE, use the asymptotic approximation for the variance.
#'
#' @return A list containing:
#' \item{gamma}{Estimated \eqn{\gamma} parameter.}
#' \item{variance}{Estimated variance (exact or approximate).}
#'
#' @export
#' @examples
#' set.seed(123)
#' X=YNM_series(T=50, gamma=1.2, dist="weibull", shape=2, scale=1)
#' X
#' # [1] 0.3984360 0.9623867 1.5793358 0.5810157 1.1990188 1.1023248 0.7779932
#' # [8] 1.7265598 1.6089237 1.8683241 1.9077715 1.1938959 2.1730424 1.8272045
#' # [15] 1.7397198 1.5702099 2.1379330 1.8103224 1.9107890 2.0386229 1.9813020
#' # [22] 2.6492839 2.3424281 2.4793168 2.4569918 2.0060300 2.2629934 2.3609622
#' # [29] 2.1392445 2.1252441 3.2127646 2.2451416 2.6506264 2.7209208 2.5084019
#' # [35] 2.8734849 3.1062846 2.5480120 2.6141440 2.5286307 2.5044793 2.8896440
#' # [43] 2.7151352 2.8182222 3.3822301 3.7379543 3.1460838 3.1213262 3.1302067
#' # [50] 3.1171138
#'
#' estimate_YNM_mle_indicator(X)
#' # $gamma
#' # [1] 1.183
#'
#' # $variance
#' # [1] 0.006915017
#'
#' estimate_YNM_mle_indicator(X, approximate = TRUE)
#' # $gamma
#' # [1] 1.183
#'
#' # $variance
#' # [1] 0.00512213
#' @export
estimate_YNM_mle_indicator <- function(X, variance = TRUE, approximate = FALSE,
                                       min = 1, max = 5, step = 0.001) {

  if (!is.numeric(X) || length(X) < 4)
    stop("X must be a numeric vector of length > 4.")

  # --- Precompute quantities ---
  T <- length(X)
  delta <- is_rec(X)
  N <- rec_count(X)

  # --- Derivative of log-likelihood wrt v ---
  dLog <- function(v) {
    k <- 2:T
    term <- delta[k] * (k - 1) * v^(k - 2) / (1 - v^(k - 1))
    -N / (1 - v) + (T - N) / v + (T * v^(T - 1)) / (1 - v^T) + sum(term)
  }

  # --- Search for optimal gamma ---
  gamma_grid <- seq(min, max, by = step)
  v_grid <- 1 / gamma_grid
  dvals <- vapply(v_grid, dLog, numeric(1))
  gamma_hat <- gamma_grid[which.min(abs(dvals))]

  # --- Variance estimation ---
  var_hat <- NA_real_

  if (isTRUE(variance)) {
    if (approximate) {
      # Approximate variance (large-sample)
      v <- 1 / gamma_hat
      var_hat <- (1 - v) / (T * v^3)
    } else {
      # Exact Fisher Information (requires rec_count_mean_YNM and rec_rate_YNM)
      ent <- rec_count_mean_YNM(T = T, gamma = gamma_hat)
      a <- (1 / (gamma_hat^2 * (gamma_hat - 1)^2)) * ent
      b <- (1 / gamma_hat^2) * (T - ent)
      c <- T * (1 + gamma_hat^T * (T - 1)) / (gamma_hat^2 * (gamma_hat^T - 1)^2)

      i <- 2:T
      d <- (i - 1) * (1 + (i - 2) * gamma_hat^(i - 1)) *
        rec_rate_YNM(gamma_hat, i) /
        (gamma_hat^2 * (gamma_hat^(i - 1) - 1)^2)

      I <- a + b - c - sum(d)
      var_hat <- 1 / I
    }
  }

  # --- Return results ---
  return(list(
    param = gamma_hat,
    variance = var_hat
  ))
}

