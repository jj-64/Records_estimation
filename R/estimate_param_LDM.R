
## Based on moments, number of records ---------------------
#' Estimation of \eqn{\theta} in LDM Process Based on Number of Records (with optional variance)
#'
#' Estimates the parameter \eqn{\theta} in the Linear Drift Model (LDM) process using
#' the number of records, if requested, it also computes the variance of the estimator.
#'
#' @details
#' An LDM process is defined as:
#' \deqn{ X_t = Y_t + \theta t, }
#' where \eqn{Y_t} are independent and identically distributed (i.i.d.) random variables.
#' The properties are studied when \eqn{Y_t} follows a \eqn{Gumbel(\alpha, \beta)} distribution.
#'
#' The estimation formula is given by:
#' \deqn{ \hat{\theta} = -\log\left(1 - \frac{N_T}{T} \right) }
#' where \eqn{N_T} is the number of records observed in the process \eqn{X_t}, and \eqn{T} is the length of the process.
#'
#'
#' The variance formula is given by:
#' \deqn{ V(\hat{\theta}) = \frac{1 - e^{-\theta / \text{scale}}}{e^{-\theta / \text{scale}}} }
#'
#' This estimator is known to be biased, and a bias-corrected version is provided in
#' \code{\link{estimate_LDM_moments_unbias}}.
#' @param X numeric vector representing the time series data.
#' @param variance Logical if to compute variance (default = TRUE)
#' @param scale positive Numeric. Default is 1. The scale parameter of the \eqn{Gumbel} underlying distribution used for variance computation.
#' @return A list with elements:
#'   \item{theta}{Estimate of \eqn{\theta}.}
#'   \item{variance}{Asymptotic variance (or NA if variance = FALSE or could not be computed).}
#' @export
#' @examples
#' Yt <- rnorm(25)
#' Xt <- Yt + 0.2 * (1:25)
#' rec_count(Xt)  # Number of records
#' estimate_LDM_moments(X = Xt)
#' estimate_LDM_moments(X = c(0.428,1.311,2.023,2.882,2.096,-0.197,1.339,
#' 1.748,1.418, 0.711, 1.999,3.598, 3.308, 3.942,2.025,3.282,4.043, 0.492,
#' 4.639, 1.408, 3.525, 5.398,  3.719, 3.741, 4.729))
#' # $theta
#' # [1] 0.4462871
#'
#' # $variance
#' # [1] 0.5625
estimate_LDM_moments <- function(X, variance = TRUE, scale=1) {
  if (!is.numeric(X) || length(X) < 4) stop("X must be a numeric vector of length >= 4.")
  # stopifnot(is.numeric(X), length(X) >= 4)
  # stopifnot(scale > 0)
  if (scale <= 0)
    stop("scale must be positive.")

	# estimator
    theta = -log(1 - rec_count(X) / length(X))
    var_theta = NA
	 # variance formula
  if(variance) {var_theta =  (1-exp(-theta/scale))/exp(-theta/scale) }
  return(list("param"= theta, "variance"=var_theta))
}


## Unbiased Moments -----------------------------

#' Bias-corrected estimator for theta in LDM (with optional variance)
#'
#' Compute a bias-corrected estimate of \eqn{\theta} for the Linear Drift Model (LDM)
#' based on the number of records (NT) estimator in \code{\link{estimate_LDM_moments}},
#' and optionally compute the
#' asymptotic variance of the bias-corrected estimator.
#' If requested (variance = TRUE), it computes the variance \eqn{\vartheta(\theta)}
#' for the unbiased estimator of \eqn{\theta} in the LDM process.
#'
#' The function calls the (unbiased) NT estimator routine "estimate_LDM_moments"
#' (which must return a list with \eqn{\theta} and optionally "variance"), applies the
#' analytical bias correction and — if requested — computes the variance using the
#' first order derivative formula and the variance of the original estimator.
#'
#' @details
#' The correction is computed based on the expected bias of the estimator.
#'
#' The variance formula is given by:
#'
#' \deqn{\vartheta(\theta) = \lambda(\theta) \left( \frac{d H_T(\theta)}{d \theta} \right)^2}
#'
#' where:
#'
#' \deqn{\frac{d H_T(\theta)}{d \theta} = 1 - \frac{d G_T(\theta)}{d \theta}}
#'
#' \deqn{= 1 - \frac{\frac{d P(\theta)}{d \theta}}{(1 - P(\theta))^2} \left[ \frac{1}{T} \sum_{t=1}^{T} P_t(\theta) - P(\theta) \right]}
#'
#' \deqn{+ \frac{1}{(1 - P(\theta))} \left[ \frac{1}{T} \sum_{t=1}^{T} \frac{d P_t(\theta)}{d \theta} - \frac{d P(\theta)}{d \theta} \right]}
#'
#' \deqn{+ \frac{\frac{d P(\theta)}{d \theta}}{(1 - P(\theta))^3} \left[ \frac{1}{T^2} \sum_{t=1}^{T} P_t(\theta) (1 - P_t(\theta)) + \left( \frac{1}{T} \sum_{t=1}^{T} P_t(\theta) - P(\theta) \right)^2 \right]}
#'
#' \deqn{+ \frac{1}{2 (1 - P(\theta))^2} \left[ \frac{1}{T^2} \sum_{t=1}^{T} \frac{d P_t(\theta)}{d \theta} - \frac{1}{T^2} \sum_{t=1}^{T} \frac{d^2 P_t(\theta)}{d \theta^2} \right]}
#'
#' \deqn{+ 2 \left( \frac{1}{T} \sum_{t=1}^{T} P_t(\theta) - P(\theta) \right) \left( \frac{1}{T} \sum_{t=1}^{T} \frac{d P_t(\theta)}{d \theta} - \frac{d P(\theta)}{d \theta} \right)}
#'
#'where \eqn{P(\theta) = 1-e^{-\theta}}, \eqn{P_t(\theta) = \frac{1-e^{-\theta}}{1-e^{-\theta t}}}, and \eqn{\lambda(\theta)} is the variance of the estimated paramter with bias in \code{\link{estimate_LDM_moments}}
#'
#' @param X Numeric vector. Observed series.
#' @param variance Logical. If TRUE (default) compute and return variance;
#' otherwise skip heavy computation.
#' @param scale Numeric. Scale parameter for the LDM (defualt = 1. theta is assumed to be in units such that
#'              the user may pass scaled theta; internally we treat theta/scale).
#' @return A list with elements:
#'   \item{theta}{Bias-corrected estimate of \eqn{\theta}.}
#'   \item{variance}{Asymptotic variance (or NA if variance = FALSE or could not be computed).}
#' @export
#' @examples
#' \dontrun{
#' Yt <- rnorm(25)
#' Xt <- Yt + 0.2 * (1:25)
#' rec_count(Xt)  # Number of records
#' res <- estimate_LDM_moments_unbias(Xt, variance = TRUE, scale = 1)
#' res$theta; res$variance
#' }
#' estimate_LDM_moments_unbias(X = c(0.428,1.311,2.023,2.882,2.096,-0.197,1.339,
#' 1.748,1.418, 0.711, 1.999,3.598, 3.308, 3.942,2.025,3.282,4.043, 0.492, 4.639,
#'  1.408, 3.525, 5.398,  3.719, 3.741, 4.729))
#'
#' # $theta
#' # [1] 0.3601306
#'
#' # $variance
#' # [1] 0.3116952
estimate_LDM_moments_unbias <- function(X, variance = TRUE, scale = 1) {

  if (!is.numeric(X) || length(X) < 4) stop("X must be a numeric vector of length >= 4.")
  if (scale <= 0)
    stop("scale must be positive.")

  T <- length(X)

  # --- 0. Obtain baseline NT estimate (assumes this function exists and returns list(theta, variance)) ---
  # Use variance = FALSE to avoid unnecessary computation inside that function if we will compute variance here
  if (!exists("estimate_LDM_moments", mode = "function")) {
    stop("Required helper function 'estimate_LDM_moments' not found in the environment.")
  }

  	## Obtain biased estimator
  Estimated <- estimate_LDM_moments(X = X, variance = TRUE, scale = scale)
  if (!is.list(Estimated) || is.null(Estimated$param)) {
    stop("estimate_LDM_moments must return a list with at least element $param.")
  }

  theta_biased <- Estimated$param
  # theta in this function is treated as theta_scaled = theta / scale (user said theta is theta/scale)
  # but to be explicit: treat input/returned theta on same scale as estimate_LDM_moments
  theta_scaled_biased <- theta_biased / scale

  # --- 1. Bias correction (vectorized) ---
  # Precompute P_t vector using your existing rec_rate_LDM implementation
  t_seq <- seq_len(T)
  P_t_vec <- rec_rate_LDM(t = t_seq, theta = theta_biased, scale = scale) # vector length T
  mean_P_t <- mean(P_t_vec)  # 1/T * sum P_t

  P_theta <- 1 - exp(-theta_scaled_biased)    # scalar
  # a and b as in your bias formula
  term1 <- (mean_P_t - P_theta) / (1 - P_theta)
  term2 <- 1 / (2 * (1 - P_theta)^2)

  # cc vector: rec_rate_LDM(theta,t)^2 + (mean_P_t - P_theta)^2
  term3 <- P_t_vec^2 + (mean_P_t - P_theta)^2
  term4<- mean_P_t / T - sum(term3) / (T^2)

  bias_est <- term1 + term2 * term4

  # bias corrected theta (on original scale returned by estimate_LDM_moments)
  theta_unbiased <- theta_biased - bias_est

  # If these theta values are supposed to be returned scaled or unscaled depends on estimate_LDM_moments contract
  # We return on the same scale as estimate_LDM_moments returned theta (consistent).

  # --- 2. Optional: compute variance of the unbiased estimator ---
  var_out <- NA_real_
  if (isTRUE(variance)) {
    # Use theta_scaled (theta / scale) for analytical derivatives
    theta_scaled <- theta_unbiased / scale

    # Basic scalars
    P_theta <- 1 - exp(-theta_scaled)
    dP_dtheta <- exp(-theta_scaled)     # d/dtheta P(theta)

    # Precompute commonly used arrays
    t_seq <- seq_len(T)
    exp_t_theta <- exp(t_seq * theta_scaled)
    denom <- exp_t_theta - 1
    # Avoid division by zero issues for very small denom
    if (any(denom == 0)) denom[denom == 0] <- .Machine$double.xmin

    # First derivative dP_t/dtheta (vectorized)
    # expression simplified from algebraic expression you provided
    dP_t <- ((exp_t_theta - t_seq * exp(theta_scaled) + t_seq - 1) *
               exp_t_theta * exp(-theta_scaled)) / (denom^2)

    # Second derivative d^2 P_t / dtheta^2 (vectorized)
    d2P_t <- -((exp(2 * t_seq * theta_scaled) +
                  (-t_seq^2 * exp(theta_scaled) + t_seq^2 + 2 * t_seq - 2) * exp_t_theta -
                  t_seq^2 * exp(theta_scaled) + t_seq^2 - 2 * t_seq + 1) *
                 exp_t_theta * exp(-theta_scaled)) / (denom^3)

    # Means: note these are 1/T * sum(...)
    mean_P_t <- mean(P_t_vec)           # 1/T * sum P_t
    mean_dP_t <- mean(dP_t)             # 1/T * sum dP_t
    mean_d2P_t <- mean(d2P_t)           # 1/T * sum d2P_t

    # dH/dθ (vectorized)
    # Compute derivative dH_T/dtheta using the formula (vectorized pieces)
    dH_T_dtheta <- 1 -
      (dP_dtheta / (1 - P_theta)^2) * (mean_P_t - P_theta) +
      (1 / (1 - P_theta)) * (mean_dP_t - dP_dtheta) +
      (dP_dtheta / (1 - P_theta)^3) * ((mean(P_t_vec * (1 - P_t_vec)) / T) + (mean_P_t - P_theta)^2) +
      (1 / (2 * (1 - P_theta)^2)) * (mean_dP_t / T - mean_d2P_t / T) +
      2 * (mean_P_t - P_theta) * (mean_dP_t - dP_dtheta)

	 # lambda = variance of biased NT estimator
    # lambda_theta: variance of the biased NT-estimator (obtained from estimate_LDM_moments if available)
    # If estimate_LDM_moments returned variance use it; otherwise fallback to analytical function if available.
    lambda_theta <- NA_real_
    if (!is.null(Estimated$variance) && is.numeric(Estimated$variance)) {
      lambda_theta <- Estimated$variance
    }  else {
      warning("Could not locate variance for the base NT estimator. Returning NA for final variance.")
      lambda_theta <- NA_real_
    }

    if (!is.na(lambda_theta)) {
      var_out <- as.numeric(lambda_theta * (dH_T_dtheta)^2)
    } else {
      var_out <- NA_real_
    }
  }

  return(list(param = theta_unbiased, variance = var_out))
}

## MLE -------------------------------------


#' Maximum Likelihood Estimation of θ Based on record indicator Series (with optional variance)
#'
#' Computes the Maximum Likelihood Estimate (MLE) of the parameter θ
#' for the Linear Drift Model (LDM), based on the record indicator series \eqn{\delta_t}.
#'
#' Optionally computes its analytical variance.
#'
#' @param X Numeric vector — the observed time series.
#' @param variance Logical. If TRUE, compute variance of θ. Default = TRUE.
#' @param scale Numeric. Scale parameter (default = 1).
#' @param min Lower bound for θ grid search space. Default = 0.0001.
#' @param max Upper bound for θ grid search space. Default = 5.
#' @param step Search grid step size (default = 0.01).
#'
#' @details
#' The log-likelihood is:
#' \deqn{
#' \log L(\theta) = N \log(1 - z) + (T - N) \log z - \log(1 - z^T)
#' - \sum_{t=2}^{T} \delta_t \log(1 - z^{t-1})
#' }
#' where \eqn{z = e^{-\theta}}, \eqn{N} is the number of records, \eqn{T} the series length,
#' and \eqn{\delta_t} are the record indicator values.
#'
#'
#'
#'
#' @return A list with:
#' \item{theta}{The MLE of θ.}
#' \item{variance}{Analytical variance (or NA if variance = FALSE).}
#'
#' @export
#' @examples
#' Yt <- rnorm(25)
#' Xt <- Yt + 0.2 * (1:25)
#' rec_count(Xt)  # Number of records
#' res <- estimate_LDM_mle_indicator(Xt, variance = TRUE, scale = 1)
#' res$theta; res$variance
#'
#' estimate_LDM_mle_indicator (X = c(0.428,1.311,2.023,2.882,2.096,-0.197,
#' 1.339,  1.748,1.418, 0.711, 1.999,3.598, 3.308, 3.942,2.025,3.282,4.043,
#' 0.492, 4.639, 1.408, 3.525, 5.398,  3.719, 3.741, 4.729))
#'
#' # $theta
#' # [1] 0.3301
#'
#' # $variance
#' # [1] 0.02139224
estimate_LDM_mle_indicator <- function(X, variance = TRUE, scale = 1, min = 0.0001,
                                       max = 5, step = 0.01) {
  if (!is.numeric(X) || length(X) < 4)
    stop("X must be a numeric vector with length > 4.")
  if (scale <= 0)
    stop("scale must be positive.")

  # --- Basic quantities ---
  T <- length(X)
  delta <- is_rec(X)        # Record indicator sequence
  N <- rec_count(X)        # Total number of records

  # --- Grid search for θ (vectorized, fast) ---
  theta_grid <- seq(min, max, by = step)
  z_grid <- exp(-theta_grid)

  # Vectorized log-likelihood computation for each z
  logL_values <- sapply(z_grid, function(z) {
    term1 <- N * log(1 - z)
    term2 <- (T - N) * log(z)
    term3 <- -log(1 - z^T)
    term4 <- -sum(delta[-1] * log(1 - z^(1:(T - 1))))
    return(term1 + term2 + term3 + term4)
  })

  # Find MLE
  theta_hat <- theta_grid[which.max(logL_values)]

  # --- Compute variance if requested ---
  var_out <- NA_real_
  if (isTRUE(variance)) {
    # Work in scaled parameterization
    theta_scaled <- theta_hat / scale

    # --- Terms for variance formula ---
    a <- exp(-2 * theta_scaled) / (1 - exp(-theta_scaled))^2
    b <- sum(rec_rate_LDM(theta = theta_scaled, t = 1:T, scale = 1))

    cc <- T * exp(-T * theta_scaled) * (T + exp(-T * theta_scaled) - 1) /
      (1 - exp(-T * theta_scaled))^2

    # Function d(t, θ)
    d_t <- (2:T - 1) * exp(-theta_scaled * (2:T - 1)) *
      ((2:T - 2) + exp(-theta_scaled * (2:T - 1))) *
      rec_rate_LDM(theta = theta_scaled, t = 2:T, scale = 1) /
      (1 - exp(-theta_scaled * (2:T - 1)))^2

    dd <- sum(d_t)

    var_out <- (a * b + T - b - cc - dd)^(-1)
  }

  return(list(param = theta_hat, variance = var_out))
}

