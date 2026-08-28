
############## Perform DTRW test ##############
#' Test for Discrete-Time Random Walk (DTRW) Structure
#'
#' Performs a composite hypothesis test to assess whether a given time series
#' follows a Discrete-Time Random Walk (DTRW) process. The test combines three
#' component tests — stationarity, independence, and symmetry — using
#' multiple-comparison correction methods (Bonferroni, Holm-Bonferoni, Holm-Sidák, or chi-square).
#'
#' @details
#' The function applies the following tests:
#' \enumerate{
#'   \item \strong{Stationarity:} Augmented Dickey–Fuller (ADF) test
#'         using \code{tseries::adf.test}.
#'         Null hypothesis (\eqn{H_0}): series is non-stationary.
#'         Desired outcome: fail to reject \eqn{H_0} (\eqn{p > \alpha}).
#'
#'   \item \strong{Independence:} Ljung–Box test (\code{Box.test}) on the increments.
#'         Null hypothesis (\eqn{H_0}): increments are independently distributed.
#'         Desired outcome: fail to reject \eqn{H_0}.
#'
#'   \item \strong{Symmetry:} Wilcoxon signed-rank test (\code{wilcox.test}) for zero median increments.
#'         Null hypothesis (\eqn{H_0}): increments are symmetric about zero.
#'         Desired outcome: fail to reject \eqn{H_0}.
#' }
#'
#' The overall decision is based on combining p-values using one of the following methods:
#' \itemize{
#'   \item \code{"Bonf"} — Bonferroni correction (default)
#'   \item \code{"Holm"} — Holm–Bonferroni sequential correction
#'   \item \code{"Sidak"} — Holm-Sidák correction
#'   \item \code{"Chisq"} — Fisher’s chi-squared combination test
#' }
#'
#' The decision rule is "DTRW" if all component tests fail to reject their nulls
#' under the chosen adjustment method; otherwise "NO".
#'
#' @param X Numeric vector of observations representing the process path.
#' @param alpha Numeric, significance level (default = 0.05).
#' @param method Character, p-value combination method: one of
#'   \code{"Bonf"}, \code{"Holm"}, \code{"Sidak"}, or \code{"Chisq"} (default = "Bonf").
#'
#' @return A list with the following elements:
#' \item{method}{Method used to combine tests.}
#' \item{p_valueStationary}{P-value from ADF stationarity test.}
#' \item{p_valueIndep}{P-value from Ljung–Box independence test.}
#' \item{p_valueSymm}{P-value from Wilcoxon symmetry test.}
#' \item{decision}{Overall decision: "DTRW" if consistent with random walk assumptions, "NO" otherwise.}
#'
#' @references
#' Ljung, G. M. and Box, G. E. P. (1978). “On a Measure of Lack of Fit in Time Series Models.”
#' \emph{Biometrika}, 65(2), 297–303.
#'
#' Holm, S. (1979). “A Simple Sequentially Rejective Multiple Test Procedure.”
#' \emph{Scandinavian Journal of Statistics}, 6(2), 65–70.
#'
#' Fisher, R. A. (1932). \emph{Statistical Methods for Research Workers.}
#' @export
#' @examples
#' set.seed(123)
#' X <- cumsum(rnorm(100))  # Simulated random walk
#' Test_DTRW_Indep(X, alpha = 0.05, method = "Bonf")
#'
Test_DTRW_Indep = function(X,alpha=0.05, method="Bonf"){

  X=X-X[1]
  increments <- diff(X)
  ## Test1: Perform ADF test H1: stationary
    #adf_test <- tseries::adf.test(X)
    #p_value1 = adf_test$p.value ## we want H0: series is not stationary, so fail to reject H0, so p_value>alpha

  ## Test2:  Perform Ljung-Box test on increments: H0:independently distributed (no autocorrelation)
    ljung_box_test <- Box.test(increments, type = "Ljung-Box")
    p_value2 = ljung_box_test$p.value  ## we want H0: fail to reject so p_value>p

  ##Test3: Wilcoxon signed-rank test for symmetry H0: symmetric
    wilcoxon_test = wilcox.test(increments, mu =0, alternative = "two.sided", exact = FALSE) #median(increments) # DTRW without drift
    p_value3 = wilcoxon_test$p.value


  ## Decision logic
    if (method == "Bonf") {
      dec <- ifelse(p_value2 > alpha / 2 & p_value3 > alpha / 2, "DTRW", "NO")

    } else if (method == "Holm") {
      pp <- sort(c(p_value2, p_value3))
      Holm <- pp > c(1 - (1 - alpha)^(1/2), 1 - (1 - alpha)^(1/2), alpha)
      dec <- ifelse(sum(Holm) == 2, "DTRW", "NO") ## Holm-Bonf is true true true

    } else if (method == "Sidak") {
      pp <- c( p_value2, p_value3)
      Sidak <- pp > (1 - (1 - alpha)^(1/2))
      dec <- ifelse(sum(Sidak) == 2, "DTRW", "NO") ## Holm-Sidak is true true true

    } else if (method == "Chisq") {
      stat <- -2 * (log(p_value1) + log(p_value2) + log(p_value3))
      dec <- ifelse(stat <= qchisq(1 - alpha, df = 6), "DTRW", "NO") #high typr I error

    } else {
      stop("method must be one of 'Bonf', 'Holm', 'Sidak', or 'Chisq'")
    }

    return(list(
      method = method,
      aleternativeIndep = "correlated",
      p_valueIndep = p_value2,
      aleternativeSymm = "NotSymmetric",
      p_valueSymm = p_value3,
      decision = dec
    ))
}


Test_DTRW_Indep2 = function(X,alpha=0.05){

  ## Test2:  Perform Ljung-Box test on increments: H0:independently distributed (no autocorrelation)
  ljung_box_test <- Box.test(diff(X), type = "Ljung-Box")
  p_value2 = ljung_box_test$p.value  ## we want H0: fail to reject so p_value>p

  ##Test3: Wilcoxon signed-rank test for symmetry H0: symmetric
  wilcoxon_test = wilcox.test(diff(X), mu =0, alternative = "two.sided", exact = FALSE) #median(increments) # DTRW without drift
  p_value3 = wilcoxon_test$p.value

  dec = ifelse(p_value2 < alpha/2 | p_value3< alpha/2 , "NO", "DTRW")
  return(list(
    aleternative = "correlated",
    p_value = p_value2,
    decision = dec
  ))
}
############## based on exact distribution of number of records ##############
#' Compute Quantile Bounds for Record Counts under a DTRW Process
#'
#' Computes lower and upper quantile bounds for the distribution of the number
#' of records in a Discrete-Time Random Walk (DTRW) of length \code{T}, using
#' theoretical probabilities from the function \code{rec_count_dist_DTRW}.
#'
#'
#' @details
#' The function constructs the empirical cumulative distribution function (CDF)
#' of the number of records observed in a DTRW and determines the quantiles
#' corresponding to the lower and upper tail probabilities \eqn{\alpha} and
#' \eqn{1 - \alpha}.
#'
#' Let \eqn{N_T} denote the number of records observed in a DTRW process of
#' length \eqn{T}. The cumulative distribution function (CDF) is computed as:
#' \deqn{
#'   F(n) = \sum_{i=1}^{n} P(N_T = i)
#' }
#' where \eqn{P(N_T = i)} is obtained from \code{rec_count_dist_DTRW(m = i, T = T)}.
#'
#' #' The function then finds the quantile indices corresponding to the
#' lower and upper tail probabilities \eqn{\alpha} and \eqn{1 - \alpha}:
#' \deqn{
#'   Q_{\alpha} = \min_n |F(n) - \alpha|, \quad
#'   Q_{1-\alpha} = \min_n |F(n) - (1 - \alpha)|
#' }
#'
#' These quantiles can be used as critical values (two-sided acceptance region) in the record-based test
#' \code{Test_DTRW_NT()}.
#'
#' @param alpha Numeric, significance level (default = 0.05).
#' @param T Integer, length of the random walk (number of time steps).
#'
#' @return A numeric vector of length 2 containing:
#' \item{1}{Lower quantile index corresponding to \eqn{\alpha}.}
#' \item{2}{Upper quantile index corresponding to \eqn{1 - \alpha}.}
#'
#' @seealso \code{\link{Test_DTRW_NT}}, \code{\link{rec_count_dist_DTRW}}
#'
#' @examples
#' Quantile_DTRW(alpha = 0.05, T = 50)
#'
#' @export
Quantile_DTRW <- function(alpha = 0.05, T) {
  Prob <- numeric(T)
  for (i in 1:T) {
    Prob[i] <- rec_count_dist_DTRW(m = i, T = T)
  }
  CDF <- cumsum(Prob)
  c(
    which.min(abs(CDF - alpha)),
    which.min(abs(CDF - (1 - alpha)))
  )
}



#' Test for Discrete-Time Random Walk (DTRW) Using Record Counts
#'
#' Performs a hypothesis test to determine whether a given time series
#' is consistent with the record statistics expected under a
#' Discrete-Time Random Walk (DTRW) process.
#'
#' @details
#' The test compares the observed number of records \eqn{N_T} in the series
#' \eqn{X_1, X_2, \ldots, X_T}. (obtained via \code{rec_count()}) against theoretical quantiles or
#' normal approximations of the record count distribution.
#'
#'For very large T, under the DTRW null hypothesis, the expected number of records and its variance
#' have asymptotic forms:
#' \deqn{
#'   E[N_T] \approx c_1 \sqrt{T}, \qquad Var[N_T] \approx c_2 T,
#' }
#' where constants \eqn{c_1} and \eqn{c_2} are model-specific.
#' A standardized test statistic is computed as:
#' \deqn{
#'   Z_T = \frac{N_T - E[N_T]}{\sqrt{Var[N_T]}}.
#' }
#'
#' Two testing modes are available:
#'
#' \itemize{
#'   \item \strong{Exact test:} (\code{approximate = FALSE}) Uses empirical quantiles from \code{Quantile_DTRW()}.
#'   \item \strong{Approximate test:} (\code{approximate = TRUE}) Uses asymptotic normal approximation
#'         with mean \code{rec_count_mean_DTRW(T)} and variance \code{rec_count_var_DTRW(T)}.
#' }
#'
#' Decision rule:
#' \itemize{
#'   \item For two-sided tests (\code{one.sided = FALSE}), accept DTRW if
#'         the observed statistic lies between the lower and upper quantiles.
#'   \item For one-sided tests (\code{one.sided = TRUE}), accept if below the upper quantile.
#' }
#'
#' @param X Numeric vector representing the observed series.
#' @param alpha Numeric, significance level (default = 0.05).
#' @param approximate Logical, if \code{TRUE} use the asymptotic normal approximation
#'   (default = \code{FALSE} for the exact quantile test).
#' @param one.sided Logical, if \code{TRUE} perform a one-sided test
#'   (default = \code{FALSE} for two-sided).
#'
#' @return A list containing:
#' \item{stat}{Observed test statistic (number of records or standardized z-score).}
#' \item{stat_theo}{Theoretical quantile bounds (if \code{approximate = FALSE}).}
#' \item{p_value}{Approximate p-value (if \code{approximate = TRUE}).}
#' \item{decision}{Character, "DTRW" if consistent with DTRW hypothesis, otherwise "NO".}
#'
#' @seealso \code{\link{Quantile_DTRW}}, \code{\link{rec_count}},
#'   \code{\link{rec_count_dist_DTRW}}, \code{\link{rec_count_mean_DTRW}}
#'
#' @examples
#' set.seed(123)
#' X <- cumsum(rnorm(100))
#'
#' # Exact quantile-based test
#' Test_DTRW_NT(X, alpha = 0.05, approximate = FALSE)
#'
#' # Asymptotic normal approximation
#' Test_DTRW_NT(X, alpha = 0.05, approximate = TRUE)
#'
#' @export
Test_DTRW_NT <- function(X, alpha = 0.05, approximate = FALSE, one.sided = FALSE) {

  obs <- rec_count(X)

  if (!approximate) {

    z_theo <- Quantile_DTRW(alpha = alpha, T = length(X))

    if (!one.sided) {
      decision <- ifelse(obs <= z_theo[2] & obs >= z_theo[1], "DTRW", "NO")
    } else {
      decision <- ifelse(obs <= z_theo[2], "DTRW", "NO")
    }

    return(list(stat = obs, stat_theo = z_theo, decision = decision))

  } else {

    z_obs <- (obs - rec_count_mean_DTRW(length(X), approx = TRUE)) / sqrt(rec_count_var_DTRW(length(X), approx = TRUE))

    if (!one.sided) {
      decision <- ifelse(abs(z_obs) <= qnorm(1 - alpha / 2), "DTRW", "NO")
    } else {
      decision <- ifelse(abs(z_obs) <= qnorm(1 - alpha), "DTRW", "NO")
    }

    p_value <- pnorm(abs(z_obs))

    return(list(stat = z_obs, p_value = p_value, decision = decision))
  }
}

######################

# Test_DTRW_ENT <- function(X, alpha= 0.05) {
#   # X: numeric time X (one realization)
#
#   T <- length(X)
#
#   # Step 1: identify records
#   records <- record_times(X)  # times of records
#   N_T <- length(records)                      # number of records
#
#   # Step 2: Fit log-log regression across subsamples
#   # Split into blocks to see growth of R_t vs t
#   block_sizes <- floor(seq(T/4, T, length.out = 5))  # subsample lengths
#   R_block <- sapply(block_sizes, function(t) {
#     length(record_values(X[1:t]))  ## number of records by blocks
#   })
#
#   df <- data.frame(logT = log(block_sizes), logR = log(R_block))
#
#   fit <- lm(logR ~ logT, data = df)
#
#   beta_hat <- coef(fit)[2]
#   se_beta  <- summary(fit)$coefficients[2,2]
#
#   # Wald test: H0: beta = 0.5
#   z <- (beta_hat - 0.5) / se_beta
#   pval <- 2 * (1 - pnorm(abs(z)))
#
#   return(list(
#     regression = summary(fit),
#     beta_hat = beta_hat,
#     se = se_beta,
#     z = z,
#     pval = pval,
#     decision = ifelse(pval < alpha, "NO", "DTRW")
#   ))
# }

#######################

## FOUR TESTS FOR THE RANDOM WALK HYPOTHESIS - Handa Test 4 - Detection rte of 10% at 5%
# Test_DTRW_LSE<- function(X, alpha= 0.05) {
#   T <- length(X)
#
#   y1=X[-1]
#   y2=X[-length(X)]
#
#   num = sum(y1*y2)
#
#   denom = sum(X^2)
#
#   alpha_hat = num/denom
#
#   T_p = T * (alpha_hat -1 )/sqrt(2)
#
#   #p_val = pnorm(T_p)
#
#   return(list(
#     z = T_p,
#     decision = ifelse(T_p <= -5.79 | T_p > 0.922, "NO", "DTRW")
#   ))
# }

## FOUR TESTS FOR THE RANDOM WALK HYPOTHESIS - Handa Test 8 - Detection rte of 10% at 5%
# Test_DTRW_tratio <- function(X, alpha= 0.05) {
#   T <- length(X)
#
#   y1=X[-1]  ##yt
#   y2=X[-length(X)]  ##yt-1
#
#   alpha_hat = sum(y1*y2)/sum(X^2)
#
#   denom = y1-alpha_hat*y2
#
#   T_p = T * sqrt(sum(X^2)) * (alpha_hat -1 )/sum(denom^2)
#
#   #p_val = pnorm(T_p)
#
#   return(list(
#     z = T_p,
#     "lb" = -5.79,
#     "ub"= 0.922,
#     decision = ifelse(T_p <= -5.79 | T_p > 0.922, "NO", "DTRW") ## Null is DTRW
#   ))
# }
