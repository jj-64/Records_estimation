#' Classical i.i.d. Test Based on Record Counts
#'
#' Tests whether a sequence of observations is i.i.d. based on the number of records,
#' following the test statistic proposed by Arnold.
#'
#' @details
#' The test statistic is:
#' \deqn{Z_T = \frac{N_T - \log(T)}{\sqrt{\log(T)}}}
#' where \eqn{N_T} is the observed number of records in the sample of size \eqn{T}.
#'
#' Under the null hypothesis (i.i.d.), \eqn{Z_T} approximately follows a standard
#' normal distribution.
#'
#' The decision rule:
#' \itemize{
#'   \item Reject H0 (NO) if \eqn{p\_value < \alpha}.
#'   \item Otherwise, fail to reject H0 (classical i.i.d.).
#' }
#'
#' @param X A numeric vector of observations.
#' @param alpha Significance level (default = 0.05).
#'
#' @return A list with:
#' \item{stat}{The observed test statistic.}
#' \item{p_value}{The p-value under the null hypothesis.}
#' \item{decision}{Character string: "Classical" or "NO"}
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' Test_iid_NT(x)
Test_iid_NT <- function(X, alpha = 0.05) {
  T <- length(X)
  obs_stat <- (rec_count(X) - log(T)) / sqrt(log(T))
  p_value <- 1 - pnorm(obs_stat, 0, 1)

  decision <- ifelse(p_value > alpha, "Classical", "NO")

  return(list(
    stat = obs_stat,
    p_value = p_value,
    decision = decision
  ))
}


#' Box–Jenkins i.i.d. Test (Ljung–Box)
#'
#' Tests whether a sequence of observations is i.i.d. using the Ljung–Box test
#' for autocorrelation.
#'
#' @details
#' The Ljung–Box test checks whether any group of autocorrelations of the series
#' is significantly different from zero. It is a common diagnostic test for i.i.d.
#' or white noise.
#'
#' @param X A numeric vector of observations.
#' @param lags Integer. Number of lags to include in the test (default = 10).
#' @param alpha Significance level (default = 0.05).
#' @param type Character. Either "Ljung-Box" (default) or "Box-Pierce".
#'
#' @return A list with:
#' \item{stat}{The Ljung–Box (or Box–Pierce) test statistic.}
#' \item{p_value}{The p-value of the test.}
#' \item{decision}{Character string: "Classical" or "NO"}
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' Test_iid_BoxJenkins(x)
Test_iid_BoxJenkins <- function(X, lags = 10, alpha = 0.05, type = "Ljung-Box") {
  bt <- Box.test(X, lag = lags, type = type)

  decision <- ifelse(bt$p.value > alpha, "Classical", "NO")

  return(list(
    stat = bt$statistic,
    p_value = bt$p.value,
    decision = decision
  ))
}
