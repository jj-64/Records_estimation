
partition_OLD <- function(X, min_expected = 1, warmup = 2, K = NULL) {
  gaps <- rec_gaps(X)

  # drop early gaps if warmup > 0
  if (length(gaps) >= (warmup + 5)) {
    gaps <- gaps[-seq_len(warmup)]
  }

  # Case 1: user specifies number of partitions
  if (!is.null(K)) {
    # use quantiles to split into K groups
    breaks <- unique(ceiling(quantile(gaps, probs = seq(0, 1, length.out = K + 1))))
    #breaks <- c(breaks, Inf)  # ensure full coverage
  } else {
    # Case 2: adaptive rule based on unique values
    num_groups <- length(unique(gaps)) + 1
    breaks <- unique(ceiling(quantile(gaps, probs = seq(0, 1, length.out = num_groups + 1))))
    breaks <- c(breaks, Inf)
  }

  # Initial grouping
  grouped_vec <- cut(gaps, breaks = breaks, include.lowest = TRUE, right = FALSE)
  freq_table <- table(grouped_vec)
  break_points <- breaks   # Track breaks explicitly


  # Adaptively merge small bins if expected counts are too low only if K not fixed
  if (is.null(K)) {
    while (any(freq_table < min_expected) && length(freq_table) > 1) {
      idx <- which.min(freq_table)

      if (idx == length(freq_table)) {
        # merge with left neighbor
        freq_table[idx - 1] <- freq_table[idx - 1] + freq_table[idx]
        freq_table <- freq_table[-idx]
        break_points <- break_points[-idx]
      } else {
        # merge with right neighbor
        freq_table[idx + 1] <- freq_table[idx + 1] + freq_table[idx]
        freq_table <- freq_table[-idx]
        break_points <- break_points[-(idx + 1)]
      }
    }
  }
  # Build interval labels
  interval_labels <- paste0("[", break_points[-length(break_points)], ",", break_points[-1], ")")

  return(list(
    j = break_points[-length(break_points)],  # partition start points
    nk = as.numeric(freq_table),              # frequencies
    labels = interval_labels                  # readable bin labels
  ))
}

#' Partition record gaps adaptively (with padding)
#'
#' Partitions record gaps into bins for Pearson-type tests.
#' Ensures at least 3 partitions by padding with zero-frequency bins if necessary.
#'
#' @param X Numeric vector (time series).
#' @param min_expected Minimum expected count per bin (default = 1).
#' @param warmup Number of initial gaps to drop (default = 2).
#' @param K Optional. If given, force exactly K partitions using quantiles.
#' @param estimated Logical. If to estimate \eqn{\gamma} through minimizing \eqn{\chi^2} (Default = TRUE)
#' @return A list with:
#'   \item{j}{Partition start points}
#'   \item{nk}{Frequencies in each partition}
#'   \item{labels}{Readable bin labels}
#' @export
partition <- function(X, min_expected = 1, warmup = NULL, K = NULL, estimated = TRUE) {
  gaps <- rec_gaps(X)

  min_K = 2+ifelse(estimated, 1, 0)

  # Drop warmup
  # if (length(gaps) >= (warmup + 5)) {
  #   gaps <- gaps[-seq_len(warmup)]
  # }
  if (is.null(warmup)){
    warmup = ceiling(0.3*length(gaps))
    gaps <- gaps[-seq_len(warmup)]
  }

  # Step 1: Build initial breakpoints
  if (!is.null(K)) {
    breaks <- unique(ceiling(quantile(gaps, probs = seq(0, 1, length.out = K + 1))))
  } else {
    num_groups <- min(length(unique(gaps)) + 1 + ifelse(estimated, 1, 0) , min_K)
    breaks <- unique(ceiling(quantile(gaps, probs = seq(0, 1, length.out = num_groups + 1))))
  }
  breaks <- c(breaks, Inf)

  # Step 2: Initial grouping
  grouped_vec <- cut(gaps, breaks = breaks, include.lowest = TRUE, right = FALSE)
  freq_table <- table(grouped_vec)
  break_points <- breaks

  # Step 3: Merge bins if too small (only when K not fixed)
  if (is.null(K)) {
    while (any(freq_table < min_expected) && length(freq_table) > 2) {
      idx <- which.min(freq_table)

      if (idx == length(freq_table)) {
        freq_table[idx - 1] <- freq_table[idx - 1] + freq_table[idx]
        freq_table <- freq_table[-idx]
        break_points <- break_points[-idx]
      } else {
        freq_table[idx + 1] <- freq_table[idx + 1] + freq_table[idx]
        freq_table <- freq_table[-idx]
        break_points <- break_points[-(idx + 1)]
      }
    }
  }

  # Step 4: Ensure at least 3 partitions by padding
  # while (length(freq_table) < 3) {
  #   unique_gaps <- sort(unique(gaps))
  #   candidate <- 1
  #   while (candidate %in% unique_gaps || candidate %in% break_points) {
  #     candidate <- candidate + 1
  #   }
  #
  #   # insert before Inf
  #   freq_table <- c(freq_table, 0)
  #   break_points <- c(break_points[-length(break_points)], candidate, Inf)
  # }

  # Step 5: Build interval labels
  interval_labels <- paste0("[", break_points[-length(break_points)], ",", break_points[-1], ")")

  return(list(
    j = break_points[-length(break_points)], # partition start points
    nk = as.numeric(freq_table),             # frequencies
    labels = interval_labels                 # readable bin labels
  ))
}


# Pearson chi-square test for YNM
Test_YNM_Pearson <- function(X, Partition = NA, gamma = NULL, K=NULL, estimated = TRUE, alpha = 0.05) {

  if (sum(is_rec(X)) <=4) {
              return(list(decision = "NO"))}

  # helper: chi-squared term
  x2_term <- function(m_1, pi, n) (n - pi * m_1)^2 / (pi * m_1)

  # helper: chi-square loss function when estimating gamma
  x2_term_g <- function(gamma, K, nk, j) {
      P_j <- numeric(K)
      if (K>1){
        for (s in 1:(K - 1)) {
          p_j <- (gamma - 1) / gamma^(j[s]:(j[s + 1] - 1))
          P_j[s] <- sum(na.omit(p_j))
        }
        P_j[K] <- 1 - sum(P_j[1:(K - 1)])  # last bin absorbs remainder
      } else {
        P_j =  (gamma - 1) / gamma^(j)
      }
      m_1 <- sum(nk)
      return(sum(x2_term(m_1, P_j, nk)))
                                         }

  # Partition handling
  if (is.na(Partition)[1]) Partition <- partition(X, K=K)

  K <- length(Partition$nk)
  nk <- Partition$nk
  j <- Partition$j
  m_1 <- sum(nk)

  if (K < (2+ifelse(estimated, 1, 0))) {print("Test cannot be performed: partitions K are not sufficient")
  return(list("decision" = NA))}

   # If gamma not provided -> estimate it
  if (is.null(gamma)) {
    gammas <- seq(1.000001, 5, by = 0.01)
    chi_values <- sapply(gammas, function(g) x2_term_g(g, K, nk, j))
    gamma <- gammas[which.min(chi_values)]
    obs_stat <-min(chi_values)
    estimated <- TRUE
                      }
    else{
      obs_stat = x2_term_g(gamma, K, nk, j)
      estimated <- FALSE
    }

  # Adjust degrees of freedom if gamma estimated
  df <- K - 1 - ifelse(estimated, 1, 0)

  # Critical value and p-value
  crit_val <- qchisq(p = 1 - alpha, df = df)
  p_value <- 1 - pchisq(q = obs_stat, df = df)

  # Decision
  decision <- ifelse(obs_stat < crit_val, "YNM", "NO")

  return(list(
    stat = obs_stat,
    p_value = p_value,
    gamma = gamma,
    df = df,
    decision = decision
  ))
}


## Test lisse
Test_YNM_Smooth <- function(X, alpha = 0.05) {

  # Meixner-like polynomials
  h2 <- function(x, a) x*(x-1) - 4*a*x + 2*a^2
  h3 <- function(x, a) x*(x-1)*(x-2) - 9*a*x*(x-1) + 18*x*a^2 - 6*a^3
  h4 <- function(x, a) 24*choose(n=x, k=4) - 96*a*choose(n=x, k=3) +
    144*a^2*choose(n=x, k=2) - 96*x*a^3 + 24*a^4
  h5 <- function(x, a) 120*choose(n=x, k=5) - 600*a*choose(n=x, k=4) +
    1200*a^2*choose(n=x, k=3) - 1200*a^3*choose(n=x, k=2) +
    600*a^4*x - 120*a^5

  # record count and mean gap
  m <- rec_count(X)
  gaps <- rec_gaps(X)
  dbar <- sum(gaps) / (m - 1)

  # helper for variance scaling
  Vr <- function(r, m, dbar, Hsum) {
    scale <- ( (m-1) * factorial(r)^2 * ((dbar-1)^2 + (dbar-1))^r )^(-0.5)
    scale * Hsum
  }

  # correction factor
  correct_factor <- function(m, dbar) {
    1 + 3.643/(m-1) - 2.314/sqrt(m-1) -
      0.447 / sqrt((m-1)*(dbar-1)/dbar)
  }

  # center gaps
  x <- gaps - 1
  a <- dbar - 1

  # compute sums of polynomials
  H2 <- sum(h2(x, a))
  H3 <- sum(h3(x, a))
  H4 <- sum(h4(x, a))
  H5 <- sum(h5(x, a))

  # compute normalized V statistics
  V2 <- Vr(2, m, dbar, H2)
  V3 <- Vr(3, m, dbar, H3)
  V4 <- Vr(4, m, dbar, H4)
  V5 <- Vr(5, m, dbar, H5)

  # observed test stat
  Sk <- V2^2 + V3^2 + V4^2 + V5^2
  obs_stat <- Sk * correct_factor(m, dbar)

  # chi-squared approximation with k=4 df
  k <- 4
  crit_val <- qchisq(1 - alpha, df = k)
  p_value <- 1 - pchisq(obs_stat, df = k)

  decision <- ifelse(obs_stat <= crit_val, "YNM", "NO")

  return(list(
    stat = obs_stat,
    p_value = p_value,
    df = k,
    decision = decision
  ))
}



#######---------------  Test_YNM_NT -------#########

Quantile_YNM=function(T,gamma, alpha= 0.05){
  Prob = 0

  for(i in 1:T){
    Prob[i]= rec_count_dist_YNM(m=i, T=T,gamma=gamma)
  }
  CDF = cumsum(Prob)  ## cumulative distribution

  #plot(x=1:T, CDF)

  return(c(which.min(abs(CDF-(alpha/2))) , which.min(abs(CDF-(1-alpha/2))) ))
}


#' YNM's Record Test (Number of Records Statistic)
#'
#' Implements YNM's test for sequences based on the distribution
#' of the number of records.
#'
#' @details
#' The test uses the number of observed records \eqn{N_T} and compares it against
#' the theoretical distribution under the null hypothesis.
#'
#' If \eqn{\gamma} is unknown, it is estimated using \code{\link{estimate_YNM_mle_indicator}}
#' with its variance.
#'
#' The quantiles (theoretical) of the record distribution are obtained from the cumulative distribution:
#' \deqn{F(m) = \sum_{i=1}^{m} P(N_T = i)} as in \code{\link{rec_count_dist_YNM}}
#'
#' @param X A numeric vector (time series).
#' @param gamma Optional. The power parameter of the YNM-Nevzorov Model. If not provided, it will be estimated.
#' @param alpha Significance level (default = 0.05).
#'
#' @return A list with:
#' \item{stat}{Observed number of records.}
#' \item{stat_theo}{Theoretical quantile interval at level \eqn{1-\alpha}.}
#' \item{gamma_hat}{Estimated or provided \eqn{\gamma}.}
#' \item{decision}{Character string: "NO" (reject null) or "YNM" (fail to reject).}
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- LDM_series(T=50, theta=0.5, dist="gumbel", loc=0, scale=1)
#' Test_YNM_NT(x)
#' y= VGAM::rgumbel(50,0,1)
#' Test_YNM_NT(y)
Test_YNM_NT <- function(X, gamma = NA, alpha = 0.05) {
  T <- length(X)
  obs <- rec_count(X)
  v_gamma = NA_real_

  # --- Estimate gamma if not provided ---
  if (is.na(gamma)) {
    estimated =  estimate_YNM_mle_indicator(X, variance = TRUE, approximate=FALSE, min = 1.01, max=5, step=0.01)
    gamma <- estimated$param
    v_gamma = estimated$variance
  }

  # --- Variance  ---
  if(is.na(v_gamma)){
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
        v_gamma <- 1 / I
      }
  }

  # --- test statistic ---
  z_theo <- (gamma-1) / sqrt(v_gamma)

  # --- First check using normal approximation: Ho: Gamma =1 ---
  if ( z_theo <= qnorm(1-alpha/2) & z_theo >= qnorm(alpha/2) ) {  ## Gamma = 1
    decision <- "NO"
    stat_theo <- z_theo
  } else {
    # --- Fallback to exact quantiles ---
    probs <- vapply(1:T, function(i) rec_count_dist_YNM(m = i, T = T, gamma = gamma), numeric(1))
    cdf <- cumsum(probs)

    lower <- which.min(abs(cdf - (alpha / 2)))  ##Quantile_YNM(T=T, gamma=gamma, alpha= alpha)[1]
    upper <- which.min(abs(cdf - (1 - alpha / 2)))  ###Quantile_YNM(T=T, gamma=gamma, alpha= alpha)[2]

    stat_theo <- c(lower, upper)
    decision <- ifelse(obs >= max(lower, 1) && obs <= upper, "YNM", "NO")
  }

  return(list(
    stat = obs,
    stat_theo = stat_theo,
    gamma_hat = gamma,
    decision = decision
  ))
}


#######---------------  Test_YNM_Geom -------------- #########

#The usual test assumes Gaps are i.i.d. Geom(p). If gaps are nonstationary (e.g. short early gaps, longer later gaps) the single-sample mean
#is dominated by early small gaps → p biased high → expected tail mass under the null is underestimated → observed long gaps look unsurprising → test fails to reject.
#So the issue is not the chi-square per se but the (false) assumption of stationarity/homogeneity of gap distribution across time.


#' Test for YNM–Nevzorov Geometric Record Gaps
#'
#' Performs a goodness-of-fit test for the YNM–Nevzorov model by examining
#' the distribution of record time gaps in a sequence of observations.
#' The procedure estimates the geometric parameter, computes a chi-squared
#' statistic on grouped frequencies, and evaluates whether the observed record
#' gaps are consistent with a geometric distribution.
#'
#' @details
#' Let \eqn{G_i} denote the observed record gaps. Under the YNM-Nevzorov
#' model, the gaps are approximately geometrically distributed with parameter
#' \eqn{\hat{p} = 1 / \bar{G}}. The test proceeds as follows:
#'
#' 0. Remove warm-up gaps (if NULL, 30% of gaps are removed)
#' 1. Estimate \eqn{p} and the implied \eqn{\hat{\gamma} = 1 / (1 - \hat{p})}.
#' 2. Compute expected frequencies for the first \eqn{K-1} categories and
#' group all larger gaps into the \eqn{K}-th bin.
#' 3. Form the Pearson chi-squared statistic:
#'    \deqn{ \chi^2 = \sum_{k=1}^K \frac{(O_k - E_k)^2}{E_k}, }
#'    with \eqn{df = (K - 1) - 1} degrees of freedom (adjusted for the
#'    estimated parameter).
#' 4. Return decision "YNM" if the null hypothesis of geometric gaps
#' is not rejected, and "NO" otherwise.
#'
#' A confidence interval for {\eqn{\gamma}} is also computed using a
#' normal approximation:
#' \deqn{ \text{Var}(\hat{\gamma}) = \frac{\hat{\gamma} (\hat{\gamma} - 1)^2}{n}. }
#'
#' @param X Numeric vector of observations.
#' @param alpha Numeric, significance level (default = 0.05).
#' @param K Integer, number of categories for chi-squared grouping.
#'   If NULL, defaults to \eqn{min(4, length(gaps))}.
#' @param warmup Integer, number of initial gaps to discard (default = NULL).
#' @param obs_type String. "all" if data provided is the whole series \eqn{X_t} or
#' "records" if the underlying series is \eqn{R_n}. In this case, the parameter
#' record_times must be provided.
#' @param record_times Numeric vector of the occurence times of records. (Default is NA).
#' Forced in case "obs_type" = "records"
#' @return A list with elements:
#' \item{obs_counts}{Observed counts per category.}
#' \item{exp_counts}{Expected counts per category under fitted geometric law.}
#' \item{stat}{Chi-squared test statistic.}
#' \item{df}{Degrees of freedom.}
#' \item{p_value}{P-value of the chi-squared test.}
#' \item{p_hat}{Estimated geometric parameter.}
#' \item{gamma_hat}{Estimated YNM parameter \eqn{\gamma}.}
#' \item{v_gamma_hat}{Estimated variance of \eqn{\hat{\gamma}}.}
#' \item{CL}{Confidence interval for \eqn{\gamma}.}
#' \item{decision}{Decision: "YNM" if not rejected, "NO" otherwise.}
#' @export
#' @examples
#' set.seed(123)
#' X <- YNM_series(T = 50, dist = "gumbel", gamma = 1.2, loc = 0, scale= 1)
#' Test_YNM_Geom(X, alpha = 0.05, K = 4, warmup=2)
#'
#' # $obs_counts
#' # [1] 3 3 1 4
#'
#' # $exp_counts
#'
#' # [1] 2.69 2.03 1.53 4.74
#'
#' # $stat
#' # [1] 0.8008871
#'
#' # $df
#' # [1] 2
#'
#' # $p_value
#' # [1] 0.6700228
#'
#' # $p_hat
#' # [1] 0.2444444
#'
#' # $gamma_hat
#' # [1] 1.323529
#'
#' # $v_gamma_hat
#' # [1] 0.01259414
#'
#' # $CL
#' # [1] 1.138938 1.508121
#'
#' # $decision
#' # [1] "YNM"
Test_YNM_Geom <- function(X, alpha=0.05, K=NULL, warmup=NULL, obs_type = c("all", "records"), record_times=NA) {

  obs_type <- match.arg(obs_type)

  ## Force record times if only records are used
  if (obs_type == "records" && length(record_times) != length(X)) {
    stop("Number of records is not the same as record times (or not provided)")
  } else if (obs_type == "records" && is.numeric(record_times[1])) {
    gaps <- diff(record_times)
  }

  if (sum(is_rec(X)) <=4) {
    return(list(decision = "NO"))}

  ## if vector Xt is provided
  if(obs_type == "all") gaps = rec_gaps(X)

  ## Warmup
  if (is.null(warmup) ){ warmup = ceiling(1/3* length(gaps)) }
  #if (warmup <2 ){ warmup =2}
  if (warmup > 0 ){ gaps <- gaps[-seq_len(warmup)]}

  # Total number of observations
  n <- length(gaps)
  if( n<=1 ) {return(list("decision"="NO")) }

  ## geometric parameter
  p_hat <- 1/mean(gaps)
  gamma=1/(1-p_hat) ## gamma
  v_gamma_hat = gamma*(gamma-1)^2/n
  CL = bounds(gamma,qnorm(1-alpha),v_gamma_hat)  ## one sided test

  # theor_var = mean(gaps)*(mean(gaps)-1)
  # obs_var = var(gaps)
  # (n-1) * obs_var/theor_var <= qchisq(1-alpha, df = n-1)

  # define bins 1:(k-1) and K = ">=k"
  if (is.null(K)) {K= max(3, length(gaps) )}  ## min(4, length(ga))}
  if(K<2) {K=3}  ## to avoid df = K-2 <0
  obs_counts <- tabulate(pmin(gaps, K), nbins = K)

  # Expected frequencies/probabilities using the geometric distribution with the estimated alpha
  exp_probs <- numeric(K)
  for (k in 1:(K-1)) exp_probs[k] <- (1 - p_hat)^(k-1) * p_hat
  exp_probs[K] <- 1 - sum(exp_probs[1:(K-1)])  # tail mass

  exp_counts <- n * exp_probs

  # If any expected < 5, warn and consider decreasing k or use MC
  if (any(exp_counts < 5)) warning("Some expected counts < 5; consider pooling or MC test.")

  chi2 <- sum((obs_counts - exp_counts)^2 / exp_counts)
  df <- (K - 1) - 1   # bins-1 - params estimated

  p_value <- pchisq(chi2, df = df, lower.tail = FALSE)
  decision = ifelse(p_value>=alpha & 1 <=CL[1], "YNM", "NO") #& CL[1]>1
  # if(p_value < alpha) {
  #   decision = "NO"
  # } else if( CL[1] >1) {
  #   decision= "YNM"
  # } else {decision = "NO"}

  return(list(obs_counts = obs_counts,
              exp_counts = round(exp_counts,2),
              "stat" = chi2,
              df = df,
              "p_value" = p_value,
              "p_hat" = p_hat,
              "gamma_hat"=gamma ,
              "v_gamma_hat" = v_gamma_hat ,
              "CL" = CL,
              "decision"=decision))
}


############### IGNORE###########################################

# dispersion test for geometric gaps: asymptotic chi-square + MC p-value
# dispersion_geom_test <- function(X, alpha=0.05, B = 5000, seed = NULL, return_sim = FALSE) {
#
#    gaps = rec_gaps (X)
#   if (!is.null(seed)) set.seed(seed)
#
#   if (any(gaps < 1)) return(list(decision= "NO"))
#
#   n <- length(gaps)
#
#   if (n < 5) warning("Small n: chi-square approx will be poor; prefer Monte Carlo.")
#
#   Gbar <- mean(gaps)
#   s2 <- var(gaps) # sample variance uses denominator (n-1)
#
#   # Prevent degenerate case Gbar <= 1 (won't happen unless all gaps=1)
#   if (Gbar <= 1) {
#     # variance under estimated p is zero or negative; handle separately
#     return(list(error = "Mean of gaps <= 1 (all gaps=1?), cannot compute dispersion test.",
#                 Gbar = Gbar, s2 = s2))
#   }
#
#   # test statistic
#   D <- (n - 1) * s2 / (Gbar * (Gbar - 1))
#
#   # asymptotic p-value using chi-square_{n-1}
#   p_chisq <- pchisq(D, df = n - 1, lower.tail = FALSE)
#
#   res <- list(n = n, Gbar = Gbar, s2 = s2, D = D,
#               p_chisq= p_chisq,
#               decision = ifelse(p_chisq>alpha, "Geom", "NO"))
#   return(res)
# }
