#####------------------  LDM_Regression Hubert ------------------  ####
#' Hypothesis Test for Linear Drift Model using Robust Regression
#'
#' Tests whether the observed sequence \code{X} follows the Linear Drift Model (LDM)
#' by fitting a robust linear regression of \code{X} on time and testing
#' if the slope (drift parameter) is significantly different from zero.
#'
#' @param X Numeric vector of observations.
#' @param alpha Numeric, significance level (default = 0.05).
#' @param RSq Numeric, minimum adjusted R-squared required to accept the LDM
#'   hypothesis (default = 0.8).
#' @param obs_type String. "all" if data provided is the whole series \eqn{X_t} or
#' "records" if the underlying series is \eqn{R_n}. In this case, the parameter
#' record_times must be provided.
#' @param record_times Numeric vector of the occurence times of records. (Default is NA).
#' Forced in case "obs_type" = "records"
#' @return A list with components:
#' \describe{
#'   \item{stat}{Estimated slope coefficient (drift parameter).}
#'   \item{p_value}{p_value for testing \eqn{H_0: \theta = 0}.}
#'   \item{RS}{Adjusted R-squared of the robust regression.}
#'   and adjusted R-squared exceeds the threshold, \code{"NO"} otherwise.}
#'   \item{MSE}{Mean squared error of the fitted model (normalized by \eqn{\sum X^2}).}
#'   \item{STD}{Estimated standard error of the slope coefficient.}
#'   \item{decision}{Decision: \code{"LDM"} if the null hypothesis is rejected
#' }
#'
#' @details
#' The test fits a robust regression model of the form:
#' \deqn{ X_t = \theta (t - \bar{t}) + \varepsilon_t }
#' where time \eqn{t} and observations \eqn{X_t} are centralized.
#' The decision rule is:
#' \enumerate{
#'   \item Reject \eqn{H_0: \theta = 0} if the p_value is below \code{alpha}.
#'   \item Accept LDM only if the adjusted R-squared exceeds \code{RSq}.
#' }
#'
#' @examples
#' set.seed(123)
#' t <- 1:50
#' X <- 0.3 * t + rnorm(50, sd = 5)  # Linear drift with noise
#' Test_LDM_Regression(X, alpha = 0.05, RSq = 0.7)
#' @export
Test_LDM_Regression <- function(X, alpha = 0.05, RSq = 0.8, obs_type = c("all"), record_times=NA) {

  ##Force to provide record times in case only records are used
  if (obs_type == "records" & length(record_times) != length(X) ){
    stop("Number of records is not the same as record times (Or not provided)")
  } else  if (obs_type == "records" & is.numeric(record_times[1])){
    t= record_times
  }

  if (obs_type == "all") {t <- seq_along(X)}

  # Centralize data
  t <- t - mean(t)
  X <- X - mean(X)

  # Fit robust regression without intercept
  train <- data.frame(X = X, t = t, rec= ifelse(is_rec(X) ==1, rec_times(X), 0))
  robust_model <- robustbase::lmrob(X ~ t - 1, data = train, method = "MM")

  # Extract coefficient and standard error
  coef_summary <- coef(summary(robust_model))
  coefficient <- robust_model$coefficients
  std_error <- coef_summary[1, 2]

  # Compute p_value for H0: theta = 0
  t_stat <- coef_summary[1, 1] / std_error
  p_value <- 1 - pt(t_stat, df = robust_model$df.residual)

  # Compute fit quality
  MSE <- sum(robust_model$residuals)^2 / sum(X^2)
  RS <- summary(robust_model)$adj.r.squared

  # Decision rule
  #decision <- ifelse(p_value < alpha & RS >= RSq, "LDM", "NO")
  decision <- ifelse(p_value < alpha & RS >= RSq, "LDM", "NO")

  return(list(
    stat = coefficient,
    p_value = p_value,
    RS = RS,
    MSE = MSE,
    STD = std_error,
    decision = decision
  ))
}

#########------------------  Test based on observed number of records ------------------ ################

#' Quantile Function for LDM Distribution of Record Numbers in case of Gumbel underlying distribution
#'
#' Computes the lower and upper quantile bounds for the distribution
#' of the number of records under the Linear Drift Model (LDM) under Gumble (µ, σ).
#'
#' @param T Integer, sample size (number of observations).
#' @param theta Numeric, drift parameter \eqn{\theta > 0}.
#' @param scale Numeric, scale parameter σ (default = 1).
#' @param alpha Numeric, significance level for the two-sided interval (default = 0.05).
#'
#' @return A numeric vector of length 2 giving the lower and upper quantile indices
#' corresponding to probabilities \eqn{\alpha/2} and \eqn{1-\alpha/2}.
#'
#' @details
#' The function computes the distribution of the number of records under the LDM,
#' using the Stirling numbers of the second kind (\code{Stirling_2nd_LDM}) and
#' the probability mass function \code{rec_count_dist_LDM}.
#' The cumulative distribution function (CDF) is then compared to the
#' desired quantile levels.
#'
#' @seealso \code{\link{Test_LDM_NT}} for hypothesis testing of the number of records.
#' @export
#' @examples
#' Quantile_LDM(T = 20, theta = 0.5, alpha = 0.05)
Quantile_LDM <- function(T, theta, scale = 1, alpha = 0.05) {
  S <- Stirling_2nd_LDM(T, theta, scale)

  # Vectorized probability computation
  Prob <- vapply(1:T, function(i) {
    rec_count_dist_LDM(m = i, T = T, theta = theta, scale = scale, s = S)
  }, numeric(1))

  # Cumulative distribution
  CDF <- cumsum(Prob)

  # Return indices closest to alpha/2 and 1 - alpha/2
  return(c(
    which.min(abs(CDF - alpha / 2)),
    which.min(abs(CDF - (1 - alpha / 2)))
  ))
}



#' Hypothesis Test for Number of Records under LDM and Gumbel underlying distribution
#'
#' Tests whether the observed number of records in a sequence \code{X}
#' is consistent with the Linear Drift Model (LDM).
#'
#' @param X Numeric vector of observations.
#' @param alpha Numeric, significance level (default = 0.05).
#'
#' @return A list with components:
#' \describe{
#'   \item{stat}{Observed number of records.}
#'   \item{stat_theo}{Theoretical expected statistic under LDM
#'   (either z-statistic or quantile interval).}
#'   \item{theta_hat}{Estimated drift parameter \eqn{\theta}.}
#'   \item{v_theta_hat}{Estimated variance for drift parameter \eqn{\theta}.}
#'   \item{decision}{Decision: \code{"LDM"} if accepted, \code{"NO"} otherwise.}
#' }
#'
#' @details
#' The test proceeds in two steps:
#' \enumerate{
#'   \item Estimate the drift parameter \eqn{\theta} using
#'   \code{Estim_theta_indicator}.
#'   \item Compute a z-test statistic \eqn{\hat{\theta} / \sqrt{\mathrm{Var}(\hat{\theta})}}.
#'   If the statistic lies outside the acceptance region (i.e it is significantly different from zero), compute the quantile
#'   interval of the number of records using \code{Quantile_LDM}.
#' }
#'
#' @seealso \code{\link{Quantile_LDM}} for quantile computation.
#' @export
#' @examples
#' set.seed(123)
#' X <- rnorm(50, mean = 0.2 * (1:50))  # Linear drift
#' Test_LDM_NT(X, alpha = 0.05)
Test_LDM_NT <- function(X, alpha = 0.05) {
  obs <- rec_count(X)
  scale <- 1

  # Estimate theta and variance
  estimated = estimate_LDM_mle_indicator(X=X, variance = TRUE, scale=1, min = 0.0001, max = 5, step = 0.01)
  theta <- estimated$param
  v_theta <- estimated$variance
  z_theo <- theta / sqrt(v_theta)

  # Decision rule
  if (z_theo <= qnorm(1 - alpha / 2) && z_theo >= qnorm(alpha / 2)) {
    decision <- "NO"
  } else {
    z_theo <- Quantile_LDM(T = length(X), theta = theta, scale = scale, alpha = alpha)
    decision <- ifelse(obs <= z_theo[2] && obs >= max(z_theo[1], 1), "LDM", "NO")
  }

  return(list(
    stat = obs,
    stat_theo = z_theo,
    theta = theta,
    variance = v_theta,
    decision = decision
  ))
}

########----------------- Two-stage based on Şen Slope  ---------------- ##########

#' Sequential Test for the Linear Drift Model (LDM)
#'
#' Implements a two-stage sequential testing procedure (based on Şen-type tests)
#' to detect the presence of a linear drift in a time series.
#'
#' The procedure partitions the data into blocks and checks whether block means
#' and slopes are consistent with a monotonic non-zero linear drift. If stage 1 accepts, stage 2
#' refines the decision by testing slope significance across partitions.
#'
#' @param X Numeric vector of observations.
#' @param time Numeric vector of the record times. (default = NA means we are using all observed X series)
#' @param alpha Numeric, significance level (default = 0.05).
#' @param pooled Logical, if \code{TRUE}, pooled variance across blocks is used
#'   in stage 1 for the noise variance estimation. Otherwise, variance is estimated directly
#'   from residuals (default = \code{FALSE}).
#' @return A list with results from both stages:
#' \describe{
#'   \item{block_size}{Size of each partitioned block.}
#'   \item{means}{Block means.}
#'   \item{slopes}{Estimated block-to-block slopes.}
#'   \item{var_x}{Estimated variance of the noise process.}
#'   \item{Z}{Test statistic from stage 1.}
#'   \item{Var_Z}{Estimated variance of \code{Z}.}
#'   \item{stat}{Standardized stage 1 test statistic.}
#'   \item{p_value}{p_value of the stage 1 test.}
#'   \item{decision}{Decision of stage 1: \code{"Sameslope"} or \code{"NO"}.}
#'   \item{z_A}{Stage 2 standardized slope for first block comparison.}
#'   \item{z_B}{Stage 2 standardized slope for third block comparison.}
#'   \item{slope_SD}{Estimated standard deviation of slope differences (stage 2).}
#'   \item{decision}{Final decision: \code{"LDM"} or \code{"NO"}.}
#' }
#'
#' @details
#' Stage 1 partitions the series into 3 blocks of size \eqn{m = \lfloor n/3 \rfloor}.
#' It computes block means and slopes and tests whether the central block mean
#' aligns with the average of the first and third block means:
#' \deqn{ Z = \frac{2 \bar{X}_2 - \bar{X}_1 - \bar{X}_3}{m} }
#' with variance estimated either via pooling (\code{pooled=TRUE}) or directly.
#'
#' Stage 2 tests whether the slopes between blocks are significantly different
#' from zero:
#' \deqn{ z_A = \frac{s_1}{\hat{\sigma}_s}, \quad z_B = \frac{s_2}{\hat{\sigma}_s} }
#' where \eqn{\hat{\sigma}_s = \sqrt{2 \hat{\sigma}_y / m^3}} and \eqn{Y} is the underlying distribution of the LDM process.
#'
#' The procedure stops at stage 1 if no drift is detected.
#'
#' @examples
#' set.seed(123)
#' t <- 1:60
#' X <- 0.1 * t + rnorm(60, sd = 1)
#' Test_LDM_Sequential(X, alpha = 0.05)
#'
#'
#' @export
Test_LDM_Sequential <- function(X, time = NA,alpha = 0.05) {
  if(length(X)<4) {return(decision = "NO")}

  # Stage 1
  res1 <- Test_LDM_Sen_stage1(X, time=time, alpha = 2 * alpha, pooled = FALSE)  # one-sided

  if (res1$decision1 == "NO") {
    names(res1)[which(names(res1) == 'decision1')] = "decision"
    return(res1)
  }

  # Stage 2
  res2 <- Test_LDM_Sen_stage2(res1, alpha = alpha)

  return(c(res1, res2))
}

#' @rdname Test_LDM_Sequential
#' @param pooled boolean (Default = FALSE), if variance is assumed pooled
#' @export
Test_LDM_Sen_stage1 <- function(X, time=NA, alpha = 0.05, pooled = FALSE) {

  if(is.na(time[1])) {time = 1:length(X)}
  if (length(X) != length(time)) stop("X and time must have the same length")

  #n <- length(X)
  n_blocks <- 3

  # Sort by time
  ord <- order(time)
  X <- X[ord]
  time <- time[ord]

  block_means <- numeric(n_blocks)
  blocks <- vector("list", n_blocks)

  # Partition data into equal blocks
  m <- floor(length(X) / n_blocks)
  for (i in 1:n_blocks) {
    start <- (i - 1) * m + 1
    end <- ifelse( i< n_blocks, i * m, length(X) )
    blocks[[i]] <- X[start:end]
    block_means[i] <- mean(blocks[[i]])}


  # Slopes between consecutive block means
  slopes <- numeric(n_blocks)
  for (i in 1:(n_blocks - 1)) {
    slopes[i + 1] <- (block_means[i + 1] - block_means[i]) / m
  }

  # Residual estimates (noise)
  x_hat <- list(
    blocks[[1]] - slopes[2] * seq_len(m),
    blocks[[2]] - slopes[2] * ((m + 1):(2 * m)),
    blocks[[3]] - slopes[3] * (((2 * m) + 1):(length(X)))
  )

  # Stage 1 statistic (equality of slopes)

  z <- (2 * block_means[2] - block_means[1] - block_means[3])

  # Variance of Z
  if (pooled) {
    s_sq <- sapply(x_hat, var)
    var_x <- (m - 1) * sum(s_sq) / (3 * m - 3)
  } else {
    var_x <- var(unlist(x_hat))
  }

  # Standardized test statistic
  var_z <- (6 * var_x) / m^2
  T_p <- z / sqrt(var_z)

  # p_value
  p_value <- if (max(time) >= 30) {
    1 - pnorm(abs(T_p))
  } else {
    1 - pt(q = abs(T_p), df = n - 3)
  }

  # Decision: same slope or reject
  decision <- ifelse(p_value > alpha / 2, "Sameslope", "NO")

  return(list(
    block_size = m,
    #block_times = block_times,
    means = block_means,
    slopes = slopes[-1],
    var_x = var_x,
    Z = z,
    Var_Z = var_z,
    stat = T_p,
    p_value = p_value,
    decision1 = decision
  ))
}

#' @rdname Test_LDM_Sequential
#' @param stage1 output of \code{\link{Test_LDM_Sen_stage1}}
#' @export
Test_LDM_Sen_stage2 <- function(stage1, alpha = 0.05) {
  m <- stage1$block_size
  s1 <- stage1$slopes[1]
  s2 <- stage1$slopes[2]
  var_x <- stage1$var_x

  # Standard deviation of slope estimator
  sigma_s <- sqrt(2 * var_x / m^3)

  # Test statistics
  z_A <- s1 / sigma_s
  z_B <- s2 / sigma_s

  # Decision
  decision <- ifelse(abs(z_A) <= qnorm(1 - alpha / 2) | abs(z_B) <= qnorm(1 - alpha / 2),
                "NO", "LDM")

  return(list(z_A = z_A, z_B = z_B, slope_SD = sigma_s, decision = decision))
}


##### ------------------  Based on Sen-Slope ------------------ #############

Test_LDM_Sen = function(X, alpha=0.05){

  n= length(X)

  m <- floor(n/2)
  y1 <- X[1:m]
  y2 <- X[(m+1):(2*m)]
  if (length(y1) != length(y2)) y2 = y2[1:length(y1)]

  # block means
  y1_bar <- mean(y1)
  y2_bar <- mean(y2)

  # slope and intercept
  s <- (y2_bar - y1_bar) / m
  a <- mean(X) - s * mean(1:n)

  ## Fit
  fit <- a + (1:n) * s
  noise <- X - fit

  # correlation:order-stat correlation
  rho <- cor(sort(y1), sort(y2))     # paired correlation between blocks: closer to one, then no trend bcz on 1:1 line

  ## Noise estimation
  sigma_s = sqrt(8 / n^3 * var(X)) #* sqrt((1 - abs(rho)))
  #sigma_s_theo = 8* var(x_hat)* (2)/(n^3)
  #CL_low_H0 = qnorm(alpha/2)* sigma_s  #-1.959 for alpha =0.05
  #CL_up_H0 = qnorm(1-alpha/2)* sigma_s
  #"CL_H0" = c(CL_low_H0, CL_up_H0)
  #decision= ifelse(s >= CL_up_H0 | s<= CL_low_H0, "LDM","NO")
  z=s/sigma_s
  decision = ifelse(abs(z) >= qnorm(1-alpha/2), "LDM", "NO")

   return(list( "stat" = z,"theta_hat"=s,"theta_SD"=sigma_s,"decision"=decision))
  #abline(a=a, b=s)
  #plot(sort(y2),sort(y1))
  #abline(a=0, b=1)
            }

# Test_LDM_Sen_stage1 <- function(X, alpha=0.05, pooled = FALSE) {
#   n <- length(X)
#   n_blocks <- 3    #floor(length(X) / m)
#   m = floor(n/n_blocks)
#
#   rho_all = numeric(n_blocks)
#   block_means <- numeric(n_blocks)
#
#   slopes <- 0
#   var_slopes = 0
#   blocks= list()
#   x_hat = list()
#
#   var_slope = function(n,m, sigma2, correlation) {2/3 * sigma2/m^3 * (1-correlation) }
#
#   for (i in 1:n_blocks) {
#     start <- (i - 1) * m + 1
#     end <- i * m
#     block <- X[start:end]
#
#     # Store block mean
#     block_means[i] <- mean(block)
#     blocks[[i]] = block
#   }
#
#   ## slopes
#   for (i in 1:(n_blocks-1)){
#     slopes[i+1] = (block_means[i+1] - block_means[i])/m
#   }
#
#   rho_all[1] =  cor(sort(blocks[[1]]), sort(blocks[[2]])  )  # 1 and 2
#   rho_all[2] =  cor(sort(blocks[[2]]), sort(blocks[[3]])  )  # 2 and 3
#   rho_all[3] =  cor(sort(blocks[[1]]), sort(blocks[[3]])  )  # 1 and 3
#
#   ## Estimation of x_hat, the noise
#   x_hat[[1]] = blocks[[1]] - slopes[2] * c(1:m)
#   x_hat[[2]] = blocks[[2]] - slopes[2] * c((m+1):(2*m))
#   x_hat[[3]] = blocks[[3]] - slopes[3] * c(((2*m)+1):(3*m))
#
#
#   # Results
#   w <- data.frame(
#     block = 1:n_blocks,
#     mean_value = block_means,
#     slopes = slopes,
#     var_slopes = var_slopes
#   )
#
#   # 3. Compute z
#   z = (2*w[2,2]-w[1,2]-w[3,2])/m
#
#   # 4. Estimate variance of z, xhat and V(x)
#
#   if(pooled){
#     # Option A: pooled variance across all three groups (assuming equal variances)
#       s1_sq <- var(x_hat[[1]]) #var(blocks[[1]])
#       s2_sq <- var(x_hat[[2]]) #var(blocks[[2]])
#       s3_sq <- var(x_hat[[3]])# var(blocks[[3]])
#
#       var_x <- (m-1) * (s1_sq + s2_sq + s3_sq) / (3*m - 3)  ## variance of X
#       var_z <- (6 *var_x) / m^3
#   # } else if (unpooled) {
#   #   # Option B: group-wise (if you don’t assume equal variances)
#   #     var_z <- (4*s2_sq + s1_sq + s3_sq) / m^3
#   } else{
#     var_x = var(unlist(x_hat))
#     var_z =  6 * var_x/m^3 #2 * (var(x_hat)/m^3) * (3- 2* rho_all[1] - 2*rho_all[2] + rho_all[3])  #
#   }
#
#   ## estimated variance of slope
#   #var_slopes[2] = var_slope(n=n, m=m, sigma2 = var_x, correlation = rho_all[1])
#   #var_slopes[3] = var_slope(n=n, m=m, sigma2 = var_x, correlation = rho_all[2])
#
#   #5.Compute test statistic
#   T_p <- z / sqrt(var_z)
#
#   #p_value <- 2 * pt(-abs(T_p), df = n - 3)
#   if(n >=30){
#   p_value = 1-pnorm(abs(T_p))
#   }else{
#   p_value = 1-pt(q=abs(T_p), df = n-3)
#   }
#
#   #decision= ifelse((T_p <= qnorm(alpha,0,1) | T_p >= qnorm(1-alpha,0,1)), "NO", "SameSlope")
#   decision = ifelse (p_value > alpha/2 , "Sameslope", "NO")
#   #decision= ifelse( z >= qnorm(alpha/2,0,1)*sqrt(var_z) &  z <=  qnorm(1-alpha/2,0,1)*sqrt(var_z), "LDM", "NO")
#   return(list("block_size"= m,"means"= w$mean_value, "slopes"= w$slopes[-1],"var_x" = var_x,"Z"=z, "Var_Z" = var_z, "stat"=  T_p,"p_value" =  p_value, "decision" = decision))
# }
#
# Test_LDM_Sen_stage2 = function(stage1, alpha=0.05){
#
#   ## as if testing Y1-bar = y3_bar or if slope from 1st and third partition is zero
#   m = stage1$block_size
#
#   ## slopes
#   s1 = stage1$slopes[1]
#   s2 = stage1$slopes[2]
#   #s_hat = mean(s1, s2)
#
#   ## var (s-hat)
#   var_x = stage1$var_x
#   sigma_s = sqrt(2*var_x / m^3)
#
#   ## Test
#   z_A=  s1/sigma_s
#   z_B=  s2/sigma_s
#   decision = ifelse(abs(z_A) <= qnorm(1-alpha/2) | abs(z_B) <= qnorm(1-alpha/2),  "NO", "LDM")
#
#   return(list( "z_A" = z_A, "z_B" = z_B, "slope_SD"=sigma_s,"decision"=decision))
# }
#
# Test_LDM_Sequential <- function(X, alpha = 0.05) {
#   t= seq_along(X)
#   m = floor(length(X)/3)
#
#   # Step 1: Run Test_LDM_Sen_stage1
#   res1 <- Test_LDM_Sen_stage1(X, alpha = 2*alpha, pooled = FALSE)  ## one sided test
#
#   if (res1$decision == "NO") {
#     # Stop immediately
#     return(res1)
#   }
#
#   # Step 2: Only run if stage1 says "LDM"
#   n= floor(length(X)/2)
#   # res2_A <- Test_LDM_Sen(X[1:n], alpha = alpha)
#   # res2_B <- Test_LDM_Sen(X[(n +1):length(X)], alpha = alpha)
#   # res2 = ifelse(res2_A$decision == "NO" |  res2_B$decision == "NO", "NO", "LDM")
#   res2 = Test_LDM_Sen_stage2(res1, alpha=alpha)
#
#   return( c(res1,res2 ))
# }

