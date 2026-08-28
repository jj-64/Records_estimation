## Wrapper -------------------
#' Record Statistics for Multiple Models
#'
#' Compute the distribution, expected number of records or the variance of the number
#' of records under a variety of stochastic models.
#' @details For model-specific functions and formulas, look at
#' \code{\link{rec_count_mean_iid}}
#' \code{\link{rec_count_var_iid}}
#' \code{\link{rec_count_dist_iid}}
#' \code{\link{rec_count_mean_DTRW}}
#' \code{\link{rec_count_var_DTRW}}
#' \code{\link{rec_count_dist_DTRW}}
#' \code{\link{rec_count_mean_LDM}}
#' \code{\link{rec_count_var_LDM}}
#' \code{\link{rec_count_dist_LDM}}
#' \code{\link{rec_count_mean_YNM}}
#' \code{\link{rec_count_var_YNM}}
#' \code{\link{rec_count_dist_YNM}}
#' @param model Character string specifying the model. Available models:
#'   \code{"iid"}, \code{"DTRW"}, \code{"LDM"}, \code{"YNM"}.
#' @param stat Either \code{"mean"} or \code{"var"} or \code{"dist"}.
#' @param T Integer. Length of the time series.
#' @param m Integer. Number < T for which we compute the probability distribution
#' @param ... Additional arguments passed to model-specific functions.
#'
#' @return A numeric value, representing either the expected number
#'   of records or its variance or distribution.
#'
#' @export
#'
#' @examples
#' rec_count_stats("iid", stat = "mean", T = 100)
#' # [1] 5.187378
#' rec_count_stats("YNM", stat = c("mean"), T = 200, gamma = 1.1)
#' # [1] 20.99712
#' rec_count_stats("YNM", stat = "var", T = 200, gamma = 1.1)
#' # [1] 17.63478
#' rec_count_stats("LDM", stat="dist", T=100, m=7, theta=0.4, scale=1)
#' # [1] 0.00000000004169792
#' rec_count_stats("LDM", stat="dist", T=25, m=c(1:5), theta=0.4, scale=1)
#' # [1] 0.00002232982 0.00030729569 0.00201438186 0.00837570638 0.02481251456
#' rec_count_stats("iid", stat="dist", T=50, m=4)
#' rec_count_stats("DTRW", stat="dist", T=200, m=10, approx=TRUE)
rec_count_stats <- function(model,
                         stat = c("mean", "var", "dist"),
                         T,
                         m = NULL,
                         ...){

  stat <- match.arg(stat)

  # check model exists
  if (is.null(record_registry[[model]]))
    stop("Unknown model: ", model)

  entry <- record_registry[[model]]

  # check required args for model
  passed <- list(...)
  missing <- setdiff(entry$required_args, names(passed))
  if (length(missing) > 0) {
    stop("Model '", model, "' requires argument(s): ",
         paste(missing, collapse=", "))
  }

  # dispatch based on stat
  if (stat == "mean") {
    return(entry$mean_fun(T = T, ...))
  }

  if (stat == "var") {
    return(entry$var_fun(T = T, ...))
  }

  if (stat == "dist") {
    if (is.null(m)) stop("Argument 'm' must be provided for stat='dist'.")
    return(entry$dist_fun(m = m, T = T, ...))
  }

  stop("Unknown stat: ", stat)
}

## iid -------------------
#'Exact Expected number of records in Classical Model
#'
#' The average expected number of records in an i.i.d process. It is distribution-free,
#' i.e. independent from the process underlying distribution. It only depends from the length of the series.
#' @param T numeric, the length of the series
#' @param approx logical, if approximated value. default = FALSE
#'
#' @returns a single value of the expected number of records
#' @export
#' @details It is \eqn{E(N_T) = \sum_{1}^{T}\frac{1}{t}} where \eqn{T} is the length of the series
#' If approximated, as the sum of \eqn{1/t} can be approximated by \eqn{log(T) + \omega} where \eqn{\omega = 0.577215} is euler mascheroni constant
#'
#' @examples
#' rec_count_mean_iid(T=25)
#' # 3.815958
#' rec_count_mean_iid(T=25, approx = TRUE)
#' # 3.796091
rec_count_mean_iid = function(T, approx = FALSE){
  if (approx){
    return(log(T) + 0.57721566490153)
  } else {
  return(sum( 1/(1:T) )) }
}


#Variance of number of records
#' Variance of number of records in Classical Model
#'
#'The variance of number of records in classical model.
#'@details It is computed as
#'\eqn{V(N_T) = \sum_{1}^{T} 1/t - \sum_{1}^{T} 1/t^2}
#'where T is the length of series
#' @param T integer. length of the time series.
#' @returns a value of the variance of number of records
#' @export
#'
#' @examples rec_count_var_iid(T=25)
#' # [1] 2.210235
#'
#' rec_count_mean_iid(T=25)
#' # [1] 3.815958
#'
#' # For a series of length 25 and following a classical model, we expect to
#' # observe around 3.81 records with a variance of 2.21
#' rec_count_dist_iid(m=1,T=25)
rec_count_var_iid = function(T){
  sum( 1/(1:T) ) -  sum( 1/(1:T)^2 )
}

# ## Fucntion needed to compute the distribution of number of records

#' Stirling function of the first kind in Classical
#'
#' Stirling function of the first kind in the case of classical process.
#' It is used later to compute the distribution of number of records
#'
#'@details
#'Each item is computed as
#'\deqn{\prod_{t=0}^{n-1} (s-t) = \sum_{k=0}^{n} s^k S(n,k)}
#'
#' @param n positive integer
#' @param k integer in 0:n.
#' @returns number of S(n,k)
#' @export
Stirling_first_kind = function(n, k){
  as.numeric ( abs(gmp::Stirling1(n=n, k=k)) )
}

#' Distribution of number of records in Classical Model
#'
#' Exact distribution of the number of records in classical model. The function
#' computes the probability of observing \eqn{m} records in a process of length
#' \eqn{T}.
#' @details
#' Computes the exact probability of observing \eqn{m} records in a process of
#' length \eqn{T},
#' given by:
#' \deqn{ P(N_T = m) = \frac{s(T, m)}{T!} }
#' where \eqn{s(T, m)} is the unsigned Stirling number of the first kind.
#'
#' @param m Integer. The number of records whose probability is to be computed.
#' @param T Integer. The length of the series.
#' @param s Optional. The Stirling number of the first kind, \eqn{s(T,m)}.
#' If not provided, it will be computed internally.
#'
#' @return A probability value (numeric) less than or equal to 1.
#'
#' @examples
#' rec_count_dist_iid(m = 1, T = 25)
#' # [1] 0.04
#'
#' rec_count_dist_iid(m = 2, T = 25)
#' # [1] 0.1510383
#'
#' rec_count_dist_iid(m = 2, T = 25, s = Stirling_first_kind(n = 25, k = 2))
#' # [1] 0.1510383
#'
#' @references
#' Nevzorov, V. B. (2001). *Records: Mathematical Theory*. American Mathematical
#'  Society.
#' \doi{10.1090/surv/172}
#'
#' @export
rec_count_dist_iid = function(m,T,s=NA){
  if(is.na(s)== TRUE) s= as.numeric ( abs(gmp::Stirling1(n=T, k=m)) )
      #Stirling_first_kind(n=T, k=m)
  s / factorial(T)
}

## DTRW ----------------------
#' Exact Expected number of records in DTRW Model
#'
#' The average exact expected number of records in a DTRW process.
#' It is distribution-free, i.e. independent from the process underlying distribution.
#' It only depends from the length of the series.
#' @details It is computed as \eqn{E(N_T) = (2T +1) \times 2^{-2T}
#' \times C_{T}^{2T}} where \eqn{T} is the length of the process
#' If approximated, it is computed as \eqn{E(N_T) = \sqrt(4T/\pi)} where \eqn{T}
#'  is the length of the process
#'
#' @param T numeric, the length of the series
#' @param approx logical, if approximated value. default = FALSE
#'
#' @returns a single value of the expected number of records
#' @export
#'
#' @examples
#' rec_count_mean_DTRW(T=25)
#' # [1] 5.726034
#'
#' rec_count_mean_DTRW(T=25, approx = TRUE)
#' # [1] 5.641896
rec_count_mean_DTRW <- function(T, approx = FALSE) {
  if(approx) {
    return(sqrt(4*T/pi))
    } else {
  return((2*T+1)*choose(2*T,T)*2^(-2*T)) }
}

## Variance of expected number of Records DTRW
#' Variance of number of records in DTRW Model
#'
#'The approximated variance of number of records in DTRW model.
#' @param T numeric, the length of the series
#' @param approx logical, if approximated value. default = FALSE
#'
#'@details Exact variance is computed as
#'\eqn{V(N_T)= 2T+2-E(N_T)-E(N_T)^2}
#'where \eqn{T} is the length of series and \eqn{E(N_T)} is the expected number
#'of records in DTRW process
#'
#' Approximated variance is computed as
#'\eqn{V(N_T)_approx = 2T \times (1-\frac{2}{\pi})}
#'
#' @returns a value of the variance of number of records
#' @export
#'
#' @examples
#' rec_count_var_DTRW(T=25, approx = TRUE)
#' # [1] 18.16901
#'
#' rec_count_mean_DTRW(T=25, approx = TRUE)
#' # [1] 5.641896
#'
#' rec_count_var_DTRW(T=25, approx = FALSE)
#' # [1] 13.4865
#'
#' rec_count_mean_DTRW(T=25, approx = FALSE)
#' # [1] 5.726034
#'
#' # For a series of length 25 and following a DTRW model, we expect to observe
#' # around 5.62 records with an approximated variance of 18.169
rec_count_var_DTRW = function(T, approx = FALSE){
  if (approx){
  return(2*(1-2/pi)*T)
  } else {
    m = rec_count_mean_DTRW(T, approx= FALSE)
    v=2*T+2-m-m^2
    return(v)
  }
}


#' Approximated Distribution of number of records in DTRW Model
#'
#'Exact and approximate distribution of the number of records in DTRW model.
#' The function computes the probability of observing \eqn{m} records in a
#' process of length \eqn{T}.
#' #' @details
#' The Exact probability of observing \eqn{m} records in a process of length
#' \eqn{T} is given by:
#' \deqn{ P(N_T = m) = 2^{-2T+m-1} \times C_{T}^{2T-m+1} }
#' where  \eqn{T} is the length of the process
#'
#' The approximated probability of observing \eqn{m} records in a process of length \eqn{T} is given by:
#' \deqn{ P(N_T = m) = \frac{e^{\frac{-m^2}{4T}}}{\sqrt(\pi T)} }
#' where  \eqn{T} is the length of the process
#'
#' @param m the number of records we are computing its probabolity (integer)
#' @param T length of the series (integer)
#' @param approx logical, if approximated value. Default = FALSE
#'
#' @returns a probability less than one
#' @export
#'
#' @examples
#' rec_count_dist_DTRW(m=1,T=25, approx = TRUE)
#' # [1]  0.1117152
#'
#' rec_count_dist_DTRW(m=1,T=25, approx = TRUE)
#' # [1] 0.1122752
#'
#' rec_count_dist_DTRW(m=1,T=25)
#' # [1] 0.1122752
#'
#' rec_count_dist_DTRW(m=5,T=25)
#' # [1] 0.09867345
rec_count_dist_DTRW=function(m,T, approx = FALSE){
  if ( approx){
  return (exp(-m^2/(4*T))/sqrt(pi*T))
  } else {
    choose(2*T-m+1,T)*2^(-2*T+m-1)
  }
}


#' Survival Probability
#'
#'Probability of \eqn{{X_1, ..., X_n}} to stay below zero
#' @param n integer, time
#'
#' @returns a probability
#' @export
#'
#' @examples
#' Survival(10)
#' # [1] 0.1761971
Survival = function(n){ 2^(-2*n) * choose(2*n,n)}

#' First Pass probability
#'
#'Probability of \eqn{{X_1, ..., X_n-1}} to stay below zero and \eqn{X_n} to pass zero for the first time
#' @param n integer, time
#'
#' @returns a probability
#' @export
#'
#' @examples
#' FirstPass(10)
#' # [1] 0.009273529
FirstPass=function(n){Survival(n-1)-Survival(n)}

## LDM --------------

### Expected number of records LDM
#'  Expected number of records in Linear Drift Model (LDM)
#'
#' The average expected number of records in a LDM process with a Gumbel underlying distribution.
#' @details
#' An LDM process is defined as \eqn{X_t = Y_t + \theta t },
#' where \eqn{Y_t} are independent and identically distributed random variable.
#'
#' An explicit formula exists when \eqn{Y_t} follow a \eqn{Gumbel(\alpha, \beta)} distribution. It is computed as:
#'
#'  \deqn{E(N_T) = \sum_{t=1}^{T} \frac{1-e^{\theta / \beta}}{1-e^{\theta t / \beta}}}
#'
#' For other increment distributions, the variance is estimated by simulation.
#'
#' @param T Integer. The length of the series.
#' @param theta Numeric. The linear drift coefficient \eqn{\theta > 0}.
#' @param dist Character, distribution name. One of:
#'   "beta", "gumbel", "weibull", "frechet", "exp", "pareto", "norm", "exp", "pareto", "uniform".
#' @param ... Additional parameters specific to the chosen distribution:
#'   \describe{
#'     \item{beta}{`shape1`, `shape2`}
#'     \item{gumbel}{`loc`, `scale`}
#'     \item{weibull}{`shape`, `scale`}
#'     \item{frechet}{`shape`, `scale`}
#'     \item{norm}{`loc`, `sd`}
#'     \item{exp}{`rate`}
#'     \item{pareto}{`scale`, `shape`}
#'     \item{unifom}{`min`, `max`}
#'   }
#' @param n_sim Numeric. Number of Monte Carlo simulations (default 1000)
#' @returns a single value of the expected number of records
#' @export
#'
#' @examples
#' rec_count_mean_LDM(T=25, theta=0.5, dist="gumbel", n_sim=100,location=0, scale=1)
#' # [1] 10.93343
#'
#' rec_count_mean_LDM(T=25, theta=0.5, dist="gumbel", n_sim=100,location=0, scale=2)
#' # [1] 7.320699
#'
#' rec_count_mean_LDM(T=25, theta=0.5, dist="norm", n_sim=100,mean=0, sd=1)
#' # [1] 13.18
rec_count_mean_LDM <- function(T, theta, dist = c("beta", "gumbel", "weibull", "frechet", "norm", "exp", "pareto", "uniform"), n_sim = 1000, ...) {
  dist <- match.arg(dist)   # enforce valid choice
  args <- list(...)

  if (dist == "gumbel") {
    ## Explicit formula: sum of record probabilities
    s <- rec_rate_LDM(t = 1:T, theta=theta, loc = 0, scale = args$scale)
    return(sum(s))

  } else {
    ## Simulation-based estimate
    recs <- numeric(n_sim)
    for (i in 1:n_sim) {
      X <- LDM_series(T = T, theta = theta, dist = dist, ...)
      recs[i] <- rec_count(X)
    }
    return(mean(recs))
  }
}



### Variance of the number of Records LDM
#'  Variance of the number of records in Linear Drift Model (LDM)
#'
#' The variance of the number of records in a LDM process with a Gumbel underlying distribution
#' or estimated by Monte Carlo simulation otherwise.
#'
#' @details
#' An LDM process is defined as \eqn{X_t = Y_t + \theta t },
#' where \eqn{Y_t} are independent and identically distributed random variable.
#'
#' the variance is computed as:
#' \deqn{Var(N_T) = \sum_{t=1}^T P_t(1 - P_t),}
#' where \eqn{P_t} is the probability of a record at time \eqn{t}.
#'
#' An explicit formula exists in the case of Gumbel \eqn{G(\alpha, \beta)} distribution It is computed as:
#'
#'  \deqn{V(N_T) = \sum_{t=1}^{T} \frac{1-e^{\theta / \beta}}{1-e^{\theta t / \beta}} - \sum_{t=1}^{T} (\frac{1-e^{\theta / \beta}}{1-e^{\theta t / \beta}})^2}

#' For other increment distributions, the variance is estimated by simulation.
#'
#' @inheritParams rec_count_mean_LDM
#' @returns a single value: the variance of the number of records
#' @export
#'
#' @examples
#' rec_count_var_LDM(T=25, theta=0.5, dist="gumbel", n_sim=100, loc=0, scale=1)
#' # [1] 5.761166
#'
#' rec_count_var_LDM(T=25, theta=0.5, dist="gumbel", n_sim=100, loc=0, scale=2)
#' # [1] 4.509757
#'
#' rec_count_var_LDM(T=25, theta=0.5, dist="norm", n_sim=100, mean=0, sd=1)
#' # [1] 4.492424
rec_count_var_LDM <- function(T,
                    theta,
                    dist = c("beta", "gumbel", "weibull", "frechet", "norm", "exp", "pareto", "uniform"),
                    n_sim = 1000,
                    ...) {
  dist <- match.arg(dist)
  args <- list(...)

  if (dist == "gumbel") {
    ## Explicit formula
    s <- rec_rate_LDM(t = 1:T, theta = theta, loc = args$location, scale = args$scale)
    return(sum(s * (1 - s)))

  } else {
    ## Simulation-based estimate
    recs <- numeric(n_sim)
    for (i in 1:n_sim) {
      X <- LDM_series(T = T, theta = theta, dist = dist, ...)
      recs[i] <- rec_count(X)
    }
    return(var(recs))
  }
}


## Apply Liapanuv conditions
# CLT_LDM = function(X){
#   theta = theta_estm2B(X)
#   t=length(X)
#   N = 1  ## first trivial record
#   for(i in 2:t) {N[i] = rec_count(X[1:i])}  ## observed number of records series
#   E=1
#   for(i in 2:t) {E[i] = rec_count_mean_LDM(T=i, theta=theta)}  ## expected number of records
#   sigm=rec_count_var_LDM(T=1,theta=theta)
#   for(i in 2:t) {sigm[i] =rec_count_var_LDM(T=i,theta=theta)} ## sum of variance
#   s=sqrt(cumsum(sigm))
#
#   Z = (N-E)/s
#   return(Z)
# }
#
# ## Fucntion needed to compute the distribution of number of records


#' Stirling function of the second kind in LDM
#'
#' Stirling function of the second kind in the case of LDM process. It is used later to compute the distribution of number of records
#'
#'@details
#'Each item of the matrix is computed as
#'\deqn{\prod_{t=1}^{T} u_t}
#'where \eqn{u_t = e^{-\theta/\beta} \times \frac{(1-e^{-\theta t/\beta})}{(1-e^{-\theta/\beta})}}

#' \eqn{T} is the length of the process, \eqn{\theta} is the LDM slope, \eqn{\beta} is the scale parameter of the \eqn{Gumbel} distribution
#' @param T the length of the series
#' @param theta the LDM parameter, slope \eqn{\theta > 0}
#' @param scale (optional, 1 by default) the scale or \eqn{\beta} parameter of the \eqn{Gumbel} distribution

#'
#' @returns a lower triangular \eqn{T \times T} matrix of all combinations of Stirling number.
#' @export
#'
#' @examples
#' Stirling_2nd_LDM(T=5, theta=0.5)
#' #       [,1]     [,2]     [,3]     [,4] [,5]
#' # [1,] 1.0000000 0.000000 0.000000 0.000000    0
#' # [2,] 0.6065307 1.000000 0.000000 0.000000    0
#' # [3,] 0.5910096 1.580941 1.000000 0.000000    0
#' # [4,] 0.7077578 2.484250 2.778481 1.000000    0
#' # [5,] 0.9433531 4.018954 6.187619 4.111357    1
Stirling_2nd_LDM = function(T,theta,scale=1){  ## compute stirling number of second kind

   u_t_LDM=function(t,theta,scale=1){
     exp(-theta/scale) * (1-exp(-theta*t/scale)) / (1-exp(-theta/scale)) }

  u = u_t_LDM(t=1:T,theta=theta,scale=scale)

  s=matrix(0,nrow=T,ncol=T)  ##S(T,m)= 0 for m>T upper diagonal
  s[1,1]=1 ## first one is always 1

  for(j in 2:T){			## rows are T, columns are m
    s[j,1]=u[j-1]*s[j-1,1] ## since s(j,0)=0
    s[j,j]=1  ## diagonal is 1

    if(j>2){
      for(k in 2 : (j-1) ){
        s[j,k]= s[j-1,k-1] + u[j-1] * s[j-1,k]
      }
    }
  }
  return(s)
}

#' Distribution of number of records in LDM process
#'
#'Distribution of the number of records in LDM model with \eqn{Gumbel} underlying distribution. The function computes the probability of observing \eqn{m} records in a process of length \eqn{T}.
#' @details
#' An LDM process is defined as \eqn{X_t = Y_t + \theta t },
#' where \eqn{Y_t} are independent and identically distributed random variable.
#' The properties are studied when \eqn{Y_t} follow a \eqn{Gumbel(\alpha, \beta)} distribution.
#'
#' It is computed as \eqn{P[N_T=m] = e^{\theta \times T / \beta} * s[T,m] / \prod_{t=1}^{T} u_t }
#'
#' where \eqn{u_t = e^{-\theta/\beta} \times \frac{(1-e^{-\theta t/\beta})}{(1-e^{-\theta/\beta})}}
#'
#' \eqn{T} is the length of the process, \eqn{\theta} is the LDM slope, \eqn{\beta} is the scale parameter of the \eqn{Gumbel} distribution
#' @param m the observed number of record
#' @param T the length of the series
#' @param theta the LDM parameter, slope \eqn{\theta > 0}
#' @param scale (optional, 1 by default) the scale or \eqn{\beta} parameter of the \eqn{Gumbel} distribution
#' @param s (NA by default) the Stirling matrix of the seconf kind
#' @returns a probability less than one
#' @export
#'
#' @examples
#' rec_count_dist_LDM(m=5,T=25, theta=0.5)
#' # 0.006915892
#'
#' rec_count_dist_LDM(m=5,T=25, theta=0.5, s = Stirling_2nd_LDM(T=25,theta=0.5,scale=1))
#' # 0.006915892
rec_count_dist_LDM = function(m,T,theta,scale=1,s=NA){  ## number of m, T, theta and Stirling matrix

  u_t_LDM=function(t,theta,scale=1){
    exp(-theta/scale) * (1-exp(-theta*t/scale)) / (1-exp(-theta/scale)) }

  ## compute stirling matrix
  if(is.na(s)[1] == TRUE) {
    s=Stirling_2nd_LDM(T=T,theta=theta,scale=scale)}

  p=(prod(u_t_LDM(t=1:T,theta=theta,scale=scale)))

  return(exp(-theta*T/scale) * s[T,m]/ p)
}

## YNM-Nevzorov -----------------

### Expected number of records in YANG-Nevzorov process
#' Expected number of records in YANG-Nevzorov process
#'
#'The expected number of records in a YNM-Nevzorov process. It is distribution free, i.e. does not depend on the underlying distribution of the process.
#'@details
#' A YNM-Nevzorov process can be expressed as
#' \deqn{X_t = max\{Y_1, \cdots, Y_{\gamma^t}\}}
#' where \eqn{t} is the time, \eqn{Y_t} are idenpendent and identically distributed random vector.
#'
#'In fact, \eqn{F_{X_t}(x) = \{F_{Y_t}(x)\}^{\gamma^t}}, where \eqn{F(.)} is the cumulative distribution function
#'
#'The expected number is computed as
#'\deqn{E(N_T) = \sum_{t=1}^{T} \frac{\gamma^t (\gamma -1)}{\gamma  (\gamma^t -1)}}
#' @param T the length of the series
#' @param gamma The power of the process, \eqn{\gamma >= 1}
#'
#' @returns a value for the expected number of records
#' @export
#'
#' @examples
#' rec_count_mean_YNM (T=25, gamma=1.1)
#' # [1] 5.000207
rec_count_mean_YNM = function(T,gamma){
  s=0
  for(k in 1:T){
    #s[k] = (gamma^k * (gamma-1)) / (gamma* (gamma^k -1))
    s[k] = rec_rate_YNM(gamma,k) }

  return(sum(s))
}

#'  Variance of number of records in YANG-Nevzorov process
#'
#' The variance of number of records in a YANG-Nevzorov process.
#'@details
#' A YNM-Nevzorov process can be expressed as
#' \deqn{X_t = max\{Y_1, \cdots, Y_{\gamma^t}\}}
#' where \eqn{t} is the time, \eqn{Y_t} are idenpendent and identically distributed random vector.
#'
#'In fact, \eqn{F_{X_t}(x) = \{F_{Y_t}(x)\}^{\gamma^t}}, where \eqn{F(.)} is the cumulative distribution function
#'
#'The variance of number of records is computed as
#'\deqn{V(N_T) = \sum_{t=1}^{T} P_t - \sum_{t=1}^{T} (P_t)^2}
#'
#'where \eqn{P_t =  \frac{\gamma^t (\gamma -1)}{\gamma  (\gamma^t -1)}}

#' @param T the length of the series
#' @param gamma The power of the process, \eqn{\gamma >= 1}
#'
#' @returns a single value of the variance of number of records
#' @export
#'
#' @examples
#' rec_count_var_YNM(T=25, gamma=1.1)
#' # [1] 3.10049
rec_count_var_YNM = function(T, gamma){
  s=0; s2=0
  for(k in 1:T){
    s[k] = rec_rate_YNM(gamma,k)
    s2[k] =s[k]^2
  }
   return(sum(s)-sum(s2))
}

## helper function for computing
u_t_YNM = function(t, gamma) {(1-gamma^t)/(gamma^t * (1- gamma))}

# 'Weighted Stirling Numbers of the Second Kind (YANG-Nevzorov version)
#'
#' Stirling function of the second kind in the case of YANG-Nevzorov (YN) process.
#'
#' It is used later to compute the distribution of number of records
#'
#' Computes a triangular table of weighted Stirling numbers of the second kind
#' adapted for YNM's model with parameter \eqn{\gamma}.
#'
#'
#'@details
#'Each item of the matrix is computed as
#'\deqn{\prod_{t=1}^{T} u_t}
#'where \eqn{u_t = \frac{1-\gamma^t}{\gamma^t \times (1-\gamma)}}

#' \eqn{T} is the length of the process, \eqn{\gamma} is the YN power
#'
#' The recursion is:
#' \deqn{S(j,k) = S(j-1, k-1) + u_{j-1} S(j-1,k),}
#' with boundary conditions:
#' - \eqn{S(1,1) = 1}
#' - \eqn{S(j,1) = u_{j-1} S(j-1,1)}
#' - \eqn{S(j,j) = 1}
#'
#'
#' @param T the length of the series
#' @param gamma The power of the process, \eqn{\gamma >= 1}
#'
#' @returns  A \eqn{T \times T} lower triangular matrix of all combinations of Stirling number.
#' @export
#'
#' @examples
#' Stirling_2nd_YNM(T=5, gamma=1.1)
#' #          [,1]      [,2]     [,3]     [,4]   [,5]
#' # [1,]  1.0000000  0.000000  0.00000 0.000000    0
#' # [2,]  0.9090909  1.000000  0.00000 0.000000    0
#' # [3,]  1.5777611  2.644628  1.00000 0.000000    0
#' # [4,]  3.9236583  8.154560  5.13148 1.000000    0
#' # [5,] 12.4374688 29.772515 24.42066 8.301346    1
Stirling_2nd_YNM = function(T,gamma){  ## compute stirling number of second kind

  # Precompute u_t for efficiency
  t = 1:T
  u <- (1-gamma^t)/(gamma^t * (1- gamma))  #u_t_YNM

  # Initialize triangular matrix
  s <- matrix(0, nrow = T, ncol = T) ##S(T,m)= 0 for m>T upper diagonal
  s[1, 1] <- 1  ## first one is always 1

  # Fill table
  for (j in 2:T) {  ## rows are T, columns are m
    s[j, 1] <- u[j - 1] * s[j - 1, 1] ## since s(j,0)=0
    s[j, j] <- 1  ## diagonal is 1

    if (j > 2) {
      for (k in 2:(j - 1)) {
        s[j, k] <- s[j - 1, k - 1] + u[j - 1] * s[j - 1, k]
      }
    }
  }
  return(s)
}



#' Distribution of number of records in YANG-Nevzorov (YN) process
#'
#'Distribution of the number of records YANG-Nevzorov (YN) process. The function computes the probability of observing \eqn{m} records in a process of length \eqn{T}.
#' @details
#' A YNM-Nevzorov process can be expressed as
#' \deqn{X_t = max\{Y_1, \cdots, Y_{\gamma^t}\}}
#' where \eqn{t} is the time, \eqn{Y_t} are idenpendent and identically distributed random vector.
#'
#'In fact, \eqn{F_{X_t}(x) = \{F_{Y_t}(x)\}^{\gamma^t}}, where \eqn{F(.)} is the cumulative distribution function
#'
#' It is computed as \eqn{P[N_T=m] = \gamma^{-T} * s[T,m] / \prod_{t=1}^{T} u_t }
#'
#'where \eqn{u_t = \frac{1-\gamma^t}{\gamma^t (1-\gamma)}}
#'
#' \eqn{T} is the length of the process, \eqn{\theta} is the LDM slope, \eqn{\beta} is the scale parameter of the \eqn{Gumbel} distribution
#' @param m the observed number of record
#' @param T the length of the series
#' @param gamma The power of the process, \eqn{\gamma >= 1}
#' @param s Optional. Precomputed Stirling matrix from
#'   \code{\link{Stirling_2nd_YNM}} (default = NA so it is recomputed internally).
#'
#' @returns Probability value \eqn{P(N_T = m)} less than one
#' @export
#'
#' @examples
#' rec_count_dist_YNM(m=5,T=25, gamma=1.1)
#'  # 0.2223667
#' rec_count_dist_YNM(m=5,T=25, gamma=1.1, s = Stirling_2nd_YNM(T=25,gamma=1.1))
#'  # 0.2223667
rec_count_dist_YNM <- function(m, T, gamma, s = NULL) {
  if (is.null(s)) {
    s <- Stirling_2nd_YNM(T = T, gamma = gamma)
  }

  # Precompute product of u_t
  t = 1:T
  u <- (1-gamma^t)/(gamma^t * (1- gamma))  #u_t_YNM
  p <- prod(u)

  return(s[T, m] / ((gamma^T) * p))
}

