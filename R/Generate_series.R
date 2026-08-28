## LDM -----------------------------
#' Generate a Linear Drift Model (LDM) Series with different Noise distribution
#'
#' Generates an LDM series where the noise follows a user defined distribution.
#'
#' @details
#' The series is generated as:
#' \deqn{x_t = \theta t + y_t, \quad y_t \sim \text{F}(\text{.}, \beta)}
#' where:
#' - \eqn{\theta} is the linear drift.
#' - \eqn{y_t} are i.i.d. distributed random variables.
#' - \eqn{\beta} are specific parameters of the distribution F of \eqn{Y} .
#'
#' @param T Integer. The length of the series.
#' @param theta Numeric. The linear drift coefficient \eqn{\theta > 0}.
#' @param dist Character, distribution name. One of:
#'   "beta", "gumbel", "weibull", "frechet", "exp", "pareto", "norm", "exp",
#'    "pareto", "uniform".
#' @param ... Additional parameters specific to the chosen distribution:
#'   \describe{
#'     \item{beta}{`shape1`, `shape2`}
#'     \item{gumbel}{`loc`, `scale`}
#'     \item{weibull}{`shape`, `scale`}
#'     \item{frechet}{`shape`, `scale`}
#'     \item{norm}{`mean`, `sd`}
#'     \item{exp}{`rate`}
#'     \item{pareto}{`scale`, `shape`}
#'     \item{unifom}{`min`, `max`}
#'   }
#' @return A numeric vector representing the LDM series.
#' @export
#' @examples
#' #Gumbel with loc, scale
#'LDM_series(100, theta = 0.1, dist = "gumbel", loc = 0, scale = 2)
#'
#' #Weibull
#'LDM_series(100, theta = 0.05, dist = "weibull", shape = 2, scale = 1)
#'
#' #Beta
#'LDM_series(100, theta = 0.2, dist = "beta", shape1 = 2, shape2 = 5)
#'
#' #Normal
#'LDM_series(100, theta = 0.1, dist = "norm", mean = 0, sd = 1)
LDM_series <- function(T, theta, dist = c("beta", "gumbel", "weibull", "frechet",
                                    "norm", "exp", "pareto", "uniform"), ...) {
  dist <- match.arg(dist)   # enforce valid choice
  args <- list(...)

  # Simulate noise y according to chosen distribution
  y <- switch(
    dist,

    beta = {
      shape1 <- args$shape1 %||% 1
      shape2 <- args$shape2 %||% 1
      if (shape1 <= 0 | shape2 <= 0) stop("Enter positive values for shape1 and shape2")
      rbeta(T, shape1, shape2)
    },

    gumbel = {
      loc   <- args$loc   %||% 0
      scale <- args$scale %||% 1
      if (scale <= 0) stop("Enter a positive value for scale")
      VGAM::rgumbel(T, location = loc, scale = scale)
    },

    weibull = {
      shape <- args$shape %||% 1
      scale <- args$scale %||% 1
      if (scale <= 0 | shape <= 0) stop("Enter positive values for shape and scale")
      rweibull(T, shape = shape, scale = scale)
    },

    frechet = {
      scale <- args$scale %||% 1
      shape <- args$shape %||% 2
      if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")
      VGAM::rfrechet(T, loc = 0, scale = scale, shape = shape)
    },

    norm = {
      mean <- args$mean %||% 0
      sd   <- args$sd   %||% 1
      if (sd <= 0) stop("Enter positive value for sd")
      rnorm(T, mean = mean, sd = sd)
    },

    exp = {
      rate <- args$rate %||% 1
      if (rate <= 0) stop("Enter positive value for rate")
      rexp(T, rate=rate)
    },

    pareto = {
      scale <- args$scale %||% 1
      shape <- args$shape %||% 2
      if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")
      VGAM::rpareto(T, scale, shape)
    },

    uniform ={
      min <- args$min %||% -1
      max <- args$max %||% 1
      runif(T, min, max)
    }

  )

  # Add deterministic linear drift
  x <- theta * (1:T) + y
  return(x)
}


# #' Generate a Linear Drift Model (LDM) Series with Beta Noise
# #'
# #' Generates an LDM series where the noise follows a Beta distribution.
# #'
# #' @details
# #' The series is generated as:
# #' \deqn{x_t = \theta t + y_t, \quad y_t \sim \text{Beta}(\text{shape1}, \text{shape2})}
# #' where:
# #' - \eqn{\theta} is the linear drift.
# #' - \eqn{y_t} are i.i.d. Beta-distributed random variables.
# #'
# #' @param T Integer. The length of the series.
# #' @param theta Numeric. The linear drift coefficient \eqn{\theta > 0}.
# #' @param shape1 Positive numeric. The first shape parameter of the Beta distribution.
# #' @param shape2 Positive numeric. The second shape parameter of the Beta distribution.
# #' @return A numeric vector representing the LDM series.
# #' @examples
# #' LDM_series_Beta(100, theta = 0.5, shape1 = 2, shape2 = 5)
# LDM_series_Beta <- function(T, theta, shape1 = 1, shape2 = 1) {
# if (shape1 <= 0 | shape2 <= 0) stop("Enter positive values for shape parameters")
# y <- rbeta(T, shape1, shape2)
# x <- theta * (1:T) + y
# return(x)
# }


# #' Generate an LDM Series with Gumbel Noise
# #'
# #' @details
# #' \deqn{x_t = \theta t + y_t, \quad y_t \sim \text{Gumbel}(\text{loc}, \text{scale})}
# #'
# #' @inheritParams LDM_series_Beta
# #' @param loc Numeric. The location parameter of the Gumbel distribution.
# #' @param scale Positive numeric. The scale parameter of the Gumbel distribution.
# #' @return A numeric vector representing the LDM series.
# LDM_series_Gumbel <- function(T, theta, loc = 0, scale = 1) {
# if (scale <= 0) stop("Enter a positive value for scale")
# y <- VGAM::rgumbel(T, loc, scale)
# x <- theta * (1:T) + y
# return(x)
# }

# #' Generate an LDM Series with Weibull Noise
# #'
# #' @details
# #' \deqn{x_t = \theta t + y_t, \quad y_t \sim \text{Weibull}(\text{scale}, \text{shape})}
# #'
# #' @inheritParams LDM_series_Beta
# #' @param scale Positive numeric. The scale parameter of the Weibull distribution.
# #' @param shape Positive numeric. The shape parameter of the Weibull distribution.
# #' @return A numeric vector representing the LDM series.
# LDM_series_Weibull <- function(T, theta, shape = 1,scale = 1) {
# if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")
# y <- rweibull(T, shape=shape, scale=scale)
# x <- theta * (1:T) + y
# return(x)
# }

# #' Generate an LDM Series with Frechet Noise
# #'
# #' @details
# #' \deqn{x_t = \theta t + y_t, \quad y_t \sim \text{Frechet}(\text{scale}, \text{shape})}
# #'
# #' @inheritParams LDM_series_Beta
# #' @param scale Positive numeric. The scale parameter of the Frechet distribution.
# #' @param shape Positive numeric. The shape parameter of the Frechet distribution.
# #' @return A numeric vector representing the LDM series.
# LDM_series_Frechet <- function(T, theta, scale = 1, shape = 2) {
# if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")
# y <- VGAM::rfrechet(T, loc = 0, scale, shape)
# x <- theta * (1:T) + y
# return(x)
# }

# #' Generate an LDM Series with Exponential Noise
# #'
# #' @details
# #' \deqn{x_t = \theta t + y_t, \quad y_t \sim \text{Exp}(\text{rate})}
# #'
# #' @inheritParams LDM_series_Beta
# #' @param rate Positive numeric. The rate parameter of the Exponential distribution.
# #' @return A numeric vector representing the LDM series.
# LDM_series_Exp <- function(T, theta, rate = 1) {
# if (rate <= 0) stop("Enter a positive value for scale:1/ rate")
# y <- rexp(T, rate=rate)
# x <- theta * (1:T) + y
# return(x)
# }

# #' Generate an LDM Series with Pareto Noise
# #'
# #' @details
# #' \deqn{x_t = \theta t + y_t, \quad y_t \sim \text{Pareto}(\text{scale}, \text{shape})}
# #'
# #' The top 10% of values are filtered out to reduce extreme values.
# #'
# #' @inheritParams LDM_series_Beta
# #' @inheritParams LDM_series_Weibull
# #' @return A numeric vector representing the LDM series.
# LDM_series_Pareto <- function(T, theta, scale = 0.5, shape = 4) {
# if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")
# y <- VGAM::rpareto(T * 1.5, scale, shape)
# y <- y[y <= quantile(y, 0.9)] # Remove top 10% of extreme values
# x <- theta * (1:T) + y[1:T]
# return(x)
# }

# #' Generate an LDM Series with Normal Noise
# #'
# #' @details
# #' \deqn{x_t = \theta t + y_t, \quad y_t \sim \mathcal{N}(\text{loc}, \text{sd})}
# #'
# #' @inheritParams LDM_series_Beta
# #' @param loc Numeric. The mean (location) of the normal distribution.
# #' @param sd Positive numeric. The standard deviation of the normal distribution.
# #' @return A numeric vector representing the LDM series.
# LDM_series_Norm <- function(T, theta, loc = 0, sd = 1) {
# if (sd <= 0) stop("Enter a positive value for standard deviation")
# y <- rnorm(T, loc, sd)
# x <- theta * (1:T) + y
# return(x)
# }

# #' Generate an LDM Series with Uniform Noise
# #'
# #' @details
# #' \deqn{x_t = \theta t + y_t, \quad y_t \sim \text{Uniform}(\text{min}, \text{max})}
# #'
# #' @inheritParams LDM_series_Beta
# #' @param min Numeric. The lower bound of the uniform distribution.
# #' @param max Numeric. The upper bound of the uniform distribution.
# #' @return A numeric vector representing the LDM series.
# LDM_series_Unif <- function(T, theta, min = -1, max = 1) {
# y <- runif(T, min, max)
# x <- theta * (1:T) + y
# return(x)
# }


## YNM --------------------
#' Generate a YNM Process
#'
#' Simulates a YNM process of length T under different underlying distributions.
#' The YNM process is defined by a parameter gamma > 0 controlling the block maxima.
#'
#' @param T Integer, length of the series.
#' @param gamma Numeric > 0, YNM process parameter.
#' @param dist Character, distribution name. One of:
#'   "beta", "gumbel", "weibull", "frechet", "exp", "pareto", "norm", "pareto_trunc".
#' @param ... Additional parameters specific to the chosen distribution:
#'   \describe{
#'     \item{beta}{`shape1`, `shape2`}
#'     \item{gumbel}{`loc`, `scale`}
#'     \item{weibull}{`shape`, `scale`}
#'     \item{frechet}{`shape`, `scale`}
#'     \item{exp}{`rate`}
#'     \item{pareto}{`scale`, `shape`}
#'     \item{norm}{`mean`, `sd`}
#'     \item{pareto_trunc}{`scale`, `shape`, `xmax`}
#'   }
#'
#' @return A numeric vector of length T, the simulated YNM process.
#' @examples
#' \dontrun{
#' YNM_series(100, gamma = 1.5, dist = "gumbel", loc = 0, scale = 1)
#' YNM_series(10, gamma = 2, dist = "beta", shape1 = 2, shape2 = 5)
#' YNM_series(100, gamma = 1.2, dist = "norm", loc = 0, sd = 1)
#' }
#' @export
YNM_series <- function(T, gamma, dist = c("beta", "gumbel", "weibull",
                                          "frechet", "exp", "pareto",
                                          "norm", "pareto_trunc"), ...) {
  dist <- match.arg(dist)
  args <- list(...)
  X <- numeric(T)

  for (i in 1:T) {
    m <- gamma^i
    u <- runif(1)

    X[i] <- switch(
      dist,

      beta = {
        shape1 <- args$shape1 %||% 1
        shape2 <- args$shape2 %||% 1
        if (shape1 <= 0 | shape2 <= 0) stop("Enter positive values for shape1 and shape2")
        max(rbeta(m, shape1, shape2))
      },

      gumbel = {
        loc   <- args$loc   %||% 0
        scale <- args$scale %||% 1
        if (scale <= 0) stop("Enter a positive value for scale")
        loc - scale * log(-log(u^(1/m)))
      },

      weibull = {
        shape <- args$shape %||% 1
        scale <- args$scale %||% 1
        if (scale <= 0 | shape <= 0) stop("Enter positive values for shape and scale")
        (-scale^shape * log(1 - u^(1/m)))^(1/shape)
      },

      frechet = {
        shape <- args$shape %||% 2
        scale <- args$scale %||% 1
        if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")
        scale * (-(1/m) * log(u))^(-1/shape)
      },

      exp = {
        rate <- args$rate %||% 1
        if (rate <= 0) stop("Enter positive value for rate")
        -log(1 - u^(1/m)) / rate
      },

      pareto = {
        scale <- args$scale %||% 0.5
        shape <- args$shape %||% 4
        if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")
        scale * (1 - u^(1/m))^(-1/shape)
      },

      norm = {
        loc <- args$mean %||% 0
        sd  <- args$sd  %||% 1
        if (sd <= 0) stop("Enter positive value for sd")
        qnorm(u^(1/m), mean = loc, sd = sd)
      },

      pareto_trunc = {
        scale <- args$scale
        shape <- args$shape
        xmax  <- args$xmax

        if (scale <= 0 || shape <= 0 || xmax <= scale)
          stop("Invalid parameters")

        # truncated Pareto inverse CDF for the maximum
        a <- (scale / xmax)^shape
        scale * (1 - u^(1/m) * (1 - a))^(-1/shape)
      }
    )
  }

  X <- X[is.finite(X)] # drop inf/nan if any
  return(X)
}

# #' Generate a YNM Series with Beta Noise
# #'
# #' Generates a YNM series where the noise follows a Beta distribution.
# #'
# #' @details
# #' The series is generated as:
# #' \deqn{X_k = \max(Y_k), \quad Y_k \sim \text{Beta}(\text{shape1},
# #' \text{shape2}), \quad \text{length}(Y_k) = \gamma^k}
# #' where:
# #' - \eqn{\gamma} controls the sample size growth.
# #' - \eqn{X_k} represents the maximum value of the Beta-distributed sample.
# #'
# #' @param T Integer. The length of the series.
# #' @param gamma Positive numeric. Growth parameter controlling sample size \eqn{\gamma \geq 1}.
# #' @param shape1 Positive numeric. The first shape parameter of the Beta distribution.
# #' @param shape2 Positive numeric. The second shape parameter of the Beta distribution.
# #' @return A numeric vector representing the YNM series.
# #' @export
# #' @examples
# #' YNM_series_Beta(100, gamma = 1.1, shape1 = 2, shape2 = 5)
# YNM_series_Beta <- function(T, gamma, shape1 = 1, shape2 = 1) {
# if (shape1 <= 0 | shape2 <= 0) stop("Enter positive values for shape parameters")

# X <- numeric(T)
# for (k in 1:T) {
# y <- rbeta(gamma^k, shape1, shape2)
# X[k] <- max(y)
# }
# return(X)
# }

# #' Generate a YNM Series with Gumbel Noise
# #'
# #' @details
# #' \deqn{X_i = \text{loc} - \text{scale} \cdot \log(-\log(U^{1/\gamma^i}))}
# #' where \eqn{U \sim \text{Uniform}(0,1)}.
# #'
# #' @inheritParams YNM_series_Beta
# #' @param loc Numeric. The location parameter of the Gumbel distribution.
# #' @param scale Positive numeric. The scale parameter of the Gumbel distribution.
# #' @return A numeric vector representing the YNM series.
# #' @export
# YNM_series_Gumbel <- function(T, gamma, loc = 0, scale = 1) {
# if (scale <= 0) stop("Enter a positive value for scale")

# y <- numeric(T)
# for (i in 1:T) {
# y[i] <- loc - scale * log(-log(runif(1)^(1/gamma^i)))
# }
# return(y)
# }

# #' Generate a YNM Series with Weibull Noise
# #'
# #' @details
# #' \deqn{X_i = (-\text{scale}^\text{shape} \cdot \log(1 - U^{1/\gamma^i}))^{1/\text{shape}}}
# #' where \eqn{U \sim \text{Uniform}(0,1)}.
# #'
# #' @inheritParams YNM_series_Beta
# #' @param scale Positive numeric. The scale parameter of the Weibull distribution.
# #' @param shape Positive numeric. The shape parameter of the Weibull distribution.
# #' @return A numeric vector representing the YNM series.
# #' @export
# YNM_series_Weibull <- function(T, gamma, scale = 1, shape = 1) {
# if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")

# X <- numeric(T)
# for (i in 1:T) {
# X[i] <- (-scale^shape * log(1 - runif(1)^(1/gamma^i)))^(1/shape)
# }
# return(X)
# }

# #' Generate a YNM Series with Frechet Noise
# #'
# #' @details
# #' \deqn{X_i = \text{scale} \cdot \left(-\frac{1}{\gamma^i} \log(U) \right)^{-1/\text{shape}}}
# #' where \eqn{U \sim \text{Uniform}(0,1)}.
# #'
# #' @inheritParams YNM_series_Beta
# #' @param scale Positive numeric. The scale parameter of the Frechet distribution.
# #' @param shape Positive numeric. The shape parameter of the Frechet distribution.
# #' @return A numeric vector representing the YNM series.
# #' @export
# YNM_series_Frechet <- function(T, gamma, shape, scale) {
# if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")

# y <- numeric(T)
# for (i in 1:T) {
# y[i] <- scale * (-(1/gamma^i) * log(runif(1)))^(-1/shape)
# }
# return(y)
# }

# #' Generate a YNM Series with Exponential Noise
# #'
# #' @details
# #' \deqn{X_i = -\frac{\log(1 - U^{1/\gamma^i})}{\text{rate}}}
# #' where \eqn{U \sim \text{Uniform}(0,1)}.
# #'
# #' @inheritParams YNM_series_Beta
# #' @param rate Positive numeric. The rate parameter of the Exponential distribution.
# #' @return A numeric vector representing the YNM series.
# #' @export
# YNM_series_Exp <- function(T, gamma, rate = 1) {
# if (rate <= 0) stop("Enter a positive value for rate")

# y <- numeric(T)
# for (i in 1:T) {
# y[i] <- -log(1 - runif(1)^(1/gamma^i)) / rate
# }
# return(y)
# }

# #' Generate a YNM Series with Pareto Noise
# #'
# #' @details
# #' \deqn{X_i = \text{scale} \cdot (1 - U^{1/\gamma^i})^{-1/\text{shape}}}
# #' where \eqn{U \sim \text{Uniform}(0,1)}.
# #'
# #' @inheritParams YNM_series_Beta
# #' @inheritParams YNM_series_Weibull
# #' @return A numeric vector representing the YNM series.
# #' @export
# YNM_series_Pareto <- function(T, gamma, scale = 0.5, shape = 4) {
# if (scale <= 0 | shape <= 0) stop("Enter positive values for scale and shape")

# X <- numeric(T)
# for (i in 1:T) {
# X[i] <- scale * (1 - runif(1)^(1/gamma^i))^(-1/shape)
# }
# return(X)
# }

# #' Generate a YNM Series with Normal Noise
# #'
# #' @details
# #' The series is generated as:
# #' \deqn{X_i = \max(Y_i), \quad Y_i \sim \mathcal{N}(\text{loc}, \text{sd}), \quad \text{length}(Y_i) = \gamma^i}
# #'
# #' @inheritParams YNM_series_Beta
# #' @param mean Numeric. The mean (location) of the normal distribution.
# #' @param sd Positive numeric. The standard deviation of the normal distribution.
# #' @return A numeric vector representing the YNM series.
# #' @export
# YNM_series_Norm <- function(T, gamma, mean = 0, sd = 1) {
# if (sd <= 0) stop("Enter a positive value for standard deviation")

# # Pre-allocate
# y <- numeric(T)

# for (i in 1:T) {
# m <- gamma^i
# u <- runif(1)                # uniform random variable
# y[i] <- qnorm(u^(1/m), mean = mean, sd = sd)
# }
# y= y[is.finite(y)]

# return(y)
# }



##  DTRW -----------------------
#' Generate a DTRW Process
#'
#' Simulates a discrete-time random Walk process of length T under different
#' underlying distributions.
#'
#' For DTRW without a drift, the distribution of the errors shall be continuous
#' and symmetric.
#'
#' @param T Integer, length of the series.
#' @param dist Character, distribution name. One of:
#'  "norm", "cauchy, "uniform"
#' @param ... Additional parameters specific to the chosen distribution:
#'   \describe{
#'     \item{norm}{`mean`, `sd`}
#'     \item{cauchy}{`loc`, `scale`}
#'     \item{uniform}{`min`, `max`}
#'   }
#'
#' @return A numeric vector of length T, the simulated DTRW process.
#' @examples
#' DTRW_series(10,  dist = "cauchy", loc = 0, scale = 1)
#' # [1] -0.6905644  2.1214308  2.7874249  4.1135190  3.3054300  2.6729198
#' # 2.5556620  1.6442579 15.6275793 15.4525462
#'
#' DTRW_series(100,  dist = "uniform", min = -1, max = 1)
#' DTRW_series(100,  dist = "norm", loc = 0, sd = 1)
#' @export
DTRW_series <- function(T, dist = c("norm", "cauchy", "uniform"), ...) {
  dist <- match.arg(dist)
  args <- list(...)
  X <- numeric(T)

  for (i in 1:T) {

    e <- switch(
      dist,

      norm = {
        loc <- args$mean %||% 0
        sd  <- args$sd  %||% 1
        if (sd <= 0) stop("Enter positive value for sd")
        rnorm(T, loc, sd)

      },

      cauchy = {
        loc <- args$loc %||% 0
        scale  <- args$scale  %||% 1
        if (scale <= 0) stop("Enter positive value for scale")
        rcauchy(T, loc, scale)
      },

      uniform = {
        min <- args$min %||% -1
        max <- args$max %||% 1
        runif(T, min, max)
      }
    )
  }

  X = cumsum(e)
  X <- X[is.finite(X)] # drop inf/nan if any
  return(X)
}

# #' Generate a Discrete-Time Random Walk (DTRW) Series with Cauchy Noise
# #'
# #' Generates a DTRW series where the increments follow a Cauchy distribution.
# #'
# #' @details
# #' The series is generated recursively as:
# #' \deqn{x_t = x_{t-1} + e_t, \quad e_t \sim \text{Cauchy}(\text{loc}, \text{scale})}
# #' where:
# #' - \eqn{T} is the length of the series.
# #' - \eqn{e_t} are i.i.d. random variables from a Cauchy distribution.
# #' - \eqn{\text{loc}} is the location parameter.
# #' - \eqn{\text{scale}} is the scale parameter.
# #'
# #' @param T Integer. The length of the series.
# #' @param loc Numeric. The location parameter of the Cauchy distribution (default = 0).
# #' @param scale Positive numeric. The scale parameter of the Cauchy distribution (default = 1).
# #' @return A numeric vector representing the DTRW series.
# #' @export
# #' @examples
# #' DTRW_series_Cauchy(100, loc = 0, scale = 1)
# DTRW_series_Cauchy <- function(T, loc = 0, scale = 1) {
  # if (scale <= 0) stop("Enter a positive value for scale")
  # e <- rcauchy(T, location = loc, scale = scale) ## Generate increments
  # x <- numeric(T)

  # for (i in 2:T) {
    # x[i] <- x[i - 1] + e[i]
  # }

  # return(x)
# }

# #' Generate a Discrete-Time Random Walk (DTRW) Series with Normal Noise
# #'
# #' Generates a DTRW series where the increments follow a Normal distribution.
# #'
# #' @details
# #' The series is generated recursively as:
# #' \deqn{x_t = x_{t-1} + e_t, \quad e_t \sim \mathcal{N}(\text{mean}, \text{sd})}
# #' where:
# #' - \eqn{T} is the length of the series.
# #' - \eqn{e_t} are i.i.d. random variables from a Normal distribution.
# #' - \eqn{\text{mean}} is the mean of the Normal distribution.
# #' - \eqn{\text{sd}} is the standard deviation of the Normal distribution.
# #'
# #' @param T Integer. The length of the series.
# #' @param mean Numeric. The mean (location) parameter of the normal distribution.
# #' @param sd Positive numeric. The standard deviation of the normal distribution.
# #' @return A numeric vector representing the DTRW series.
# #' #@export
# #' @examples
# #' DTRW_series_Norm(100, mean = 0, sd = 1)
# DTRW_series_Norm <- function(T, mean, sd) {
#  if (sd <= 0) stop("Enter a positive value for standard deviation")
#  # e <- rnorm(T, mean = mean, sd = sd) ## Generate increments
#  # x <- numeric(T)
#
#  x <- cumsum(rnorm(T, mean = mean, sd = sd))
#  #x <- x - x[1] # Normalize the series so it starts at 0
#  return(x)
# }

