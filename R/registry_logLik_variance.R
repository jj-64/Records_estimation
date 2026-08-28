var_logLik_records <- function(model, obs_type, dist, param_name) {
  # lookup
  f <- var_registry[[model]][[obs_type]][[dist]][[param_name]]

  if (is.null(f))
    stop("variance expression not registered for this (model, obs_type, dist, param_name).")

  #f(data, params)
  return(f)
}
library(MASS)
## Helper ---------------------------------
# empty container
var_registry <- new.env(parent = emptyenv())

## Function to register likelihood expressions later
register_var <- function(model, obs_type, dist, param_name, fun) {
  if (!exists(model, envir = var_registry))
    var_registry[[model]] <- list()

  if (!obs_type %in% names(var_registry[[model]]))
    var_registry[[model]][[obs_type]] <- list()

  var_registry[[model]][[obs_type]][[dist]][[param_name]] <- fun
}

## Classical, Xt ---------------

# register_var( model = "iid", obs_type = "all", dist = "frechet_inv_scale", param_name = "shape",
#               fun = function(data, params) {
#
#                 if(!is.numeric(data)) stop("data should be a numerical vector")
#                 n <- length(data)
#
#                 if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present.")
#
#                 shape <- as.numeric(params$shape)
#                 scale <- 1/as.numeric(params$scale)
#
#                 s1=-n/(shape^2)
#
#                 s2 = - sum( log(scale*data)^2/(scale*data)^shape )
#
#                 var_a = -1 / (s1+s2)
#
#                 return(var_a)
#
#
#               }
# )
#
# register_var( model = "iid", obs_type = "all", dist = "frechet_inv_scale", param_name = "scale",
#               fun = function(data, params) {
#
#                 if(!is.numeric(data)) stop("data should be a numerical vector")
#                 n <- length(data)
#
#                 if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present.")
#
#                 shape <- as.numeric(params$shape)
#                 scale <- 1/as.numeric(params$scale)
#
#                 s1 = (shape*n)/scale^2
#
#                 s2 = -((shape^2+shape)/scale^2) * sum((data*scale)^(-shape))
#
#                 var_A <- -1 / (s1 + s2)
#                 var_scale = var_A / scale^4
#
#                 return(var_scale)
#
#               }
# )
# # without hessian matrix
# register_var( model = "iid", obs_type = "all", dist = "frechet_inv_scale", param_name = "all",
#     fun <- function(data, params) {
#
#       # --- Checks ---
#       if (!is.numeric(data)) stop("data should be a numerical vector")
#
#       if (!all(c("shape", "scale") %in% names(params))) {
#         stop("parameters shape and scale should be present.")
#       }
#
#       # --- Inputs ---
#       n <- length(data)
#
#       shape <- as.numeric(params$shape)
#       scale_orig <- as.numeric(params$scale)  # original scale
#       scale <- 1 / scale_orig                # inverted scale (your parametrization)
#
#       # ________________________
#       # Variance of shape
#       # ________________________
#       s1_shape <- -n / (shape^2)
#       s2_shape <- -sum((log(scale * data)^2) / (scale * data)^shape)
#
#       var_shape <- -1 / (s1_shape + s2_shape)
#
#       # ________________________
#       # Variance of scale
#       # ________________________
#       s1_scale <- (shape * n) / scale^2
#       s2_scale <- -((shape^2 + shape) / scale^2) * sum((data * scale)^(-shape))
#
#       var_A <- -1 / (s1_scale + s2_scale)
#
#       # Transform back to original scale parameter
#       var_scale <- var_A / scale^4
#
#       # --- Safety check ---
#       if (any(is.nan(c(var_shape, var_scale))) || any(!is.finite(c(var_shape, var_scale)))) {
#         return(c(shape = -Inf, scale = -Inf))
#       }
#
#       # --- Return named vector ---
#       return(c(shape = var_shape, scale = var_scale))
#     }
# )

register_var( model = "iid", obs_type = "all", dist = "frechet_inv_scale", param_name = "all",
  fun <- function(data, params) {

    # --- Checks ---
    if (!all(c("shape", "scale") %in% names(params))) {
      stop("parameters shape and scale should be present")
    }

    # --- Extract data ---
    x <- data
    if(!is.numeric(x)) stop("data should be a numerical vector")
    T_val <- length(x)

    alpha <- as.numeric(params$shape)
    A <- as.numeric(params$scale)

    # --- Precomputations ---
    Ax <- A * x
    logAx <- log(Ax)

    # ________________________====
    # Hessian components
    # ________________________====

    # d²ℓ / dA²
    d2_A2 <- (alpha * T_val) / (A^2) -
      alpha * (alpha + 1) * sum(x^2 * Ax^(-alpha - 2))

    # d²ℓ / dα²
    d2_alpha2 <- -T_val / (alpha^2) -
      sum(Ax^(-alpha) * logAx^2)

    # cross derivative
    d2_Aalpha <- -T_val / A +
      sum(x * Ax^(-alpha - 1) * (1 - alpha * logAx))

    # --- Hessian matrix ---
    H <- matrix(c(
      d2_A2,       d2_Aalpha,
      d2_Aalpha,   d2_alpha2
    ), nrow = 2, byrow = TRUE)

    colnames(H) <- rownames(H) <- c("scale", "shape")

    # ____ Fisher Information (observed) ____
    Fisher <- -H
    #vcov_mat = solve(Fisher)
    vcov_mat <- tryCatch({
      solve(Fisher)
    }, error = function(e) {
      warning("Matrix singular. Falling back to Moore-Penrose pseudo-inverse.")
      # ginv handles near-singular matrices smoothly
      MASS::ginv(Fisher)
    })

    #var_location = -1/d2_mu2, var_scale = -1/d2_sigma2,
    return(list(Hessian = H, vcov = vcov_mat , var_scale = diag(vcov_mat)[1], var_shape = diag(vcov_mat)[2]))
  }
)


##Gumbel
register_var( model = "iid", obs_type = "all", dist = "gumbel", param_name = "all",
  fun <- function(data, params) {

    if (!all(c("location", "scale") %in% names(params))) {
      stop("parameters must include location and scale")
    }

    # --- Extract ---
    x <- data
    if(!is.numeric(x)) stop("data should be a numerical vector")
    mu <- params[["location"]]
    sigma <- params[["scale"]]

    if (sigma <= 0) stop("scale must be > 0")

    Tn <- length(x)
    y <- x - mu
    z <- y / sigma
    e <- exp(-z)

    # ________________________
    # Gradient
    # ________________________

    d_mu <- (Tn / sigma) - sum(e) / sigma

    d_sigma <- -Tn / sigma +
      sum(y) / sigma^2 -
      sum(y * e) / sigma^2

    gradient <- c(location = d_mu,
                  scale = d_sigma)

    # ________________________
    # Hessian
    # ________________________

    d2_mu2 <- -Tn / sigma^2 -
      sum(e) / sigma^2

    d2_sigma2 <- Tn / sigma^2 -
      2 * sum(y) / sigma^3 +
      2 * sum(y * e) / sigma^3 +
      sum((y^2) * e) / sigma^4

    d2_mu_sigma <- -Tn / sigma^2 +
      sum(e) / sigma^2 -
      sum(y * e) / sigma^3

    H <- matrix(c(
      d2_mu2, d2_mu_sigma,
      d2_mu_sigma, d2_sigma2
    ), nrow = 2, byrow = TRUE)

    rownames(H) <- c("location", "scale")
    colnames(H) <- c("location", "scale")

    Fisher = - H
    #vcov_mat = solve(Fisher)
    vcov_mat <- tryCatch({
      solve(Fisher)
    }, error = function(e) {
      warning("Matrix singular. Falling back to Moore-Penrose pseudo-inverse.")
      # ginv handles near-singular matrices smoothly
      MASS::ginv(Fisher)
    })

    return(list(
      gradient = gradient,
      hessian = H,
      vcov = vcov_mat,
      var_location = diag(vcov_mat)[1],
      var_scale = diag(vcov_mat)[2]
    ))
  }
)
## Classical, Rn ---------------

## Frechet
register_var( model = "iid", obs_type = "records", dist = "frechet", param_name = "scale",
              fun = function(data, params) {
                if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                Rn <- data$rec_values
                Ln <- data$rec_times
                T_val <- data$time[1]
                m = length(Rn)

                if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present.")

                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                s1=-(m-1)/(scale^2)

                s2 = -((shape^2+shape)/scale^2) * sum((Rn[-1]*scale)^(-shape))

                s3 = (shape+1)*(m-1)/scale^2

                # s4a=0
                # for (i in 1:(m - 1)) {
                #   s4a[i] = (Ln[i+1]-Ln[i]-1)*(scale*Rn[i])^(-shape)
                # }
                s4a <- (Ln[-1] - Ln[-m] - 1) * (scale * Rn[-m])^(-shape)
                s4 = -(shape^2+shape)*sum(s4a)/scale^2

                ## case where last record is not last observation
                s5=0
                if( (Ln[m] < T_val) == TRUE) {
                  s5 = -(shape*(shape+1))*(T_val-Ln[m])*(scale*Rn[m])^(-shape)/scale^2
                }

                #var_A <- -1 / (s1 + s2 + s3 + s4 + s5)
                info <- s1 + s2 + s3 + s4 + s5
                if (info >= 0) warning("Observed information is non-negative; variance may be invalid")

                var_A <- -1 / info
                var_scale = var_A / scale^4
                return(var_scale)
              }
)
register_var( model = "iid", obs_type = "records", dist = "frechet", param_name = "shape",
              fun = function(data, params) {
                if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                Rn <- data$rec_values
                Ln <- data$rec_times
                T_val <- data$time[1]

                if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present.")

                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                m = length(Rn)

                s1 = -(m - 1) / shape^2
                s2 = -sum((log(scale * Rn[-1]))^2 / (scale * Rn[-1])^shape)

                s3a = 0
                for (i in 1:(m - 1)) {
                  s3a[i] = (Ln[i + 1] - Ln[i] - 1) * (scale * Rn[i])^(-shape) * (log(scale * Rn[i]))^2
                }
                s3 = -sum(s3a)

                ## Case where last record is not the last observation
                s4 = 0
                if (Ln[m] < T_val) {
                  s4 = -(T_val - Ln[m]) * (scale * Rn[m])^(-shape) * (log(scale * Rn[m]))^2
                }

                var_a = -1 / (s1 + s2 + s3 + s4)

                return(var_a)
              }
)

register_var( model = "iid", obs_type = "records", dist = "frechet", param_name = "all",
              fun <- function(data, params) {

                # --- Checks ---
                if (!all(c("rec_values", "rec_times", "time") %in% names(data))) {
                  stop("data must contain rec_values, rec_times, and time")
                }

                if (!all(c("shape", "scale") %in% names(params))) {
                  stop("parameters shape and scale should be present")
                }

                # --- Extract data ---
                Rn <- data$rec_values
                Ln <- data$rec_times
                T_val <- data$time[1]
                m <- length(Rn)

                shape <- as.numeric(params$shape)
                scale_orig <- as.numeric(params$scale)
                scale <- 1 / scale_orig   # same parametrization as your code

                # ________________________
                # Variance of scale
                # ________________________
                s1_sc <- -(m - 1) / (scale^2)
                s2_sc <- -((shape^2 + shape) / scale^2) * sum((Rn[-1] * scale)^(-shape))
                s3_sc <- (shape + 1) * (m - 1) / scale^2

                s4a <- (Ln[-1] - Ln[-m] - 1) * (scale * Rn[-m])^(-shape)
                s4_sc <- -((shape^2 + shape) / scale^2) * sum(s4a)

                s5_sc <- 0
                if (Ln[m] < T_val) {
                  s5_sc <- -(shape * (shape + 1)) * (T_val - Ln[m]) *
                    (scale * Rn[m])^(-shape) / scale^2
                }

                info_scale <- s1_sc + s2_sc + s3_sc + s4_sc + s5_sc

                if (info_scale >= 0) warning("Scale information is non-negative")

                var_scale_internal <- -1 / info_scale    # variance in inverted scale
                # If you want back to original scale, uncomment:
                # var_scale <- var_scale_internal / scale^4
                var_scale <- var_scale_internal


                # ________________________
                # Variance of shape
                # ________________________
                s1_sh <- -(m - 1) / shape^2

                s2_sh <- -sum((log(scale * Rn[-1]))^2 / (scale * Rn[-1])^shape)

                s3a <- (Ln[-1] - Ln[-m] - 1) *
                  (scale * Rn[-m])^(-shape) *
                  (log(scale * Rn[-m]))^2
                s3_sh <- -sum(s3a)

                s4_sh <- 0
                if (Ln[m] < T_val) {
                  s4_sh <- -(T_val - Ln[m]) *
                    (scale * Rn[m])^(-shape) *
                    (log(scale * Rn[m]))^2
                }

                var_shape <- -1 / (s1_sh + s2_sh + s3_sh + s4_sh)

                # --- Safety check ---
                if (any(is.nan(c(var_shape, var_scale))) ||
                    any(!is.finite(c(var_shape, var_scale)))) {
                  return(c(shape = -Inf, scale = -Inf))
                }

                # --- Output ---
                return(c(shape = var_shape, scale = var_scale))
              }
)


# register_var( model = "iid", obs_type = "records", dist = "frechet_inv_scale", param_name = "scale",
#               fun = function(data, params) {
#                 if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")
#
#                 Rn <- data$rec_values
#                 Ln <- data$rec_times
#                 T_val <- data$time[1]
#                 m = length(Rn)
#
#                 if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present.")
#
#                 shape <- as.numeric(params$shape)
#                 scale <- 1/as.numeric(params$scale)
#
#                 s1=-(m-1)/(scale^2)
#
#                 s2 = -((shape^2+shape)/scale^2) * sum((Rn[-1]*scale)^(-shape))
#
#                 s3 = (shape+1)*(m-1)/scale^2
#
#                 # s4a=0
#                 # for (i in 1:(m - 1)) {
#                 #   s4a[i] = (Ln[i+1]-Ln[i]-1)*(scale*Rn[i])^(-shape)
#                 # }
#                 s4a <- (Ln[-1] - Ln[-m] - 1) * (scale * Rn[-m])^(-shape)
#                 s4 = -(shape^2+shape)*sum(s4a)/scale^2
#
#                 ## case where last record is not last observation
#                 s5=0
#                 if( (Ln[m] < T_val) == TRUE) {
#                   s5 = -(shape*(shape+1))*(T_val-Ln[m])*(scale*Rn[m])^(-shape)/scale^2
#                 }
#
#                 #var_A <- -1 / (s1 + s2 + s3 + s4 + s5)
#                 info <- s1 + s2 + s3 + s4 + s5
#                 if (info >= 0) warning("Observed information is non-negative; variance may be invalid")
#
#                 var_A <- -1 / info
#                 #var_scale = var_A / scale^4
#                 return(var_A)
#               }
# )

## without Hessian matrix
# register_var( model = "iid", obs_type = "records", dist = "frechet_inv_scale", param_name = "all",
#   fun <- function(data, params) {
#
#     # --- Checks ---
#     if (!all(c("rec_values", "rec_times", "time") %in% names(data))) {
#       stop("data must contain rec_values, rec_times, and time")
#     }
#
#     if (!all(c("shape", "scale") %in% names(params))) {
#       stop("parameters shape and scale should be present")
#     }
#
#     # --- Extract data ---
#     Rn <- data$rec_values
#     Ln <- data$rec_times
#     T_val <- data$time[1]
#     m <- length(Rn)
#
#     shape <- as.numeric(params$shape)
#     scale_orig <- as.numeric(params$scale)
#     scale <- 1 / scale_orig   # same parametrization as your code
#
#
#     # ________________________
#     # Variance of scale
#     # ________________________
#
#     s1_sc <- -(m - 1) / (scale^2)
#     s2_sc <- -((shape^2 + shape) / scale^2) * sum((Rn[-1] * scale)^(-shape))
#     s3_sc <- (shape + 1) * (m - 1) / scale^2
#
#     s4a <- (Ln[-1] - Ln[-m] - 1) * (scale * Rn[-m])^(-shape)
#     s4_sc <- -((shape^2 + shape) / scale^2) * sum(s4a)
#
#     s5_sc <- 0
#     if (Ln[m] < T_val) {
#       s5_sc <- -(shape * (shape + 1)) * (T_val - Ln[m]) *
#         (scale * Rn[m])^(-shape) / scale^2
#     }
#
#     info_scale <- s1_sc + s2_sc + s3_sc + s4_sc + s5_sc
#
#     if (info_scale >= 0) warning("Scale information is non-negative")
#
#     var_scale_internal <- -1 / info_scale    # variance in inverted scale
#     # If you want back to original scale, uncomment:
#     # var_scale <- var_scale_internal / scale^4
#     var_scale <- var_scale_internal
#
#
#     # ________________________
#     # Variance of shape
#     # ________________________
#     s1_sh <- -(m - 1) / shape^2
#
#     s2_sh <- -sum((log(scale * Rn[-1]))^2 / (scale * Rn[-1])^shape)
#
#     s3a <- (Ln[-1] - Ln[-m] - 1) *
#       (scale * Rn[-m])^(-shape) *
#       (log(scale * Rn[-m]))^2
#     s3_sh <- -sum(s3a)
#
#     s4_sh <- 0
#     if (Ln[m] < T_val) {
#       s4_sh <- -(T_val - Ln[m]) *
#         (scale * Rn[m])^(-shape) *
#         (log(scale * Rn[m]))^2
#     }
#
#     var_shape <- -1 / (s1_sh + s2_sh + s3_sh + s4_sh)
#
#     # --- Safety check ---
#     if (any(is.nan(c(var_shape, var_scale))) ||
#         any(!is.finite(c(var_shape, var_scale)))) {
#       return(c(shape = -Inf, scale = -Inf))
#     }
#
#     # --- Output ---
#     return(c(shape = var_shape, scale = var_scale))
#   }
# )

register_var( model = "iid", obs_type = "records", dist = "frechet_inv_scale", param_name = "all",
              fun <- function(data, params) {

                # --- Checks ---
                if (!all(c("rec_values", "rec_times", "time") %in% names(data))) {
                  stop("data must contain rec_values, rec_times, and time")
                }

                if (!all(c("shape", "scale") %in% names(params))) {
                  stop("parameters shape and scale should be present")
                }

                # --- Extract data ---
                Rn <- data$rec_values      # record values
                Ln <- data$rec_times       # record times
                T_val <- data$time[1]
                m <- length(Rn)

                alpha <- as.numeric(params$shape)
                A <- as.numeric(params$scale)

                # --- Build sequences ---
                # x_{l_{n+1}} and x_{l_n}
                x_ln1 <- Rn[-1]
                x_ln  <- Rn[-length(Rn)]

                # Delta_{n+1}
                delta_n1 <- Ln[-1] - Ln[-length(Ln)]

                # last time
                l_NT <- Ln[length(Ln)]

                # weights
                w_n <- delta_n1 - 1
                w_T <- (T_val - l_NT)

                NTm1 <- length(x_ln1)

                # --- Precomputations ---
                Ax1 <- A * x_ln1
                Axn <- A * x_ln
                AxT <- A * Rn[length(Rn)]

                log1 <- log(Ax1)
                logn <- log(Axn)
                logT <- log(AxT)

                # ________________________====
                # Hessian components
                # ________________________====

                # ____ d²ℓ / dA² ____
                d2_A2 <- (alpha * NTm1) / (A^2) -
                  alpha * (alpha + 1) * sum(x_ln1^2 * Ax1^(-alpha - 2)) -
                  alpha * (alpha + 1) * sum(w_n * x_ln^2 * Axn^(-alpha - 2)) -
                  alpha * (alpha + 1) * w_T * Rn[length(Rn)]^2 * AxT^(-alpha - 2)

                # ____ d²ℓ / dα² ____
                d2_alpha2 <- -NTm1 / (alpha^2) -
                  sum(Ax1^(-alpha) * log1^2) -
                  sum(w_n * Axn^(-alpha) * logn^2) -
                  w_T * AxT^(-alpha) * logT^2

                # ____ cross derivative ____
                term1 <- x_ln1 * Ax1^(-alpha - 1) * (1 - alpha * log1)
                termn <- x_ln * Axn^(-alpha - 1) * (1 - alpha * logn)
                termT <- Rn[length(Rn)] * AxT^(-alpha - 1) * (1 - alpha * logT)

                d2_Aalpha <- -NTm1 / A +
                  sum(term1) +
                  sum(w_n * termn) +
                  w_T * termT

                # ________________________====
                # Hessian matrix
                # ________________________====

                H <- matrix(c(
                  d2_A2,       d2_Aalpha,
                  d2_Aalpha,   d2_alpha2
                ), nrow = 2, byrow = TRUE)

                colnames(H) <- rownames(H) <- c("scale", "shape")

                # ____ Fisher Information (observed) ____
                Fisher <- -H
                #vcov_mat = solve(Fisher)
                vcov_mat <- tryCatch({
                  solve(Fisher)
                }, error = function(e) {
                  warning("Matrix singular. Falling back to Moore-Penrose pseudo-inverse.")
                  # ginv handles near-singular matrices smoothly
                  MASS::ginv(Fisher)
                })

                #var_location = -1/d2_mu2, var_scale = -1/d2_sigma2,
                return(list(Hessian = H, vcov = vcov_mat , var_scale = diag(vcov_mat)[1], var_shape = diag(vcov_mat)[2]))
              }
)


## Gumbel
register_var(model = "iid", obs_type = "records", dist = "gumbel", param_name = "all",
  fun <- function(data, params) {

    # --- Checks ---
    if (!all(c("rec_values", "rec_times", "time") %in% names(data))) {
      stop("data must contain rec_values, rec_times, and time")
    }

    if (!all(c("location", "scale") %in% names(params))) {
      stop("parameters location and scale should be present")
    }

    # --- Extract data ---
    Rn <- data$rec_values     # x_{l_n}
    Ln <- data$rec_times      # l_n
    T_val <- data$time[1]
    m <- length(Rn)

    mu <- params[["location"]]
    sigma <- params[["scale"]]

    # --- Derived vectors ---
    x_next <- Rn[-1]          # x_{l_{n+1}}
    x_curr <- Rn[-m]          # x_{l_n}
    delta  <- diff(Ln)        # Delta_{n+1}

    last_weight <- T_val - Ln[m]

    # --- Standardized variables ---
    z_next <- (x_next - mu) / sigma
    z_curr <- (x_curr - mu) / sigma
    z_last <- (Rn[m] - mu) / sigma

    e_next <- exp(-z_next)
    e_curr <- exp(-z_curr)
    e_last <- exp(-z_last)

    # ________________________
    # Second derivatives
    # ________________________

    # ---- d2 / d mu^2 ____
    d2_mu2 <- #- (m - 1) / sigma^2
      - sum(e_next) / sigma^2 -
      sum((delta - 1) * e_curr) / sigma^2 -
      last_weight * e_last / sigma^2

    # ____ d2 / d sigma^2 ____
    d2_sigma2 <- (m - 1) / sigma^2 -
      2 * sum(x_next - mu) / sigma^3 +
      2 * sum((x_next - mu) * e_next) / sigma^3 -
      sum((x_next - mu)^2 * e_next) / sigma^4 +
      2 * sum((delta - 1) * (x_curr - mu) * e_curr) / sigma^3 -
      sum((delta - 1) * (x_curr - mu)^2 * e_curr) / sigma^4 +
      2 * last_weight * (Rn[m] - mu) * e_last / sigma^3 -
      last_weight * (Rn[m] - mu)^2 * e_last / sigma^4

    # ____ cross derivative d2 / d mu d sigma ____
    d2_mu_sigma <- - (m - 1) / sigma^2 +
      sum(e_next) / sigma^2 -
      sum((x_next - mu) * e_next) / sigma^3 +
      sum((delta - 1) * e_curr) / sigma^2 -
      sum((delta - 1) * (x_curr - mu) * e_curr) / sigma^3 +
      last_weight * e_last / sigma^2 -
      last_weight * (Rn[m] - mu) * e_last / sigma^3

    # --- Hessian matrix ---
    H <- matrix(c(
      d2_mu2,        d2_mu_sigma,
      d2_mu_sigma,   d2_sigma2
    ), nrow = 2, byrow = TRUE)

    rownames(H) <- c("location", "scale")
    colnames(H) <- c("location", "scale")


    # ____ Fisher Information (observed) ____
    Fisher <- -H
    #vcov_mat = solve(Fisher)
    vcov_mat <- tryCatch({
      solve(Fisher)
    }, error = function(e) {
      warning("Matrix singular. Falling back to Moore-Penrose pseudo-inverse.")
      # ginv handles near-singular matrices smoothly
      MASS::ginv(Fisher)
    })
    #var_location = -1/d2_mu2, var_scale = -1/d2_sigma2,
    return(list(Hessian = H, vcov = vcov_mat , var_location = diag(vcov_mat)[1], var_scale = diag(vcov_mat)[2]))
  }
)
## DTRW, Xt -------------
  ##norm
register_var( model = "dtrw", obs_type = "all", dist = "norm", param_name = "sd",
              fun = function(data, params) {
                if(!is.numeric(data)) stop("data should be a numerical vector")

                if( all(c("sd") %in% names(params)) == FALSE ) stop("parameter sd should be present.")

                n = length(data)-1
                sd <- as.numeric(params$sd) #1/

                #var_est =  2 * sd^4 / n  # variance of sigma2
                var_est = sd^2/(2 * n)
                return(var_est)
              }
)

register_var( model = "dtrw", obs_type = "all", dist = "unif", param_name = "max",
              fun = function(data, params, biased = TRUE) {
                if(!is.numeric(data)) stop("data should be a numerical vector")

                if( all(c("max") %in% names(params)) == FALSE ) stop("parameter max should be present.")

                n = length(data)-1
                # max_est = max(abs(diff(data)))
                # max_est_corr = (n+1)/n * max_est
                max_est = params[["max"]]
                if(biased){
                  var_est = n / ((n + 1)^2 * (n + 2)) * max_est^2
                } else {
                  var_est =  max_est^2 / (n * (n + 2))
                }
                return(var_est )
              }
)
## DTRW, Rn -------------
  ##norm
register_var( model = "dtrw", obs_type = "records", dist = "norm", param_name = "sd",
              fun = function(data, params) {
                if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                Rn <- data$rec_values
                Ln <- data$rec_times
                T <- data$time

                if( all(c("sd") %in% names(params)) == FALSE ) stop("parameter sd should be present.")

                mean = 0
                sd <- as.numeric(params$sd) #1/

                m = length(Rn)
                #return(sd^2/(2*sqrt(m)))
                #var_est = 2*sd^4/(m-1) # variance of sigma2
                var_est = sd^2/(2*(m-1))
                return(var_est )
              }
)

## uniform
register_var( model = "dtrw", obs_type = "records", dist = "unif", param_name = "max",
              fun = function(data, params, biased = TRUE) {
                if( all(c("rec_values") %in% names(data)) == FALSE ) stop("a list of rec_values should be present.")

                Rn <- data$rec_values
                #Ln <- data$rec_times
                #T <- data$time

                if( all(c("max") %in% names(params)) == FALSE ) stop("parameter max should be present.")

                m <- length(Rn)
                #b_hat <- max(abs(diff(Rn)))
                #b_hat_corr <- m/(m-1) * b_hat
                b_hat = params[["max"]]
                if(biased){
                var_est = (m - 1) / (m^2 * (m + 1)) * b_hat^2
                } else {
                var_est = b_hat^2/ ( (m-1) * (m+1))
                }

                return(var_est )
              }
)

## LDM, Xt --------------

  ##frechet
register_var( model = "LDM", obs_type = "all", dist = "frechet", param_name = "theta",
              fun = function(data, params) {
                if(!is.numeric(data)) stop("data should be a numerical vector")
                n = length(data)

                if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")

                theta <- as.numeric(params$theta)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                y = data-theta*(1:n) ## frechet iid
                s1 = scale^(-shape)*shape*(-shape-1) * sum((1:n)^2 * y^(-shape-2))
                s2 = (shape+1)*sum(((1:n)/y)^2)
                fisher= s1+s2
                return(-1/fisher)

              }
)

register_var( model = "LDM", obs_type = "all", dist = "frechet", param_name = "scale",
              fun = function(data, params) {
                if(!is.numeric(data)) stop("data should be a numerical vector")
                n = length(data)

                if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")

                theta <- as.numeric(params$theta)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                y = data-theta*(1:n) ## frechet iid
                s1 = -(shape+1)*shape*scale^(-shape-2)*sum(y^-shape)
                s2=shape*n/scale^2
                fisher= s1+s2
                return(-1/fisher)

              }
)

register_var( model = "LDM", obs_type = "all", dist = "frechet", param_name = "shape",
              fun = function(data, params) {
                if(!is.numeric(data)) stop("data should be a numerical vector")
                n = length(data)

                if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")

                theta <- as.numeric(params$theta)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                y = data-theta*(1:n) ## frechet iid
                p1 = y^(-shape)*(log(scale*y))^2
                s1 = -scale^(-shape)*sum(p1)
                s2 = -n/(shape^2)
                fisher =  s1+s2
                return(-1/fisher)

              }
)

register_var(model = "LDM",obs_type = "all",dist = "frechet", param_name = "all",   # 👈 unified
  fun = function(data, params) {

    # ________________________
    # 1. Checks
    # ________________________
    if (!is.numeric(data)) stop("data should be numeric")

    required <- c("theta", "shape", "scale")
    if (!all(required %in% names(params))) {
      stop("parameters theta, shape, scale must be present")
    }

    # ________________________
    # 2. Extract inputs
    # ________________________
    n <- length(data)

    theta <- as.numeric(params$theta)
    shape <- as.numeric(params$shape)
    scale <- 1 / as.numeric(params$scale)

    t_seq <- 1:n
    y <- data - theta * t_seq

    # safety check
    if (any(y <= 0, na.rm = TRUE)) return(Inf)

    # ________________________
    # 3. Compute Fisher info
    # ________________________

    #theta
        s1_t <- scale^(-shape) * shape * (-shape - 1) *
          sum(t_seq^2 * y^(-shape - 2))

        s2_t <- (shape + 1) * sum((t_seq / y)^2)

        var_theta = -1 /( s1_t + s2_t)

     # scale
        s1_A <- -(shape + 1) * shape * scale^(-shape - 2) *
          sum(y^(-shape))

        s2_A <- shape * n / scale^2

        var_scale = -1 / (s1_A + s2_A)


      # "shape"
        p1_a <- y^(-shape) * (log(scale * y))^2

        s1_a <- -scale^(-shape) * sum(p1_a)
        s2_a <- -n / (shape^2)

        var_shape = -1 / (s1_a + s2_a)

    # ________________________
    # 4. Return variance
    # ________________________

    #if (!is.finite(fisher) || fisher == 0) return(Inf)

    ## --- 8. Return ---
    return(c(
      var_theta = var_theta,
      var_shape = var_shape,
      var_scale = var_scale
    ))
  }
)

register_var(model = "LDM",obs_type = "all",dist = "frechet_inv_scale", param_name = "all",   # 👈 unified
  fun <- function(data, params) {

    # --- Checks ---
    if (!all(c("theta", "shape", "scale") %in% names(params))) {
      stop("params must contain theta, shape, scale")
    }

    # --- Extract ---
    x <- data
    if(!is.numeric(x)) stop("data should be a numerical vector")
    T_val <- length(x)

    theta <- params$theta
    alpha <- params$shape
    A <- params$scale

    # --- Precompute ---
    t <- seq_len(T_val)
    y <- x - theta * t
    z <- A * y
    logz <- log(z)

    #________________________
    # Second derivatives
    #________________________

    # theta-theta
    d2_tt <- +(1+alpha)*sum(t^2 / y^2) -
      alpha*(alpha+1)*A^2 * sum(t^2 * z^(-alpha-2))

    # A-A
    d2_AA <- alpha*T_val/A^2 -
      alpha*(alpha+1)*sum(y^2 * z^(-alpha-2))

    # alpha-alpha
    d2_aa <- -T_val/alpha^2 -
      sum(z^(-alpha) * logz^2)

    # theta-A
    d2_tA <- -alpha*sum(t * z^(-alpha-1)) +
      alpha*(alpha+1)*A*sum(t * y * z^(-alpha-2))

    # theta-alpha
    d2_ta <- sum(t / y) -
      A*sum(t * z^(-alpha-1)) +
      alpha*A*sum(t * z^(-alpha-1) * logz)

    # A-alpha
    d2_Aa <- -T_val/A +
      sum(y * z^(-alpha-1) * (1 - alpha*logz))

    # --- Hessian matrix ---
    H <- matrix(c(
      d2_tt, d2_tA, d2_ta,
      d2_tA, d2_AA, d2_Aa,
      d2_ta, d2_Aa, d2_aa
    ), nrow = 3, byrow = TRUE)

    colnames(H) <- rownames(H) <- c("theta", "scale", "shape")

    # ____ Fisher Information (observed) ____
    Fisher <- -H
    #vcov_mat = solve(Fisher)
    vcov_mat <- tryCatch({
      solve(Fisher)
    }, error = function(e) {
      warning("Matrix singular. Falling back to Moore-Penrose pseudo-inverse.")
      # ginv handles near-singular matrices smoothly
      MASS::ginv(Fisher)
    })

    #var_location = -1/d2_mu2, var_scale = -1/d2_sigma2,
    return(list(Hessian = H, vcov = vcov_mat , var_theta = diag(vcov_mat)[1],var_scale = diag(vcov_mat)[2], var_shape = diag(vcov_mat)[3]))
  }
)

## Gumbel
register_var(model = "LDM",obs_type = "all",dist = "gumbel", param_name = "all",   # 👈 unified
  fun <- function(data, params) {

    if (!all(c("theta", "location", "scale") %in% names(params))) {
      stop("params must contain theta, location, scale")
    }

    # --- Extract ---
    x <- data
    if(!is.numeric(x)) stop("data should be a numerical vector")
    T_val <- length(x)

    theta <- params[["theta"]]
    A <- params[["location"]]
    sigma <- params[["scale"]]

    t <- seq_len(T_val)
    y <- x - theta * t
    u <- (y - A) / sigma
    e <- exp(-u)

    # _____________________________
    # Hessian components
    # _____________________________

    # theta-theta
    d2_tt <- -sum(t^2 * e) / sigma^2

    # A-A
    d2_AA <- -sum(e) / sigma^2

    # sigma-sigma
    d2_aa <- T_val/sigma^2 -
      2*sum(y - A)/sigma^3 +
      2*sum((y - A)*e)/sigma^3 -
      sum((y - A)^2 * e)/sigma^4

    # theta-A
    d2_tA <- -sum(t * e) / sigma^2

    # theta-sigma
    d2_ta <- -sum(t)/sigma^2 +
      sum(t * e)/sigma^2 -
      sum(t * (y - A) * e)/sigma^3

    # A-sigma
    d2_Aa <- -T_val/sigma^2 +
      sum(e)/sigma^2 -
      sum((y - A)*e)/sigma^3

    # --- Hessian ---
    H <- matrix(c(
      d2_tt, d2_tA, d2_ta,
      d2_tA, d2_AA, d2_Aa,
      d2_ta, d2_Aa, d2_aa
    ), nrow = 3, byrow = TRUE)

    colnames(H) <- rownames(H) <- c("theta", "location", "scale")

    # ____ Fisher Information (observed) ____
    Fisher <- - H
    #vcov_mat = solve(Fisher)
    vcov_mat <- tryCatch({
      solve(Fisher)
    }, error = function(e) {
      warning("Matrix singular. Falling back to Moore-Penrose pseudo-inverse.")
      # ginv handles near-singular matrices smoothly
      MASS::ginv(Fisher)
    })

    #var_location = -1/d2_mu2, var_scale = -1/d2_sigma2,
    return(list(Hessian = H, vcov = vcov_mat , var_theta = diag(vcov_mat)[1] ,  var_location = diag(vcov_mat)[2],var_scale = diag(vcov_mat)[3]) )
  }
)
## LDM, Rn --------------

  ##frechet
register_var( model = "LDM", obs_type = "records", dist = "frechet", param_name = "theta",
                   fun = function(data, params) {
                     if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                     Rn <- data$rec_values
                     Ln <- data$rec_times
                     n <- data$time[1]

                     if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")

                     theta <- as.numeric(params$theta)
                     shape <- as.numeric(params$shape)
                     scale <- 1/as.numeric(params$scale)
                     m = length(Rn)

                     x= Rn - theta * Ln

                     ## 1) theta Variance
                     s1 =  (shape + 1) * sum(Ln^2 / x^2)
                     s2 = sum( Ln^2 * x^(-shape - 2))
                     s3a=0
                     for (i in 1:(m - 1)) {
                       # Inner summation from (l[i] + 1) to (l[i+1] - 1)
                       for (t in (Ln[i] + 1):(Ln[i + 1] - 1)) {
                         if( ((Ln[i]+1) <=(Ln[i+1]-1) ) == TRUE ){
                           s3a <- s3a + t^2 * (Rn[i] - theta * t)^(-shape - 2) }}
                     }

                        ## case where last record is not last observation
                     s3b=0
                     if( (Ln[m] < n) == TRUE) {
                       for (i in (Ln[m]+1):n) {
                         s3b[i] <- (i^2) * ((Rn[m] - theta * i)^(-shape - 2))
                       } }
                     s3=s3a+sum(na.omit(s3b))
                     s4 = -shape * (shape + 1) * scale^(-shape) * (s3+s2)

                     var_theta = -1 / (s1 + s4)

                   return(var_theta)
                   }
)

register_var( model = "LDM", obs_type = "records", dist = "frechet", param_name = "shape",
              fun = function(data, params) {
                if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                Rn <- data$rec_values
                Ln <- data$rec_times
                n <- data$time[1]

                if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")

                theta <- as.numeric(params$theta)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)
                m = length(Rn)

                x= Rn - theta * Ln

                ## variance of shape
                s1_h = - m / (shape^2)

                s2_h= sum((log(scale*x))^2 * x^(-shape))

                s3a_h=0
                for (i in 1:(m - 1)) {
                  # Inner summation from (l[i] + 1) to (l[i+1] - 1)
                  for (t in (Ln[i] + 1):(Ln[i + 1] - 1)) {
                    if( ((Ln[i]+1) <=(Ln[i+1]-1) ) == TRUE ){
                      term <- scale * (Rn[i] - theta * t)
                      s3a_h <- s3a_h + (term^(-shape)) * (log(term)^2)  }}
                }

                ## case where last record is not last observation
                s3b_h=0
                if( (Ln[m] < n) == TRUE) {
                  for (i in (Ln[m]+1):n) {
                    term <- scale * (Rn[m] - theta * t)
                    s3b_h[i] <-(term^(-shape)) * (log(term)^2)
                  } }

                s3_h=s3a_h+sum(na.omit(s3b_h))

                s4_h <- - (scale^(-shape)) * (s2_h + s3_h)

                var_shape <- -1 / (s1_h+s4_h)

                return(var_shape )
              }
)

register_var( model = "LDM", obs_type = "records", dist = "frechet", param_name = "scale",
              fun = function(data, params) {
                if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                Rn <- data$rec_values
                Ln <- data$rec_times
                n <- data$time[1]

                if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")
                theta <- as.numeric(params$theta)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)
                m = length(Rn)

                x= Rn - theta * Ln


                ## 1/Scale variance

                s1_s = m * shape / (scale^2)
                s2_s = sum(x^(-shape))

                s3a_s=0
                for (i in 1:(m - 1)) {
                  # Inner summation from (l[i] + 1) to (l[i+1] - 1)
                  for (t in (Ln[i] + 1):(Ln[i + 1] - 1)) {
                    if( ((Ln[i]+1) <=(Ln[i+1]-1) ) == TRUE ){
                      s3a_s <- s3a_s + (Rn[i] - theta * t)^(-shape ) }}
                }

                ## case where last record is not last observation
                s3b_s=0
                if( (Ln[m] < n) == TRUE) {
                  for (i in (Ln[m]+1):n) {
                    s3b_s[i] <- ((Rn[m] - theta * i)^(-shape))
                  } }

                s3_s=s3a_s+sum(na.omit(s3b_s))

                s4_s <- -shape * (shape + 1) * (scale^(-shape - 2)) * (s2_s + s3_s)

                var_scale <- -1 / (s1_s + s4_s)

                return(var_scale)
              }
)

register_var( model = "LDM", obs_type = "records", dist = "frechet", param_name = "all",
              fun =  function(data, params) {

                ## --- 1. Input checks ---
                required_data <- c("rec_values", "rec_times", "time")
                if (!all(required_data %in% names(data))) {
                  stop("data must contain: rec_values, rec_times, time")
                }

                required_params <- c("theta", "shape", "scale")
                if (!all(required_params %in% names(params))) {
                  stop("params must contain: theta, shape, scale")
                }

                ## --- 2. Extract data ---
                R <- data$rec_values
                L <- data$rec_times
                T <- data$time[1]

                theta <- as.numeric(params["theta"])
                a     <- as.numeric(params["shape"])
                s     <- as.numeric(params["scale"])   # original scale
                A     <- 1 / s                         # inverted scale

                m <- length(R)

                if (m == 0 || !is.finite(theta) || !is.finite(a) || !is.finite(s) ||
                    a <= 0 || s <= 0) {
                  return(c(var_theta = -Inf, var_shape = -Inf, var_scale = -Inf))
                }

                ## --- 3. Transform ---
                x <- R - theta * L
                if (any(x <= 0)) {
                  return(c(var_theta = -Inf, var_shape = -Inf, var_scale = -Inf))
                }

                ## --- 4. Common sums ---
                s3_theta <- 0
                s3_scale <- 0
                s3_shape <- 0

                if (m > 1) {
                  for (j in 1:(m - 1)) {

                    start <- L[j] + 1
                    end   <- L[j + 1] - 1

                    if (start <= end) {
                      t_seq <- start:end
                      val   <- R[j] - theta * t_seq

                      if (any(val <= 0)) {
                        return(c(var_theta = -Inf, var_shape = -Inf, var_scale = -Inf))
                      }

                      ## theta part
                      s3_theta <- s3_theta + sum(t_seq^2 * val^(-a - 2))

                      ## scale part
                      s3_scale <- s3_scale + sum(val^(-a))

                      ## shape part
                      term <- val / s   # since A*x = x/s
                      s3_shape <- s3_shape + sum((term^(-a)) * (log(term)^2))
                    }
                  }
                }

                ## --- tail contribution ---
                if (L[m] < T) {

                  t_seq <- (L[m] + 1):T
                  val   <- R[m] - theta * t_seq

                  if (any(val <= 0)) {
                    return(c(var_theta = -Inf, var_shape = -Inf, var_scale = -Inf))
                  }

                  s3_theta <- s3_theta + sum(t_seq^2 * val^(-a - 2))
                  s3_scale <- s3_scale + sum(val^(-a))

                  term <- val / s
                  s3_shape <- s3_shape + sum((term^(-a)) * (log(term)^2))
                }

                ## --- 5. Variance of theta ---
                s1_theta <- (a + 1) * sum(L^2 / x^2)
                s2_theta <- sum(L^2 * x^(-a - 2))

                ## A^(-a) = s^a
                s4_theta <- -a * (a + 1) * s^(a) * (s2_theta + s3_theta)

                var_theta <- -1 / (s1_theta + s4_theta)

                ## --- 6. Variance of scale (original s) ---
                s1_scale <- m * a / (s^2)
                s2_scale <- sum(x^(-a))

                ## A^(-a-2) = s^(a+2)
                s4_scale <- -a * (a + 1) * s^(a + 2) * (s2_scale + s3_scale)

                var_scale <- -1 / (s1_scale + s4_scale)

                ## --- 7. Variance of shape ---
                s1_shape <- -m / (a^2)
                s2_shape <- sum((log(x / s))^2 * x^(-a))

                ## A^(-a) = s^a
                s4_shape <- -s^(a) * (s2_shape + s3_shape)

                var_shape <- -1 / (s1_shape + s4_shape)

                ## --- 8. Return ---
                return(c(
                  var_theta = var_theta,
                  var_shape = var_shape,
                  var_scale = var_scale
                ))
              }
)

# register_var( model = "LDM", obs_type = "records", dist = "frechet_inv_scale_nohessian", param_name = "all",
#             fun = function(data, params) {
#       ## --- 1. Input checks ---
#       required_data <- c("rec_values", "rec_times", "time")
#       if (!all(required_data %in% names(data))) {
#         stop("data must contain: rec_values, rec_times, time")
#       }
#
#       required_params <- c("theta", "shape", "scale")
#       if (!all(required_params %in% names(params))) {
#         stop("params must contain: theta, shape, scale")
#       }
#
#       ## --- 2. Extract data ---
#       R <- data$rec_values
#       L <- data$rec_times
#       T <- data$time[1]
#
#       theta <- as.numeric(params["theta"])
#       a     <- as.numeric(params["shape"])   # shape
#       A     <- as.numeric(params["scale"])   # scale (no inversion here)
#
#       m <- length(R)
#
#       if (m == 0 || !is.finite(theta) || !is.finite(a) || !is.finite(A) ||
#           a <= 0 || A <= 0) {
#         return(c(var_theta = -Inf, var_shape = -Inf, var_scale = -Inf))
#       }
#
#       ## --- 3. Transform ---
#       x <- R - theta * L
#       if (any(x <= 0)) {
#         return(c(var_theta = -Inf, var_shape = -Inf, var_scale = -Inf))
#       }
#
#       ## --- 4. Common sums ---
#       s3_theta <- 0
#       s3_A     <- 0
#       s3_a     <- 0
#
#       if (m > 1) {
#         for (j in 1:(m - 1)) {
#
#           start <- L[j] + 1
#           end   <- L[j + 1] - 1
#
#           if (start <= end) {
#             t_seq <- start:end
#             val   <- R[j] - theta * t_seq
#
#             if (any(val <= 0)) {
#               return(c(var_theta = -Inf, var_shape = -Inf, var_scale = -Inf))
#             }
#
#             ## for theta
#             s3_theta <- s3_theta + sum(t_seq^2 * val^(-a - 2))
#
#             ## for A
#             s3_A <- s3_A + sum(val^(-a))
#
#             ## for a
#             term <- A * val
#             s3_a <- s3_a + sum((term^(-a)) * (log(term)^2))
#           }
#         }
#       }
#
#       ## --- tail contribution ---
#       if (L[m] < T) {
#
#         t_seq <- (L[m] + 1):T
#         val   <- R[m] - theta * t_seq
#
#         if (any(val <= 0)) {
#           return(c(var_theta = -Inf, var_shape = -Inf, var_scale = -Inf))
#         }
#
#         s3_theta <- s3_theta + sum(t_seq^2 * val^(-a - 2))
#         s3_A     <- s3_A     + sum(val^(-a))
#
#         term <- A * val
#         s3_a <- s3_a + sum((term^(-a)) * (log(term)^2))
#       }
#
#       ## --- 5. Variance of theta ---
#       s1_theta <- (a + 1) * sum(L^2 / x^2)
#       s2_theta <- sum(L^2 * x^(-a - 2))
#       s4_theta <- -a * (a + 1) * A^(-a) * (s2_theta + s3_theta)
#
#       var_theta <- -1 / (s1_theta + s4_theta)
#
#       ## --- 6. Variance of A (scale) ---
#       s1_A <- m * a / (A^2)
#       s2_A <- sum(x^(-a))
#       s4_A <- -a * (a + 1) * A^(-a - 2) * (s2_A + s3_A)
#
#       var_scale <- -1 / (s1_A + s4_A)
#
#       ## --- 7. Variance of a (shape) ---
#       s1_a <- -m / (a^2)
#       s2_a <- sum((log(A * x))^2 * x^(-a))
#       s4_a <- -A^(-a) * (s2_a + s3_a)
#
#       var_shape <- -1 / (s1_a + s4_a)
#
#       ## --- 8. Return all ---
#       return(c(
#         var_theta = var_theta,
#         var_shape = var_shape,
#         var_scale = var_scale
#       ))
#     }
# )

register_var(model = "LDM",obs_type = "records",dist = "frechet_inv_scale", param_name = "all",   # 👈 unified
             fun <- function(data, params) {

               # --- Checks ---
               if (!all(c("rec_values", "rec_times", "time") %in% names(data))) {
                 stop("data must contain rec_values, rec_times, and time")
               }
               if (!all(c("theta", "shape", "scale") %in% names(params))) {
                 stop("params must contain theta, shape, scale")
               }

               # --- Extract ---
               Rn <- data$rec_values
               Ln <- data$rec_times
               T_val <- data$time[1]

               theta <- params$theta
               alpha <- params$shape
               A <- params$scale

               NT <- length(Rn)

               # --- Build all (n,t) pairs ---
               y_all <- c()
               t_all <- c()

               for (n in 1:(NT-1)) {
                 t_seq <- Ln[n]:(Ln[n+1]-1)
                 y_all <- c(y_all, Rn[n] - theta * t_seq)
                 t_all <- c(t_all, t_seq)
               }

               # last block
               t_seq <- Ln[NT]:T_val
               y_all <- c(y_all, Rn[NT] - theta * t_seq)
               t_all <- c(t_all, t_seq)

               z_all <- A * y_all
               logz <- log(z_all)

               # --- terms at records ---
               y_rec <- Rn - theta * Ln

               # ________________________
               # Hessian components
               # ________________________

               # θθ
               d2_tt <- (1+alpha) * sum(Ln^2 / y_rec^2) -
                 alpha * (alpha+1) * A^2 * sum(t_all^2 * z_all^(-alpha-2))

               # AA
               d2_AA <- alpha * NT / A^2 -
                 alpha*(alpha+1)*sum(y_all^2 * z_all^(-alpha-2))

               # αα
               d2_aa <- -NT / alpha^2 -
                 sum(z_all^(-alpha) * logz^2)

               # θA
               d2_tA <- -alpha * sum(t_all * z_all^(-alpha-1)) +
                 alpha*(alpha+1)*A*sum(t_all * y_all * z_all^(-alpha-2))

               # θα
               d2_ta <- sum(Ln / y_rec) -
                 A * sum(t_all * z_all^(-alpha-1)) +
                 alpha*A*sum(t_all * z_all^(-alpha-1) * logz)

               # Aα
               d2_Aa <- -NT / A +
                 sum(y_all * z_all^(-alpha-1) * (1 - alpha * logz))

               # ________________________
               # Hessian matrix
               # ________________________

               H <- matrix(c(
                 d2_tt, d2_tA, d2_ta,
                 d2_tA, d2_AA, d2_Aa,
                 d2_ta, d2_Aa, d2_aa
               ), nrow = 3, byrow = TRUE)

               colnames(H) <- rownames(H) <- c("theta", "scale", "shape")

               # ____ Fisher Information (observed) ____
               Fisher <- -H
               #vcov_mat = solve(Fisher)
               vcov_mat <- tryCatch({
                 solve(Fisher)
               }, error = function(e) {
                 warning("Matrix singular. Falling back to Moore-Penrose pseudo-inverse.")
                 # ginv handles near-singular matrices smoothly
                 MASS::ginv(Fisher)
               })

               #var_location = -1/d2_mu2, var_scale = -1/d2_sigma2,
               return(list(Hessian = H, vcov = vcov_mat , var_theta = diag(vcov_mat)[1],var_scale = diag(vcov_mat)[2], var_shape = diag(vcov_mat)[3]))
             }
)

## Gumbel
register_var(model = "LDM",obs_type = "records",dist = "gumbel", param_name = "all",   # 👈 unified
  fun <- function(data, params) {

    Rn <- data$rec_values
    Ln <- data$rec_times
    T_val <- data$time[1]

    theta <- params$theta
    A <- params$location
    scale <- params$scale

    NT <- length(Rn)

    # accumulators
    S_e <- 0; S_te <- 0; S_t2e <- 0
    S_ye <- 0; S_tye <- 0; S_y2e <- 0

    sum_t_rec <- sum(Ln)
    sum_y_rec <- sum(Rn - theta * Ln - A)

    for (n in 1:(NT-1)) {
      t_seq <- Ln[n]:(Ln[n+1]-1)
      y <- Rn[n] - theta * t_seq
      yt <- y - A
      u <- yt / scale
      e <- exp(-u)

      S_e <- S_e + sum(e)
      S_te <- S_te + sum(t_seq * e)
      S_t2e <- S_t2e + sum(t_seq^2 * e)
      S_ye <- S_ye + sum(yt * e)
      S_tye <- S_tye + sum(t_seq * yt * e)
      S_y2e <- S_y2e + sum(yt^2 * e)
    }

    # last block
    t_seq <- Ln[NT]:T_val
    y <- Rn[NT] - theta * t_seq
    yt <- y - A
    u <- yt / scale
    e <- exp(-u)

    S_e <- S_e + sum(e)
    S_te <- S_te + sum(t_seq * e)
    S_t2e <- S_t2e + sum(t_seq^2 * e)
    S_ye <- S_ye + sum(yt * e)
    S_tye <- S_tye + sum(t_seq * yt * e)
    S_y2e <- S_y2e + sum(yt^2 * e)

    # ____________________________
    # Hessian
    # ____________________________

    d2_tt <- -S_t2e / scale^2

    d2_AA <- -S_e / scale^2

    d2_aa <- NT/scale^2 -
      2*sum_y_rec/scale^3 +
      2*S_ye/scale^3 -
      S_y2e/scale^4

    d2_tA <- -S_te / scale^2

    d2_ta <- -sum_t_rec/scale^2 +
      S_te/scale^2 -
      S_tye/scale^3

    d2_Aa <- -NT/scale^2 +
      S_e/scale^2 -
      S_ye/scale^3

    H <- matrix(c(
      d2_tt, d2_tA, d2_ta,
      d2_tA, d2_AA, d2_Aa,
      d2_ta, d2_Aa, d2_aa
    ), nrow = 3, byrow = TRUE)

    colnames(H) <- rownames(H) <- c("theta", "location", "scale")

    # ____ Fisher Information (observed) ____
    Fisher <- -H
    #vcov_mat = solve(Fisher)
    vcov_mat <- tryCatch({
      solve(Fisher)
    }, error = function(e) {
      warning("Matrix singular. Falling back to Moore-Penrose pseudo-inverse.")
      # ginv handles near-singular matrices smoothly
      MASS::ginv(Fisher)
    })

    #var_location = -1/d2_mu2, var_scale = -1/d2_sigma2,
    return(list(Hessian = H, vcov = vcov_mat , var_theta = diag(vcov_mat)[1], var_location = diag(vcov_mat)[2], var_scale = diag(vcov_mat)[3]))
  }
)
## YNM, Xt -----------------

  ## frechet
register_var( model = "YNM", obs_type = "all", dist = "frechet_inv_scale", param_name = "gamma",
              fun = function(data, params) {
                if(!is.numeric(data)) stop("data should be a numerical vector")
                n = length(data)
                if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

                gamma <- as.numeric(params$gamma)
                shape <- as.numeric(params$shape)
                scale <- as.numeric(params$scale)

                sum_val <- sum((1:n) * (0:(n-1)) * (gamma^((1:n)-2)) * (scale * data)^(-shape))
                fisher <- (-n * (n + 1) / (2 * gamma^2)) - sum_val
                return(-1 / fisher)

              }
)

register_var( model = "YNM", obs_type = "all", dist = "frechet_inv_scale", param_name = "scale",
              fun = function(data, params) {
                if(!is.numeric(data)) stop("data should be a numerical vector")
                n = length(data)
                if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

                gamma <- params[["gamma"]]
                scale <- params[["scale"]]
                shape <- params[["shape"]]

                sum_val = sum( (gamma^(1:n))*(data^-shape)  )
                fisher =  (shape*n/(scale^2)) -  shape*(shape+1) * scale^(-shape-2)*  sum_val

                return(-1 / fisher)

              }
)

register_var( model = "YNM", obs_type = "all", dist = "frechet_inv_scale", param_name = "shape",
              fun = function(data, params) {
                if(!is.numeric(data)) stop("data should be a numerical vector")
                n = length(data)
                if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

                gamma <- as.numeric(params$gamma)
                shape <- as.numeric(params$shape)
                scale <- as.numeric(params$scale)

                sum_val = sum(  (gamma^(1:n))  *  (log(scale *data))^2    * (scale*data)^(-shape)  )
                fisher =  (-n / (shape^2)) - sum_val
                return(-1 / fisher)
              }
)

register_var( model = "YNM", obs_type = "all", dist = "frechet_inv_scale", param_name = "all",
  fun <- function(data, params) {

    # ___ checks ___
    if (!all(c("gamma", "scale", "shape") %in% names(params))) {
      stop("params must contain gamma, scale, shape")
    }

    x <- data
    if(!is.numeric(x)) stop("data should be a numerical vector")
    T_val <- length(x)
    t <- seq_len(T_val)

    gamma <- params[["gamma"]]
    A <- params[["scale"]]
    alpha <- params[["shape"]]

    # --- base terms ---
    Ax <- A * x
    z <- Ax^(-alpha)
    logAx <- log(Ax)

    w <- gamma^t

    # ______________________________
    # Hessian components
    # ______________________________

    # AA
    d2_AA <- alpha*T_val/A^2 -
      alpha*(alpha+1)*sum(w * x^2 * Ax^(-alpha-2))

    # alpha-alpha
    d2_aa <- -T_val/alpha^2 -
      sum(w * z * logAx^2)

    # gamma-gamma
    d2_gg <- -sum(t/gamma^2) -
      sum(t*(t-1)*gamma^(t-2) * z)

    # A-alpha
    d2_Aa <- -T_val/A +
      sum(w * x * Ax^(-alpha-1) * (1 - alpha*logAx))

    # A-gamma
    d2_Ag <- alpha * sum(t * gamma^(t-1) * x * Ax^(-alpha-1))

    # alpha-gamma
    d2_ag <- sum(t * gamma^(t-1) * z * logAx)

    # ______________________________
    # Hessian matrix
    # ______________________________

    H <- matrix(c(
      d2_gg, d2_Ag, d2_ag,
      d2_Ag, d2_AA, d2_Aa,
      d2_ag, d2_Aa, d2_aa
    ), nrow = 3, byrow = TRUE)

    colnames(H) <- rownames(H) <- c("gamma", "scale", "shape")

    # ____ Fisher Information (observed) ____
    Fisher <- -H
    #vcov_mat = solve(Fisher)
    vcov_mat <- tryCatch({
      solve(Fisher)
    }, error = function(e) {
      warning("Matrix singular. Falling back to Moore-Penrose pseudo-inverse.")
      # ginv handles near-singular matrices smoothly
      MASS::ginv(Fisher)
    })

    #var_location = -1/d2_mu2, var_scale = -1/d2_sigma2,
    return(list(Hessian = H, vcov = vcov_mat ,var_gamma = diag(vcov_mat)[1] ,var_scale = diag(vcov_mat)[2], var_shape = diag(vcov_mat)[3]))
  }
)

# Gumbel
register_var( model = "YNM", obs_type = "all", dist = "gumbel", param_name = "all",
              fun <- function(data, params) {

                x <- data
                if(!is.numeric(x)) stop("data should be a numerical vector")
                T_val <- length(x)
                t <- seq_len(T_val)

                if( all(c("gamma", "location", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, location, scale should be present.")

                gamma <- params[["gamma"]]
                mu <- params[["location"]]
                sigma <- params[["scale"]]

                y <- x - mu
                u <- y / sigma
                e <- exp(-u)
                w <- gamma^t

                # ___ sums___
                S_e   <- sum(w * e)
                S_ye  <- sum(w * y * e)
                S_y2e <- sum(w * y^2 * e)

                sum_y <- sum(y)
                S_t <- sum(t)
                dg_w  <- t * gamma^(t - 1)
                d2g_w <- t * (t - 1) * gamma^(t - 2)

                #______________________________
                # Hessian components
                #______________________________

                d2_mm <- - S_e/sigma^2

                d2_ss <- #T_val/sigma^2 - 2 * sum_y / sigma^3 + 2 * S_ye / sigma^3 - S_y2e / sigma^4
                (T_val -
                    2 * sum(u) +
                    sum(w * u * (2 - u) * e)
                ) / sigma^2

                d2_ms <- (-T_val +sum(w * e * (1 - u)) ) / sigma^2
                #-T_val/sigma^2 +S_e/sigma^2 -S_ye/sigma^3

                d2_gg <- - S_t / gamma^2 - sum(d2g_w * e) #-sum(t/gamma^2) -sum(t*(t-1)*gamma^(t-2) * e)

                d2_gm <- -sum(dg_w * e) / sigma #- sum(t * gamma^(t-1) * e) / sigma

                d2_gs <- -sum(dg_w * u * e) / sigma #- sum(t * gamma^(t-1) * y * e) / sigma^2

                # ___ Hessian ___
                H <- matrix(c(
                  d2_gg, d2_gm, d2_gs,
                  d2_gm, d2_mm, d2_ms,
                  d2_gs, d2_ms, d2_ss
                ), nrow = 3, byrow = TRUE)

                colnames(H) <- rownames(H) <- c("gamma","location","scale")

                # ____ Fisher Information (observed) ____
                Fisher <- -H
                #vcov_mat = solve(Fisher)
                vcov_mat <- tryCatch({
                  solve(Fisher)
                }, error = function(e) {
                  warning("Matrix singular. Falling back to Moore-Penrose pseudo-inverse.")
                  # ginv handles near-singular matrices smoothly
                  MASS::ginv(Fisher)
                })

                #var_location = -1/d2_mu2, var_scale = -1/d2_sigma2,
                return(list(Hessian = H, vcov = vcov_mat , var_gamma = diag(vcov_mat)[1], var_location = diag(vcov_mat)[2], var_scale = diag(vcov_mat)[3]))
              }
)
## YNM, Rn -----------------

## frechet
register_var( model = "YNM", obs_type = "records", dist = "frechet_inv_scale", param_name = "gamma",
              fun = function(data, params) {
                if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                Rn <- data$rec_values
                Ln <- data$rec_times
                n <- data$time[1]

                if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

                v_Rn_Exp_YNM_gamma <- function(a, b, x) {
                  # Terms involving b
                  b_term <- x^b * ((b^2 - 3*b + 2) * x^2 + (4*b - 2*b^2) * x + b^2 - b)

                  # Terms involving a
                  a_term <- x^a * ((-a^2 + 3*a - 2) * x^2 + (2*a^2 - 4*a) * x - a^2 + a)

                  # Denominator
                  denominator <- (x - 1)^3 * x^2

                  # Final result
                  result <- (b_term + a_term) / denominator

                  return(result)
                }

                gamma <- as.numeric(params$gamma)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                m = length(Rn)

                term1 = -sum(Ln)/gamma^2

                term2 = numeric(m)
                for(i in 1:(m-1)){
                  term2[i]=(scale*Rn[i])^(-shape) * v_Rn_Exp_YNM_gamma(a=Ln[i],b=Ln[i+1],x=gamma)
                }
                if(Ln[m]<n){
                  term2[m] = (scale*Rn[m])^(-shape)* v_Rn_Exp_YNM_gamma(a=Ln[m],b=(n+1),x=gamma)}

                fisher= term1-sum(term2)
                return(-1/fisher)
              }
)

register_var( model = "YNM", obs_type = "records", dist = "frechet_inv_scale", param_name = "scale",
              fun = function(data, params) {
                if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                Rn <- data$rec_values
                Ln <- data$rec_times
                n <- data$time[1]

                if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

                gamma <- as.numeric(params$gamma)
                shape <- as.numeric(params$shape)
                scale <- 1/as.numeric(params$scale)

                m = length(Rn)

                term1=shape*m/scale^2

                term2 = numeric(m)
                for(i in 1:(m-1)){
                  term2[i] = (gamma^Ln[i]-gamma^Ln[i+1])*Rn[i]^(-shape)/(1-gamma)
                }

                if(Ln[m]<n){
                  term2[m] = (gamma^Ln[m]-gamma^(n+1))*Rn[m]^(-shape)/(1-gamma) }

                term2=term2*(shape^2+shape) * scale^(-shape-2)

                fisher=term1-sum(term2)

                return(-1/fisher)
              }
)

register_var( model = "YNM", obs_type = "records", dist = "frechet_inv_scale", param_name = "shape",
              fun = function(data, params) {
                if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                Rn <- data$rec_values
                Ln <- data$rec_times
                n <- data$time[1]

                if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

                gamma <- as.numeric(params$gamma)
                shape <- as.numeric(params$shape)
                scale <- as.numeric(params$scale)

                m = length(Rn)

                term1= -m/shape^2

                term2=numeric(m)
                for(i in 1:(m-1)){
                  term2[i] = (gamma^Ln[i]-gamma^Ln[i+1])/(1-gamma)  * (log(scale*Rn[i]))^2 * (scale*Rn[i])^(-shape)
                }

                if(Ln[m]<n){
                  term2[m] = (gamma^Ln[m]-gamma^(n+1))/(1-gamma) * (log(scale*Rn[m]))^2 * (scale*Rn[m])^(-shape) }

                fisher = term1 - sum(term2)

                -NT/shape^2 -
                  sum(w * z * logAx^2)

                return(-1/fisher)
              }
)

register_var(model = "YNM", obs_type = "records", dist = "frechet_inv_scale", param_name = "all",
  fun <- function(data, params) {

      if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

      if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

    # ___ extract ___
    Rn <- data$rec_values
    Ln <- data$rec_times
    T_val <- data$time[1]
    NT <- length(Rn)

    gamma <- params[["gamma"]]
    A <- params[["scale"]]
    shape <- params[["shape"]]

    # ___ base quantities ___
    Ax <- A * Rn
    z <- (Ax)^(-shape)
    logAx <- log(Ax)

    # ___- build weights ___-
    w <- rep(1, NT)

    # (γ^{l_n}-1)
    w <- w + (gamma^Ln - 1)

    # inter-record terms
    for (n in 1:(NT-1)) {
      w[n] <- w[n] +
        (gamma^(Ln[n]+1) - gamma^(Ln[n+1])) / (1 - gamma)
    }

    # last term
    if (Ln[NT] < T_val) {
      w[NT] <- w[NT] +
        (gamma^(Ln[NT]+1) - gamma^(T_val+1)) / (1 - gamma)
    }

    # ____________________________
    # Hessian components
    # ____________________________

    # AA
    d2_AA <- shape*NT/A^2 -
      shape*(shape+1)*sum(w * Rn^2 * Ax^(-shape-2))

    # shape-shape
    d2_aa <- -NT/shape^2 -
      sum(w * z * logAx^2)

    # A-shape
    d2_Aa <- -NT/A +
      sum(w * Rn * Ax^(-shape-1) * (1 - shape*logAx))

    # γ-part (numerical for stability)
    # d2w_gamma <- function(gamma, Ln, T_val) {
    #
    #   NT <- length(Ln)
    #   d2w <- numeric(NT)
    #
    #   for (n in seq_len(NT)) {
    #
    #     # upper limit of geometric sum
    #     if (n < NT) {
    #       b <- Ln[n + 1]
    #     } else if (Ln[n] < T_val) {
    #       b <- T_val + 1
    #     } else {
    #       # no geometric tail contribution
    #       d2w[n] <- Ln[n] * (Ln[n] - 1) * gamma^(Ln[n] - 2)
    #       next
    #     }
    #
    #     A  <- gamma^(Ln[n] + 1) - gamma^b
    #
    #     Ap <- (Ln[n] + 1) * gamma^(Ln[n]) -
    #       b * gamma^(b - 1)
    #
    #     App <- (Ln[n] + 1) * Ln[n] * gamma^(Ln[n] - 1) -
    #       b * (b - 1) * gamma^(b - 2)
    #
    #     d2w[n] <-
    #       Ln[n] * (Ln[n] - 1) * gamma^(Ln[n] - 2) +
    #       (
    #         App * (1 - gamma)^2 +
    #           2 * Ap * (1 - gamma) +
    #           2 * A
    #       ) / (1 - gamma)^3
    #   }
    #
    #   d2w
    # }
    # dw_gamma <- function(gamma, Ln, T_val) {
    #
    #   NT <- length(Ln)
    #   dw <- numeric(NT)
    #
    #   for (n in seq_len(NT)) {
    #
    #     if (n < NT) {
    #       b <- Ln[n + 1]
    #     } else if (Ln[n] < T_val) {
    #       b <- T_val + 1
    #     } else {
    #       dw[n] <- Ln[n] * gamma^(Ln[n] - 1)
    #       next
    #     }
    #
    #     A  <- gamma^(Ln[n] + 1) - gamma^b
    #
    #     Ap <- (Ln[n] + 1) * gamma^(Ln[n]) -
    #       b * gamma^(b - 1)
    #
    #     dw[n] <-
    #       Ln[n] * gamma^(Ln[n] - 1) +
    #       (Ap * (1 - gamma) + A) / (1 - gamma)^2
    #   }
    #
    #   dw
    # }

    # dw = dw_gamma(gamma, Ln, T_val)
    # d2w = d2w_gamma(gamma, Ln, T_val)
    # d2_gg = -sum(Ln) / gamma^2 - sum(Ln * (Ln - 1) * gamma^(Ln - 2))-
    #         sum( - z * d2w)
    #
    # d2_ga <- -sum(e * dw) / sigma
    # d2_gA <- -sum(e * u * dw) / sigma

    logLik_fun_rec <- loglik_registry[["YNM"]][["records"]][["frechet_inv_scale"]]
    eps <- 1e-6
    f_gamma <- function(g) {
      params2 <- params
      params2[["gamma"]] <- g
      logLik_fun_rec(data, params2)
    }

    d2_gg <- (f_gamma(gamma+eps) - 2*f_gamma(gamma) + f_gamma(gamma-eps)) / eps^2

    ## Helper function for d2_gA (gamma and A)
    # f_gamma_A <- function(g, A_val) {
    #   params2 <- params
    #   params2[["gamma"]] <- g
    #   params2[["A"]]     <- A_val
    #   logLik_fun_rec(data, params2)
    # }
    #
    # ## Compute mixed partial for gamma and A
    # d2_gA <- (f_gamma_A(gamma + eps, A + eps) -
    #             f_gamma_A(gamma + eps, A - eps) -
    #             f_gamma_A(gamma - eps, A + eps) +
    #             f_gamma_A(gamma - eps, A - eps)) / (4 * eps^2)
    #
    #
    # ## Helper function for d2_ga (gamma and a)
    # f_gamma_a <- function(g, a_val) {
    #   params2 <- params
    #   params2[["gamma"]] <- g
    #   params2[["a"]]     <- a_val
    #   logLik_fun_rec(data, params2)
    # }
    #
    # ## Compute mixed partial for gamma and a
    # d2_ga <- (f_gamma_a(gamma + eps, shape + eps) -
    #             f_gamma_a(gamma + eps, shape - eps) -
    #             f_gamma_a(gamma - eps, shape + eps) +
    #             f_gamma_a(gamma - eps, shape - eps)) / (4 * eps^2)

    d2_ga <- 0
    d2_gA <- 0

    # ____________________________
    # Hessian matrix
    # ____________________________

    H <- matrix(c(
      d2_gg, d2_gA, d2_ga,
      d2_gA, d2_AA, d2_Aa,
      d2_ga, d2_Aa, d2_aa
    ), nrow = 3, byrow = TRUE)

    colnames(H) <- rownames(H) <- c("gamma", "scale", "shape")

    # ____ Fisher Information (observed) ____
    Fisher <- -H
    #vcov_mat = solve(Fisher)
    vcov_mat <- tryCatch({
      solve(Fisher)
    }, error = function(e) {
      warning("Matrix singular. Falling back to Moore-Penrose pseudo-inverse.")
      # ginv handles near-singular matrices smoothly
      MASS::ginv(Fisher)
    })

    #var_location = -1/d2_mu2, var_scale = -1/d2_sigma2,
    return(list(Hessian = H, vcov = vcov_mat ,var_gamma = diag(vcov_mat)[1] ,var_scale = diag(vcov_mat)[2], var_shape = diag(vcov_mat)[3]))
  }
)

register_var(model = "YNM", obs_type = "records", dist = "gumbel", param_name = "all",
             fun <- function(data, params) {

    if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

    if( all(c("gamma", "location", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, location, scale should be present.")


    # ___ extract ___
    Rn <- data$rec_values
    Ln <- data$rec_times
    T_val <- data$time[1]
    NT <- length(Rn)

    gamma <- params[["gamma"]]
    mu <- params[["location"]]
    sigma <- params[["scale"]]

    # ____ Gumbel terms ____
    y <- Rn - mu
    u <- y / sigma
    e <- exp(-u)

    # ____ weights ____
    w <- rep(1, NT)
    w <- w + (gamma^Ln - 1)

    for (n in 1:(NT-1)) {
      w[n] <- w[n] +
        (gamma^(Ln[n]+1) - gamma^(Ln[n+1])) / (1 - gamma)
    }

    if (Ln[NT] < T_val) {
      w[NT] <- w[NT] +
        (gamma^(Ln[NT]+1) - gamma^(T_val+1)) / (1 - gamma)
    }

    # _____________________________
    # common sums
    # _____________________________

    S_w_e <- sum(w * e)
    S_ye  <- sum(w * y * e)
    S_y2e <- sum(w * y^2 * e)

    sum_y <- sum(y)
    sum_l <- sum(Ln)

    # _____________________________
    # Hessian components
    # _____________________________

    # mu-mu
    d2_mm <-  - S_w_e / sigma^2 #-sum(w * e) / sigma^2

    # sigma-sigma
    d2_ss <- (
      NT - 2 * sum(u) + sum(w * e * u * (2 - u))
    ) / sigma^2

    # mu-sigma
    d2_ms <-  (-NT + sum(w * e * (1 - u))) / sigma^2

    # _____________________________
    # gamma part (exact)
    # _____________________________

    d2w_gamma <- function(gamma, Ln, T_val) {

      NT <- length(Ln)
      d2w <- numeric(NT)

      for (n in seq_len(NT)) {

        # upper limit of geometric sum
        if (n < NT) {
          b <- Ln[n + 1]
        } else if (Ln[n] < T_val) {
          b <- T_val + 1
        } else {
          # no geometric tail contribution
          d2w[n] <- Ln[n] * (Ln[n] - 1) * gamma^(Ln[n] - 2)
          next
        }

        A  <- gamma^(Ln[n] + 1) - gamma^b

        Ap <- (Ln[n] + 1) * gamma^(Ln[n]) -
          b * gamma^(b - 1)

        App <- (Ln[n] + 1) * Ln[n] * gamma^(Ln[n] - 1) -
          b * (b - 1) * gamma^(b - 2)

        d2w[n] <-
          Ln[n] * (Ln[n] - 1) * gamma^(Ln[n] - 2) +
          (
            App * (1 - gamma)^2 +
              2 * Ap * (1 - gamma) +
              2 * A
          ) / (1 - gamma)^3
      }

      d2w
    }
    dw_gamma <- function(gamma, Ln, T_val) {

      NT <- length(Ln)
      dw <- numeric(NT)

      for (n in seq_len(NT)) {

        if (n < NT) {
          b <- Ln[n + 1]
        } else if (Ln[n] < T_val) {
          b <- T_val + 1
        } else {
          dw[n] <- Ln[n] * gamma^(Ln[n] - 1)
          next
        }

        A  <- gamma^(Ln[n] + 1) - gamma^b

        Ap <- (Ln[n] + 1) * gamma^(Ln[n]) -
          b * gamma^(b - 1)

        dw[n] <-
          Ln[n] * gamma^(Ln[n] - 1) +
          (Ap * (1 - gamma) + A) / (1 - gamma)^2
      }

      dw
    }

    dw = dw_gamma(gamma, Ln, T_val)
    d2w = d2w_gamma(gamma, Ln, T_val)
    d2_gg = -sum(Ln) / gamma^2 - sum(e * d2w)

    d2_gm <- -sum(e * dw) / sigma
    d2_gs <- -sum(e * u * dw) / sigma

    # _____________________________
    # gamma part (numerical)
    # _____________________________
    # logLik_fun_rec <- loglik_registry[["YNM"]][["records"]][["gumbel"]]
    # eps <- 1e-15
    # f_gamma <- function(g) {
    #   params2 <- params
    #   params2[["gamma"]] <- g
    #   logLik_fun_rec(data, params2)
    # }
    #
    # d2_gg <- (f_gamma(gamma+eps) - 2*f_gamma(gamma) + f_gamma(gamma-eps)) / eps^2

    # Helper function for d2_gA (gamma and scale s)
    # f_gamma_s <- function(g, s_val) {
    #   params2 <- params
    #   params2[["gamma"]] <- g
    #   params2[["scale"]]     <- s_val
    #   logLik_fun_rec(data, params2)
    # }
    #
    # ## Compute mixed partial for gamma and A
    # d2_gs <- (f_gamma_s(gamma + eps, sigma + eps) -
    #             f_gamma_s(gamma + eps, sigma - eps) -
    #             f_gamma_s(gamma - eps, sigma + eps) +
    #             f_gamma_s(gamma - eps, sigma - eps)) / (4 * eps^2)


    ## Helper function for d2_ga (gamma and location m)
    # f_gamma_m <- function(g, m_val) {
    #   params2 <- params
    #   params2[["gamma"]] <- g
    #   params2[["location"]]     <- m_val
    #   logLik_fun_rec(data, params2)
    # }

    ## Compute mixed partial for gamma and a
    # d2_gm <- (f_gamma_m(gamma + eps, mu + eps) -
    #             f_gamma_m(gamma + eps, mu - eps) -
    #             f_gamma_m(gamma - eps, mu + eps) +
    #             f_gamma_m(gamma - eps, mu - eps)  ) / (4 * eps^2)

    # _____________________________
    # Hessian matrix
    # _____________________________

    H <- matrix(c(
      d2_gg, d2_gm, d2_gs,
      d2_gm, d2_mm, d2_ms,
      d2_gs, d2_ms, d2_ss
    ), nrow = 3, byrow = TRUE)

    colnames(H) <- rownames(H) <- c("gamma", "location", "scale")

    # ____ Fisher Information (observed) ____
    Fisher <- -H
    #vcov_mat = solve(Fisher)
    vcov_mat <- tryCatch({
      solve(Fisher)
    }, error = function(e) {
      warning("Matrix singular. Falling back to Moore-Penrose pseudo-inverse.")
      # ginv handles near-singular matrices smoothly
      MASS::ginv(Fisher)
    })

    #var_location = -1/d2_mu2, var_scale = -1/d2_sigma2,
    return(list(Hessian = H, vcov = vcov_mat ,var_gamma = diag(vcov_mat)[1] ,var_location = diag(vcov_mat)[2], var_scale = diag(vcov_mat)[3]))
  }
)
