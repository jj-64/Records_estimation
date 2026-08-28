logLik_records <- function(model, obs_type, dist, data, params) {
  # lookup
  f <- loglik_registry[[model]][[obs_type]][[dist]]

  if (is.null(f))
    stop("Likelihood expression not registered for this (model, obs_type, dist).")

  f(data, params)
}

## Helper ---------------------------------
# empty container
loglik_registry <- new.env(parent = emptyenv())

## Function to register likelihood expressions later
register_loglik <- function(model, obs_type, dist, fun) {
  if (!exists(model, envir = loglik_registry))
    loglik_registry[[model]] <- list()

  if (!obs_type %in% names(loglik_registry[[model]]))
    loglik_registry[[model]][[obs_type]] <- list()

  loglik_registry[[model]][[obs_type]][[dist]] <- fun
}

## Classical, Xt -------------------------

  ##Normal
register_loglik( model = "iid", obs_type = "all", dist = "norm",
  fun = function(data, params) {
    if(!is.numeric(data)) stop("data should be a numerical vector")
    if( all(c("mean", "sd") %in% names(params)) == FALSE ) stop("parameters mean, sd should be present in a list.")
    mean <- as.numeric(params['mean'])
    sd <- as.numeric(params["sd"])
    sum(dnorm(data, mean =mean , sd = sd, log = TRUE))
  }
)
  ## Frechet
register_loglik( "iid", "all", "frechet",
  fun = function(data, params) {
    if(!is.numeric(data)) stop("data should be a numerical vector")
    if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present in a list.")
    location = ifelse("location"  %in% names(params), as.numeric(params["location"]) ,(min(data)-1e-6))

    sum(VGAM::dfrechet(x = data, location = location, shape = as.numeric(params["shape"]), scale = params["scale"] ,log = TRUE))
  }
)

## Frechet
register_loglik( "iid", "all", "frechet_inv_scale",
                 fun = function(data, params) {
                   if(!is.numeric(data)) stop("data should be a numerical vector")
                   if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present in a list.")
                   location = ifelse("location"  %in% names(params), as.numeric(params["location"]) ,(min(data)-1e-6))

                   T_val = length(data)
                   shape = params["shape"]
                   scale = params["scale"]

                   ll = T_val * log(scale * shape) - (shape +1) * sum( log(scale * data) ) - sum( (scale * data)^(-shape) )
                   return(ll)
                   #sum(VGAM::dfrechet(x = data, location = location, shape = as.numeric(params["shape"]), scale = 1/params["scale"] ,log = TRUE))
                 }
)
  ## Gumbel
register_loglik( "iid", "all", "gumbel",
                 fun = function(data, params) {
                   if(!is.numeric(data)) stop("data should be a numerical vector")
                   if( all(c("location", "scale") %in% names(params)) == FALSE ) stop("parameters location, scale should be present in a list.")

                   sum(VGAM::dgumbel(x = data, location = as.numeric(params["location"]), scale = as.numeric(params["scale"]) ,log = TRUE))
                 }
)

register_loglik( "iid", "all", "gumbel_explicit",
  fun <- function(data, params) {

    # --- Checks ---
    if (!all(c("location", "scale") %in% names(params))) {
      stop("parameters must include location and scale")
    }

    # --- Extract ---
    x = data
    if(!is.numeric(x)) stop("data should be a numerical vector")

    mu <- params["location"]
    sigma <- params["scale"]

    # --- Safety ---
    if (sigma <= 0) return(-Inf)

    # --- Standardized variable ---
    z <- (x - mu) / sigma

    # --- Log-likelihood ---
    ll <- -length(x) * log(sigma) -
      sum(z) -
      sum(exp(-z))

    return(ll)
  }

)
  ## Weibull
register_loglik( "iid", "all", "weibull",
                 fun = function(data, params) {
                   if(!is.numeric(data)) stop("data should be a numerical vector")
                   if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present in a list.")

                   return(sum(dweibull(x = data, shape = as.numeric(params["shape"]), scale = as.numeric(params["scale"]), log= TRUE )))
                 }
)
## Classical, Rn ------------------------

register_loglik("iid", "records", "norm",
  fun = function(data, params) {
    if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

    Rn <- data$rec_values
    Ln <- data$rec_times
    n  <- data$time[1]

    if( all(c("mean", "sd") %in% names(params)) == FALSE ) stop("parameters mean, sd should be present in a list.")

    mean <- as.numeric(params['mean'])
    sd <- as.numeric(params["sd"])
    m = length(Rn) #m= rec_count(y)

    ## sum of log(f_rn)
    s1 <- sum(dnorm(Rn[-1], mean = mean, sd = sd, log = TRUE))
    if (is.nan(s1) || !is.finite(s1)) return(-Inf)
    #Compute s2 using vectorized approach
    intervals <- diff(Ln)-1
    s2b <- intervals * pnorm(Rn[-m], mean = mean, sd = sd, log = TRUE)
    s2 <- sum(s2b)
    if (is.nan(s2) || !is.finite(s2)) return(-Inf)

    ##Compute s3 only if needed
    s3 <- if (Ln[m] < n) {
      (n - Ln[m]) * pnorm(Rn[m], mean = mean, sd = sd, log = TRUE)
    } else {
      0
    }
    if (is.nan(s3) || !is.finite(s3)) return(-Inf)

    #Return total log-likelihood
    return(s1 + s2 + s3)
  }
)

register_loglik("iid", "records", "frechet",
                fun = function(data, params) {
                  if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                  Rn <- data$rec_values
                  Ln <- data$rec_times
                  n  <- data$time[1]
                  m = length(Rn) #m= rec_count(y)

                  if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present in a list.")

                  scale = as.numeric(params["scale"])
                  shape = as.numeric(params["shape"])
                  location = ifelse("location"  %in% names(params), as.numeric(params["location"]) ,(min(Rn)-1e-6))


                  ## sum of log(f_rn)
                  s1 <- sum(VGAM::dfrechet(Rn[-1],  location = location, shape = shape, scale= scale, log = TRUE))
                  if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                  #Compute s2 using vectorized approach
                  intervals <- diff(Ln)-1
                  s2b <- intervals * VGAM::pfrechet(Rn[-m], location = location, shape = shape, scale= scale, log = TRUE)
                  s2 <- sum(s2b)
                  if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                  ##Compute s3 only if needed
                  s3 <- if (Ln[m] < n) {
                    (n - Ln[m]) * VGAM::pfrechet(Rn[m], location = location,  shape = shape, scale= scale, log = TRUE)
                  } else {
                    0
                  }
                  if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                  #Return total log-likelihood
                  return(s1 + s2 + s3)
                }
)

register_loglik("iid", "records", "frechet_inv_scale",
                fun = function(data, params) {
                  if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                  Rn <- data$rec_values
                  Ln <- data$rec_times
                  n  <- data$time[1]
                  m = length(Rn) #m= rec_count(y)

                  if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present in a list.")

                  scale = 1/as.numeric(params["scale"])
                  shape = as.numeric(params["shape"])
                  location = ifelse("location"  %in% names(params), as.numeric(params["location"]) ,(min(Rn)-1e-6))


                  ## sum of log(f_rn)
                  s1 <- sum(VGAM::dfrechet(Rn[-1],  location = location, shape = shape, scale= scale, log = TRUE))
                  if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                  #Compute s2 using vectorized approach
                  intervals <- diff(Ln)-1
                  s2b <- intervals * VGAM::pfrechet(Rn[-m], location = location, shape = shape, scale= scale, log = TRUE)
                  s2 <- sum(s2b)
                  if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                  ##Compute s3 only if needed
                  s3 <- if (Ln[m] < n) {
                    (n - Ln[m]) * VGAM::pfrechet(Rn[m], location = location,  shape = shape, scale= scale, log = TRUE)
                  } else {
                    0
                  }
                  if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                  #Return total log-likelihood
                  return(s1 + s2 + s3)
                }
)

register_loglik("iid", "records", "frechet_inv_scale_explicit",
    fun<- function(data, params) {

      # --- Check inputs ---
      if (!all(c("rec_values", "rec_times", "time") %in% names(data))) {
        stop("data must contain rec_values, rec_times, and time")
      }

      if (!all(c("scale", "shape") %in% names(params))) {
        stop("params must contain scale and shape")
      }

      # --- Extract data ---
      x  <- data$rec_values     # R_n values
      l  <- data$rec_times      # l_n times
      Tt <- data$time[1]        # T
      m <- length(x)           # number of records

      scale     <- as.numeric(params["scale"])
      shape <- as.numeric(params["shape"])

      # --- Basic checks ---
      if (scale <= 0 || shape <= 0) return(-Inf)

      # --- Term 1 ---
      term1 <- (m - 1) * log(scale * shape)

      # --- Term 2 ---
      term2 <- -(1 + shape) * sum(log(scale * x[-1]))

      # --- Term 3 ---
      term3 <- -sum((scale * x[-1])^(-shape))

      # --- Term 4 ---
      intervals <- diff(l) - 1
      term4 <- -sum(intervals * (scale * x[-m])^(-shape))

      # --- Term 5 ---
      term5 <- -(Tt - l[m]) * (scale * x[m])^(-shape)

      # --- Total log-likelihood ---
      ll <- term1 + term2 + term3 + term4 + term5

      # --- Safety check ---
      if (is.nan(ll) || !is.finite(ll)) return(-Inf)

      return(ll)
    }
)


register_loglik("iid", "records", "gumbel",
                fun = function(data, params) {
                  if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                  Rn <- data$rec_values
                  Ln <- data$rec_times
                  n  <- data$time[1]

                  if( all(c("location", "scale") %in% names(params)) == FALSE ) stop("parameters location, scale should be present in a list.")

                  scale = as.numeric(params["scale"])
                  location = as.numeric(params["location"])
                  m = length(Rn) #m= rec_count(y)

                  ## sum of log(f_rn)
                  s1 <- sum(VGAM::dgumbel(Rn[-1], location, scale, log = TRUE))
                  if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                  #Compute s2 using vectorized approach
                  intervals <- diff(Ln)-1
                  s2b <- intervals * VGAM::pgumbel(Rn[-m], location, scale, log = TRUE)
                  s2 <- sum(s2b)
                  if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                  ##Compute s3 only if needed
                  s3 <- if (Ln[m] < n) {
                    (n - Ln[m]) * VGAM::pgumbel(Rn[m], location, scale, log = TRUE)
                  } else {
                    0
                  }
                  if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                  #Return total log-likelihood
                  return(s1 + s2 + s3)
                }
)

register_loglik("iid", "records", "gumbel_explicit",
  fun <- function(params, data) {

    if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

    if( all(c("location", "scale") %in% names(params)) == FALSE ) stop("parameters location, scale should be present in a list.")

    location <- params[["location"]]
    sigma <-  params[["scale"]]

    # --- Data ---
    Rn <- data$rec_values
    Ln <- data$rec_times
    T_val <- data$time[1]

    m <- length(Rn)

    x_next <- Rn[-1]
    x_curr <- Rn[-m]
    delta  <- diff(Ln)
    last_weight <- T_val - Ln[m]

    # Standardized variables
    z_next <- (x_next - location) / sigma
    z_curr <- (x_curr - location) / sigma
    z_last <- (Rn[m] - location) / sigma

    # Log-likelihood
    ll <- -(m - 1)*log(sigma) -
      sum(z_next) -
      sum(exp(-z_next)) -
      sum((delta - 1) * exp(-z_curr)) -
      last_weight * exp(-z_last)

    return(ll)  # NEGATIVE for minimization
  }
)


register_loglik("iid", "records", "weibull",
                fun = function(data, params) {
                  if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                  Rn <- data$rec_values
                  Ln <- data$rec_times
                  n  <- data$time[1]

                  if( all(c("shape", "scale") %in% names(params)) == FALSE ) stop("parameters shape, scale should be present in a list.")

                  scale = as.numeric(params["scale"])
                  shape = as.numeric(params["shape"])
                  m = length(Rn) #m= rec_count(y)

                  ## sum of log(f_rn)
                  s1 <- sum(dweibull(Rn[-1], shape, scale, log = TRUE))
                  if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                  #Compute s2 using vectorized approach
                  intervals <- diff(Ln)-1
                  s2b <- intervals * pweibull(Rn[-m], shape, scale, log = TRUE)
                  s2 <- sum(s2b)
                  if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                  ##Compute s3 only if needed
                  s3 <- if (Ln[m] < n) {
                    (n - Ln[m]) * pweibull(Rn[m], shape, scale, log = TRUE)
                  } else {
                    0
                  }
                  if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                  #Return total log-likelihood
                  return(s1 + s2 + s3)
                }
)


## DTRW, Xt ---------------------------------

  ## Normal
register_loglik("DTRW", "all", "norm",
  fun = function(data, params) {
    if(!is.numeric(data)) stop("data should be a numerical vector")

    if( all(c("mean", "sd") %in% names(params)) == FALSE ) stop("parameters mean, sd should be present in a list.")

    mean <- as.numeric(params['mean'])
    sd <- as.numeric(params["sd"])

    sum(dnorm(x = diff(data), mean =mean , sd = sd, log = TRUE))
  }
)

  ##cauchy
register_loglik("DTRW", "all", "cauchy",
                fun = function(data, params) {
                  if(!is.numeric(data)) stop("data should be a numerical vector")
                  if( all(c("loc", "scale") %in% names(params)) == FALSE ) stop("parameters loc, scale should be present in a list.")

                  sum(dcauchy(x = diff(data), location =as.numeric(params["loc"]) , scale = as.numeric(params["scale"]) , log = TRUE))
                }
)

## uniform
register_loglik("DTRW", "all", "unif",
                fun = function(data, params) {
                  if(!is.numeric(data)) stop("data should be a numerical vector")

                  if( all(c("max") %in% names(params)) == FALSE ) stop("parameters min, max should be present in a list.")

                  max <- as.numeric(params["max"])
                  min <- -max

                  sum(dunif(x = diff(data), min = min , max = max, log = TRUE))
                }
)

## DTRW, Rn -------------------------------------
  ## Normal
register_loglik( "DTRW", "records", "norm",
  fun = function(data, params) {
    if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

    Rn <- data$rec_values
    Ln <- data$rec_times
    n <- data$time[1]
    if( all(c( "sd") %in% names(params)) == FALSE ) stop("parameters sd should be present in a list.")

    mean = 0
    variance = as.numeric(params["sd"])^2

    # Example: probability that S_1, S_2, ..., S_5 all < 0
    # sigma_cov = function(k, sigma){  ## k should be less than 1000, sigma is variance
    #   x=matrix(0, nrow=k, ncol=k)
    #   for(i in 1:k){
    #     for(j in 1:k){
    #       x[i,j] = min(i,j)*sigma
    #     }
    #   }
    #   return(x)
    # }

    m= length(Rn)  # Number of observed time points

    if (m < 2) {
      stop("The vector l must contain at least two indices.")
    }

    s1b = dnorm(diff(Rn), mean = 0, sd = sqrt(diff(Ln)*variance), log = TRUE )
    s1=sum(s1b)
    if (is.nan(s1) || !is.finite(s1)) return(-Inf)

    ## case where m<T
    s2=0
    if (Ln[m] < n) {
      #s2=log(mvtnorm::pmvnorm(lower=-Inf, upper=rep(0,(n - Ln[m])), sigma = sigma_cov((n - Ln[m]),params))[1])
      s2 = -0.5*log(pi * (n - Ln[m]))
    }
    if (is.nan(s2) || !is.finite(s2)) return(-Inf)

    return(s1+s2)
  }
)

  ## Cauchy
register_loglik( "DTRW", "records", "cauchy",
                 fun = function(data, params) {
                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n <- data$time[1]

                   if( all(c("scale") %in% names(params)) == FALSE ) stop("parameters  scale should be present in a list.")

                   loc = 0
                   scale = as.numeric(params["scale"])
                   m= length(Rn)  # Number of observed time points

                   if (m < 2) {
                     stop("The vector l must contain at least two indices.")
                   }

                   s1b = dcauchy(diff(Rn), location = 0, scale = sqrt(diff(Ln)*scale), log = TRUE)
                   s1=sum(s1b)
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)
                   ## case where m<T
                   s2=0
                   if (Ln[m] < n) { s2 = -0.5*log(pi * (n - Ln[m])) }
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                  ## Return
                   return(s1+s2)
                 }
)

## Uniform
register_loglik( "DTRW", "records", "unif",
                 fun = function(data, params) {
                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n <- data$time[1]

                   if( all(c("max") %in% names(params)) == FALSE ) stop("parameters max should be present in a list.")

                   max = as.numeric(params["max"])
                   min = - max

                   m= length(Rn)  # Number of observed time points

                   if (m < 2) {
                     stop("The vector l must contain at least two indices.")
                   }

                   s1b = dunif(diff(Rn), min = min, max = max, log = TRUE)
                   s1=sum(s1b)
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)
                   ## case where m<T
                   s2=0
                   if (Ln[m] < n) { s2 = -0.5*log(pi * (n - Ln[m])) }
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   ## Return
                   return(s1+s2)
                 }
)
## YNM, Xt -----------------------------------------

  ##Frechet
register_loglik( "YNM", "all", "frechet",
                 fun = function(data, params) {
                   x=data
                   if(!is.numeric(x)) stop("data should be a numerical vector")

                   n= length(as.numeric(x))

                   if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

                   gamma <- as.numeric(params["gamma"])
                   shape <- as.numeric(params["shape"])
                   scale <- as.numeric(params["scale"])
                   location = ifelse("location"  %in% names(params), as.numeric(params["location"]) ,(min(x)-1e-6))

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (shape <=0 | scale <= 0 | gamma <= 1) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {VGAM::dfrechet(x, location = par$location, shape=par$shape, scale=par$scale, log = TRUE)}
                   cdf=function(x,par) {VGAM::pfrechet(x, location = par$location, shape=par$shape, scale=par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 =  sum( (1:n) * log(gamma))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2 = sum( (-1+gamma^(1:n))  * cdf(x, par=list(location = location, shape= shape, scale= scale))  )
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   s3 = sum(pdf(x, par=list(location = location, shape= shape, scale= scale)))
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)
                   # s1 <- n* log(shape * scale^(-shape) )
                   # s2 <- -(shape + 1) * sum(log(x))
                   # s3 <- (n* (n+ 1) / 2) * log(gamma)
                   # s4 <- -sum( gamma^(1:n) * (scale * x)^(-shape) )
                   return(s1+s2+s3)
                 }
)

register_loglik( "YNM", "all", "frechet_inv_scale",
                fun <- function(data, params) {

                   # ___ checks ___

                   if (!all(c("gamma", "scale", "shape") %in% names(params))) {
                     stop("params must contain gamma, scale, shape")
                   }

                   # ___ extract ___
                   x <- data
                   if(!is.numeric(x)) stop("data should be a numerical vector")

                   T_val <- length(x)
                   t <- seq_len(T_val)

                   gamma <- as.numeric(params[["gamma"]])
                   A <- as.numeric(params[["scale"]])
                   alpha <- as.numeric(params[["shape"]])

                   # ___ basic checks ___
                   if (A <= 0 || alpha <= 0 || gamma <= 1) {
                     return(-Inf)
                   }

                   # ___ compute terms ___
                   Ax <- A * x
                   logAx <- log(Ax)
                   z <- Ax^(-alpha)         # (A x_t)^(-alpha)
                   w <- gamma^t            # gamma^t

                   # ___ log-likelihood ___
                   ll <- sum(t * log(gamma)) +
                     T_val * log(A * alpha) -
                     (1 + alpha) * sum(logAx) -
                     sum(w * z)

                   return(ll)
                 }
)

##Gumbel

register_loglik( "YNM", "all", "gumbel",
                 fun = function(data, params) {
                   x=data
                   if(!is.numeric(x)) stop("data should be a numerical vector")

                   n= length(as.numeric(x))
                   if( all(c("gamma", "loc", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, loc, scale should be present.")
                   gamma <- as.numeric(params["gamma"])
                   loc <- as.numeric(params["loc"])
                   scale <- as.numeric(params["scale"])
                   tvec = seq_len(n)

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or loc <= 0, we return -Inf)
                   if (scale <= 0 || gamma <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {VGAM::dgumbel(x, location=par$loc, scale=par$scale, log = TRUE)}
                   cdf=function(x,par) {VGAM::pgumbel(x, location=par$loc, scale=par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum( tvec * log(gamma) )
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2 = sum( (gamma^tvec - 1) * cdf(x, par=list(loc= loc, scale= scale) ) )
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   s3 = sum(pdf(x, par=list(loc= loc, scale= scale)))
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   return(s1+s2+s3)
                 }
)

register_loglik( "YNM", "all", "gumbel_explicit",
                 fun <- function(data, params) {

                   x <- data
                   if(!is.numeric(x)) stop("data should be a numerical vector")

                   T_val <- length(x)
                   t <- seq_len(T_val)

                   gamma <- params[["gamma"]]
                   mu <- params[["location"]]
                   sigma <- params[["scale"]]

                   if (sigma <= 0 || gamma <= 0) return(-Inf)

                   u <- (x - mu) / sigma
                   e <- exp(-u)
                   w <- gamma^t

                   ll <- sum(t * log(gamma)) -
                     T_val * log(sigma) -
                     sum(u) -
                     sum(w * e)

                   return(ll)
                 }
)
  ##norm
register_loglik( "YNM", "all", "norm",
                 fun = function(data, params) {
                   x=data
                   if(!is.numeric(x)) stop("data should be a numerical vector")

                   n= length(as.numeric(x))
                   if( all(c("gamma", "mean", "sd") %in% names(params)) == FALSE ) stop("parameters gamma, mean, sd should be present.")
                   gamma <- as.numeric(params["gamma"])
                   mean <- as.numeric(params["mean"])
                   sd <- as.numeric(params["sd"])

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or mean <= 0, we return -Inf)
                   if (sd <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {dnorm(x, mean=par$mean, sd=par$sd, log = TRUE)}
                   cdf=function(x,par) {pnorm(x, mean=par$mean, sd=par$sd, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum( (1:n) * log(gamma))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2 = sum((-1+gamma^(1:n)) * cdf(x, par=list(mean= mean, sd= sd)) )
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   s3 = sum( pdf(x, par=list(mean= mean, sd= sd)) )
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)


                   return(s1+s2+s3)
                 }
)

##weibull
register_loglik( "YNM", "all", "weibull",
                 fun = function(data, params) {
                   x=data
                   if(!is.numeric(x)) stop("data should be a numerical vector")

                   n= length(as.numeric(x))
                   if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, scale, shape should be present.")
                   gamma <- as.numeric(params["gamma"])
                   shape <- as.numeric(params["shape"])
                   scale <- as.numeric(params["scale"])

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (shape <=0 || scale <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {dweibull(x, shape=par$shape, scale=par$scale, log = TRUE)}
                   cdf=function(x,par) {pweibull(x, shape=par$shape, scale=par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   tvec = seq_len(n)
                   s1 = sum( tvec * log(gamma) )
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2 = sum( (gamma^tvec - 1) * cdf(x, par=list(shape= shape, scale= scale) ) )
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   s3 = sum(pdf(x, par=list(shape= shape, scale= scale)))
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   return(s1+s2+s3)
                 }
)
## YNM, Rn -----------------------------------------

  ## Frechet
register_loglik( "YNM", "records", "frechet",
  fun = function(data, params) {
    if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

    Rn <- data$rec_values
    Ln <- data$rec_times
    n <- data$time[1]

    if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

    gamma <- as.numeric(params["gamma"])
    shape <- as.numeric(params["shape"])
    scale <- as.numeric(params["scale"])
    m = length(Rn) #m= rec_count(y)  ## number of records
    location = ifelse("location"  %in% names(params), as.numeric(params["location"]) ,(min(Rn)-1e-6))

    # Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
    if (scale <= 0 || shape<=0 || gamma <0) {
      return(-Inf)  # Invalid parameters, return a large negative value
    }

    ## pdf
    pdf=function(x,par) {VGAM::dfrechet(x=x, location = par$location, shape=par$shape, scale=par$scale, log = TRUE)}
    ## cdf
    cdf=function(x,par) {VGAM::pfrechet(q=x, location = par$location, shape= par$shape, scale=par$scale, log = TRUE)}


    # your exact YNM record-pair likelihood
    s1 = sum(Ln) *log(gamma)
    if (is.nan(s1) || !is.finite(s1)) return(-Inf)

    s2 = sum( pdf(Rn, par = list(location = location, shape=shape, scale=scale)) + (-1+gamma^Ln) * cdf(Rn, par = list(location = location, shape=shape, scale=scale)))
    if (is.nan(s2) || !is.finite(s2)) return(-Inf)
    s3b = 0
    for( i in 1:(m-1)  ){
      if((Ln[i]+1 <= Ln[i+1]) == TRUE){ ## we have non-records in between
        s3b[i]= sum( ( (gamma^(Ln[i]+1) - gamma^(Ln[i+1]))/(1-gamma)) * cdf(Rn[i], par = list(location = location, shape=shape, scale=scale)) )
      }
    }
    s3 = sum(na.omit(s3b))
    if (is.nan(s3) || !is.finite(s3)) return(-Inf)

    s4=0
    if (Ln[m] < n) {
      s4 = sum( ((gamma^(Ln[m]+1) - gamma^(n+1))/(1-gamma) ) * cdf(Rn[m], par = list(location = location, shape=shape, scale=scale) ))
    }
    if (is.nan(s4) || !is.finite(s4)) return(-Inf)

     return(s1+s2+s3+s4)
    }
)

register_loglik( "YNM", "records", "frechet_inv_scale",
  fun <- function(data, params) {

    # --- checks ---
    if (!all(c("rec_values", "rec_times", "time") %in% names(data))) {
      stop("data must contain rec_values, rec_times, time")
    }
    if (!all(c("gamma", "scale", "shape") %in% names(params))) {
      stop("params must contain gamma, scale, shape")
    }

    Rn <- data$rec_values
    Ln <- data$rec_times
    T_val <- data$time[1]

    gamma <- params[["gamma"]]
    A <- params[["scale"]]
    alpha <- params[["shape"]]

    NT <- length(Rn)

    Ax <- A * Rn
    z <- Ax^(-alpha)
    logAx <- log(Ax)

    # --- initialize weights ---
    w <- rep(1, NT)

    # (γ^{l_n} - 1)
    w <- w + (gamma^Ln - 1)

    # inter-record terms
    for (n in 1:(NT-1)) {
      w[n] <- w[n] +
        (gamma^(Ln[n] + 1) - gamma^(Ln[n+1])) / (1 - gamma)
    }

    # last term
    if (Ln[NT] < T_val) {
      w[NT] <- w[NT] +
        (gamma^(Ln[NT] + 1) - gamma^(T_val + 1)) / (1 - gamma)
    }

    # --- log-likelihood ---
    ll <- sum(Ln) * log(gamma) +
      NT * log(A * alpha) -
      (1 + alpha) * sum(logAx) -
      sum(w * z)

    return(ll)
  }
)
  ## Gumbel
register_loglik( "YNM", "records", "gumbel",
  fun = function(data, params) {
                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n = data$time[1]
                   if( all(c("gamma", "location", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, location, scale should be present.")

                   gamma <- as.numeric(params["gamma"])
                   location <- as.numeric(params["location"])
                   scale <- as.numeric(params["scale"])
                   m = length(Rn) #m= rec_count(y)  ## number of records

                   # Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (scale <= 0) {
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {VGAM::dgumbel(x=x, location=par$location, scale=par$scale, log = TRUE)}
                   ## cdf
                   cdf=function(x,par) {VGAM::pgumbel(q=x, location=par$location, scale=par$scale, log = TRUE)}


                   # your exact YNM record-pair likelihood
                   s1 = sum(Ln) * log(gamma)
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2 = sum(pdf(Rn, par = list(location = location, scale=scale)) + (-1+gamma^Ln) * cdf(Rn, par = list(location=location, scale=scale)))
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   s3b = 0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s3b[i]= sum(( (gamma^(Ln[i]+1) - gamma^(Ln[i+1]))/(1-gamma)) * cdf(Rn[i], par = list(location = location, scale=scale)) )
                     }
                   }
                   s3 = sum(na.omit(s3b))
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   s4=0
                   if (Ln[m] < n) {
                     s4 = sum( ( (gamma^(Ln[m]+1) - gamma^(n+1))/(1-gamma) ) * cdf(Rn[m], par = list(location = location, scale=scale)) )
                   }
                   if (is.nan(s4) || !is.finite(s4)) return(-Inf)

                   return(s1+s2+s3+s4)
                 }
)

register_loglik( "YNM", "records", "gumbel_explicit",
  fun = function(data, params) {

    # --- checks ---
    if (!all(c("rec_values", "rec_times", "time") %in% names(data))) {
      stop("data must contain rec_values, rec_times, time")
    }
    if (!all(c("gamma", "location", "scale") %in% names(params))) {
      stop("params must contain gamma, location, scale")
    }


    Rn <- data$rec_values
    Ln <- data$rec_times
    T_val <- data$time[1]

    gamma <- as.numeric(params["gamma"])
    mu <- as.numeric(params["location"])
    sigma <- as.numeric(params["scale"])

    NT <- length(Rn)

    # --- Gumbel terms ---
    u <- (Rn - mu) / sigma
    e <- exp(-u)

    # --- weights ---
    w <- rep(1, NT)

    # gamma^{l_n} - 1
    w <- w + (gamma^Ln - 1)

    # inter-record contributions
    for (n in 1:(NT - 1)) {
      w[n] <- w[n] +
        (gamma^(Ln[n] + 1) - gamma^(Ln[n + 1])) / (1 - gamma)
    }

    # last term
    if (Ln[NT] < T_val) {
      w[NT] <- w[NT] +
        (gamma^(Ln[NT] + 1) - gamma^(T_val + 1)) / (1 - gamma)
    }

    # --- log-likelihood ---
    ll <- sum(Ln) * log(gamma) -
      NT * log(sigma) -
      sum(u) -
      sum(w * e)

    return(ll)
  }
)

  ## Norm
register_loglik( "YNM", "records", "norm",
                 fun = function(data, params) {
                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n = data$time[1]
                   if( all(c("gamma", "mean", "sd") %in% names(params)) == FALSE ) stop("parameters gamma, mean, sd should be present.")

                   gamma <- as.numeric(params["gamma"])
                   mean <- as.numeric(params["mean"])
                   sd <- as.numeric(params["sd"])
                   m = length(Rn) #m= rec_count(y)  ## number of records

                   # Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (sd <= 0) {
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {dnorm(x=x, mean=par$mean, sd=par$sd, log = TRUE)}
                   ## cdf
                   cdf=function(x,par) {pnorm(q=x, mean=par$mean, sd=par$sd, log = TRUE)}


                   # your exact YNM record-pair likelihood
                   s1 = sum(Ln) * log(gamma)
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2 = sum(pdf(Rn, par = list(mean = mean, sd=sd)) + (-1+gamma^Ln) *  cdf(Rn, par = list(mean=mean, sd=sd)))
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   s3b = 0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s3b[i]= sum( ( (gamma^(Ln[i]+1) - gamma^(Ln[i+1]-1))/(1-gamma)) * cdf(Rn[i], par = list(mean = mean, sd=sd)) )
                     }
                   }
                   s3 = sum(na.omit(s3b))
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   s4=0
                   if (Ln[m] < n) {
                     s4 = sum(( (gamma^(Ln[m]+1) - gamma^(n+1))/(1-gamma) ) * cdf(Rn[m], par = list(mean = mean, sd=sd)) )
                   }
                   if (is.nan(s4) || !is.finite(s4)) return(-Inf)

                   return(s1+s2+s3+s4)
                 }
)

  ## weibull
register_loglik( "YNM", "records", "weibull",
                 fun = function(data, params) {
                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n <- data$time[1]
                   if( all(c("gamma", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters gamma, shape, scale should be present.")

                   gamma <- as.numeric(params["gamma"])
                   shape <- as.numeric(params["shape"])
                   scale <- as.numeric(params["scale"])
                   m = length(Rn) #m= rec_count(y)  ## number of records

                   # Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (scale <= 0 | shape<=0) {
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {dweibull(x=x, shape=par$shape, scale=par$scale, log = TRUE)}
                   ## cdf
                   cdf=function(x,par) {pweibull(q=x, shape= par$shape, scale=par$scale, log = TRUE)}


                   # your exact YNM record-pair likelihood
                   s1 = log(gamma^(sum(Ln)))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2 = sum( pdf(Rn, par = list(shape=shape, scale=scale)) + (-1 + gamma^Ln) * cdf(Rn, par = list(shape=shape, scale=scale)) )
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   s3b = 0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]) == TRUE){ ## we have non-records in between
                       s3b[i]= sum( ( (gamma^(Ln[i]+1) - gamma^(Ln[i+1]))/(1-gamma)) * cdf(Rn[i], par = list(shape=shape, scale=scale) ) )
                     }
                   }
                   s3 = sum(na.omit(s3b))
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   s4=0
                   if (Ln[m] < n) {
                     s4 = sum( ((gamma^(Ln[m]+1) - gamma^(n+1))/(1-gamma) ) * cdf(Rn[m], par = list(shape=shape, scale=scale)) )
                   }
                   if (is.nan(s4) || !is.finite(s4)) return(-Inf)

                   return(s1+s2+s3+s4)
                 }
)

## exp
register_loglik( "YNM", "records", "exp",
                 fun = function(data, params) {
                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n <- data$time[1]
                   if( all(c("gamma", "rate") %in% names(params)) == FALSE ) stop("parameters gamma, rate should be present.")

                   gamma <- as.numeric(params["gamma"])
                   rate <- as.numeric(params["rate"])

                   m = length(Rn) #m= rec_count(y)  ## number of records

                   # Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (rate <=0) {
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## pdf
                   pdf=function(x,par) {dexp(x=x, rate=par$rate , log = TRUE)}
                   ## cdf
                   cdf=function(x,par) {pexp(q=x, rate=par$rate, log = TRUE)}


                   # your exact YNM record-pair likelihood
                   s1 = log(gamma^(sum(Ln)))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2 = sum( pdf(Rn, par = list(rate = rate)) + (-1 + gamma^Ln) * cdf(Rn, par = list(rate = rate) ) )
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   s3b = 0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]) == TRUE){ ## we have non-records in between
                       s3b[i]= sum( ( (gamma^(Ln[i]+1) - gamma^(Ln[i+1]))/(1-gamma)) * cdf(Rn[i], par = list(rate = rate) ) )
                     }
                   }
                   s3 = sum(na.omit(s3b))
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   s4=0
                   if (Ln[m] < n) {
                     s4 = sum( ((gamma^(Ln[m]+1) - gamma^(n+1))/(1-gamma) ) * cdf(Rn[m], par = list(rate = rate)) )
                   }
                   if (is.nan(s4) || !is.finite(s4)) return(-Inf)

                   return(s1+s2+s3+s4)
                 }
)
## LDM, Xt -----------------------------------------

  ##Norm
register_loglik( "LDM", "all", "norm",
                 fun = function(data, params) {
                   y=data
                   if(!is.numeric(y)) stop("data should be a numerical vector")

                   if( all(c("theta", "mean", "sd") %in% names(params)) == FALSE ) stop("parameters theta, mean, sd should be present.")

                   theta <- as.numeric(params["theta"])
                   mean <- as.numeric(as.numeric(params["mean"]))
                   sd <- as.numeric(params["sd"])

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if ( sd <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## transform
                   x = y - theta * (1:length(y)) ## y-theta*t

                   ## pdf
                   pdf=function(x,par) {dnorm(x=x, mean=par$mean, sd=par$sd, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(pdf(x=x, par=list(mean = mean, sd = sd)))

                   return(s1)
                 }
)
  ##Frechet
register_loglik( "LDM", "all", "frechet",
                 fun = function(data, params) {
                   y=data
                   if(!is.numeric(y)) stop("data should be a numerical vector")

                   if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")

                   theta <- as.numeric(params["theta"])
                   shape <- as.numeric(params["shape"])
                   scale <- as.numeric(params["scale"])

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (!is.finite(theta) || !is.finite(shape) || !is.finite(scale) ||
                       shape <= 0 || scale <= 0) {
                     return(-Inf)
                   }

                   ## transform
                   x = y - theta * (1:length(y)) ## y-theta*t
                   location = ifelse("location"  %in% names(params), as.numeric(params["location"]) ,(min(x)-1e-6))

                   ## pdf
                   pdf=function(x,par) {VGAM::dfrechet(x=x, location = par$location, shape=par$shape, scale=par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(pdf(x=x, par=list(location = location, shape= shape, scale= scale)))

                   # scale = 1/params["scale"]
                   # if (any(x <= 0)) return(-1000)
                   # A <- (sum(x^-scale) / n)^(1/scale)
                   # s1 <- n * log(scale * A^(-scale))
                   # s2 <- -(scale + 1) * sum(log(x))
                   # s3 <- -sum((A * x)^(-scale))
                   # return (s1 + s2 + s3)
                   return(s1)
                 }
)

##Frechet inversted scale
register_loglik( "LDM", "all", "frechet_inv_scale",
                 fun = function(data, params) {
                   y=data
                   if(!is.numeric(y)) stop("data should be a numerical vector")

                   if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")

                   theta <- as.numeric(params["theta"])
                   shape <- as.numeric(params["shape"])
                   scale <- as.numeric(params["scale"])

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (!is.finite(theta) || !is.finite(shape) || !is.finite(scale) ||
                       shape <= 0 || scale <= 0) {
                     return(-Inf)
                   }

                   ## 4. Transform
                   x = y - theta * (1:length(y)) ## y-theta*t
                   location = ifelse("location"  %in% names(params), as.numeric(params["location"]) ,0 )
                                     #(min(x)-1e-6))

                   ## pdf
                   pdf=function(x,par) {VGAM::dfrechet(x=x, location = par$location, shape=par$shape, scale= 1/par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(pdf(x=x, par=list(location = location, shape= shape, scale= scale)))

                   # scale = 1/params["scale"]
                   # if (any(x <= 0)) return(-1000)
                   # A <- (sum(x^-scale) / n)^(1/scale)
                   # s1 <- n * log(scale * A^(-scale))
                   # s2 <- -(scale + 1) * sum(log(x))
                   # s3 <- -sum((A * x)^(-scale))
                   # return (s1 + s2 + s3)
                   return(s1)
                 }
)
  ##gumbel
register_loglik( "LDM", "all", "gumbel",
                 fun = function(data, params) {
                   y=data
                   if(!is.numeric(y)) stop("data should be a numerical vector")

                   if( all(c("theta", "location", "scale") %in% names(params)) == FALSE ) stop("parameters theta, location, scale should be present.")
                   theta <- as.numeric(params["theta"])
                   location <- as.numeric(params["location"])
                   scale <- as.numeric(params["scale"])

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (scale <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## transform
                   x = y - theta * (1:length(y)) ## y-theta*t

                   ## pdf
                   pdf=function(x,par) {VGAM::dgumbel(x=x, location=par$location, scale=par$scale ,log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(pdf(x=x, par=list(location = location, scale= scale)))
                   return(s1)
                 }
)

register_loglik( "LDM", "all", "gumbel_explicit",
                 fun = function(data, params) {
                   y=data
                   if(!is.numeric(y)) stop("data should be a numerical vector")

                   if( all(c("theta", "location", "scale") %in% names(params)) == FALSE ) stop("parameters theta, location, scale should be present.")
                   theta <- as.numeric(params["theta"])
                   location <- as.numeric(params["location"])
                   scale <- as.numeric(params["scale"])
                   T_val = length(y)

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (scale <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## transform
                   x_shifted = y - theta * (1:length(y)) ## y-theta*t
                   x_scaled = (x_shifted - location)/scale

                   ## pdf
                   ll = - T_val * log(scale) - sum(x_scaled) - sum(exp(- x_scaled))

                   return(ll)
                 }
)
##weibull
register_loglik( "LDM", "all", "weibull",
                 fun = function(data, params) {
                   y=data
                   if(!is.numeric(y)) stop("data should be a numerical vector")

                   if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")
                   theta <- as.numeric(params["theta"])
                   shape <- as.numeric(params["shape"])
                   scale <- as.numeric(params["scale"])

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if (scale <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## transform
                   x = y - theta * (1:length(y)) ## y-theta*t

                   ## pdf
                   pdf=function(x,par) {dweibull(x=x, shape=par$shape, scale=par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(pdf(x=x, par=list(shape = shape, scale= scale)))
                   return(s1)
                 }
)
## LDM, Rn -----------------------------------------

    ##Gumbel
register_loglik( "LDM", "records", "gumbel",
                 fun = function(data, params) {

                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n <- data$time[1]

                   if( all(c("theta", "location", "scale") %in% names(params)) == FALSE ) stop("parameters theta, location, scale should be present.")
                   theta <- as.numeric(params["theta"])
                   location <- as.numeric(params["location"])
                   scale <- as.numeric(params["scale"])
                   m = length(Rn)

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if ( scale <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## transform
                   x = Rn - theta * Ln ## Rn-theta*ln

                   ## pdf
                   pdf=function(x,par) {VGAM::dgumbel(x=x, location=par$location, scale=par$scale , log = TRUE)}
                   ## cdf
                   cdf=function(x,par) {VGAM::pgumbel(q=x, location=par$location, scale=par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(pdf(x=x, par=list(location = location, scale = scale)))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2b=0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s2b[i]=sum(cdf( Rn[i]-theta * (Ln[i]+1):(Ln[i+1]-1), par=list(location = location, scale = scale) )  )
                     }
                   }
                   s2 = sum(na.omit(s2b))
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   ## case where last record is not last observation
                   s3=0
                   if( (Ln[m] < n) == TRUE  ) {
                     s3a=cdf( Rn[m]-theta * (Ln[m]+1):n, par=list(location = location, scale = scale) )
                     s3 = sum(s3a[is.finite(s3a)])
                   }
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   return(s1+s2+s3)
                 }
)

register_loglik( "LDM", "records", "gumbel_explicit",
            fun <- function(data, params) {

    # --- checks ---
    if (!all(c("rec_values", "rec_times", "time") %in% names(data))) {
      stop("data must contain rec_values, rec_times, time")
    }
    if (!all(c("theta", "location", "scale") %in% names(params))) {
      stop("params must contain theta, location, scale")
    }

    Rn <- data$rec_values
    Ln <- data$rec_times
    T_val <- data$time[1]

    theta <- as.numeric(params["theta"])
    A <- as.numeric(params["location"])
    scale <- as.numeric(params["scale"])

    NT <- length(Rn)

    # log(f)
    sum_linear <- 0
    sum_exp <- 0

    y <- Rn - theta * Ln
    u <- (y - A) / scale
    sum_exp <- sum(exp(-u))
    sum_linear <- sum(u)

    s1 = NT * log(scale) + sum_linear + sum_exp

    #log(F)
    sum_exp <- 0

    for (n in 1:(NT-1)) {
      if((Ln[n]+1 <= Ln[n+1]-1) == TRUE){
        t_seq <- (Ln[n]+1):(Ln[n+1]-1)
        y <- Rn[n] - theta * t_seq
        u <- (y - A) / scale

        sum_exp = sum_exp + sum(exp(-u))
     }
      }

    s2 = sum_exp

    # last block
    s3 <- 0

    if( (Ln[NT] < T_val) == TRUE  ) {
      t_seq <- (Ln[NT] +1 ):T_val
      y <- Rn[NT] - theta * t_seq
      u <- (y - A) / scale

      s3 =  sum(exp(-u))
    }


    loglik <- - (s1 + s2 + s3)

    return(loglik)
  }
)

    ##normal
register_loglik( "LDM", "records", "norm",
                 fun = function(data, params) {
                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n <- data$time[1]

                   if( all(c("theta", "mean", "sd") %in% names(params)) == FALSE ) stop("parameters theta, mean, sd should be present.")
                   theta <- as.numeric(params["theta"])
                   mean <- as.numeric(params["mean"])
                   sd <- as.numeric(params["sd"])
                   m = length(Rn)

                   ## Check for invalid parameter values
                   if (sd <= 0) {
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## transform
                   x = Rn - theta * Ln ## Rn-theta*ln

                   ## pdf
                   pdf=function(x,par) {dnorm(x=x, mean=par$mean, sd=par$sd, log = TRUE)}
                   ## cdf
                   cdf=function(x,par) {pnorm(q=x, mean=par$mean, sd=par$sd, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum(pdf(x=x, par=list(mean = mean, sd = sd)))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2b=0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s2b[i]=sum(cdf( Rn[i]-theta * (Ln[i]+1):(Ln[i+1]-1), par=list(mean = mean, sd = sd) )  )
                     }
                   }
                   s2 = sum(na.omit(s2b))
                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   ## case where last record is not last observation
                   s3=0
                   if( (Ln[m] < n) == TRUE  ) {
                     s3a= cdf( Rn[m]-theta * (Ln[m]+1):n, par=list(mean = mean, sd = sd) )
                     s3 = sum(s3a[is.finite(s3a)])
                     }
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   return(s1+s2+s3)
                 }
)

##frechet
register_loglik( "LDM", "records", "frechet",
                 fun = function(data, params) {
                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n <- data$time[1]

                   if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")
                   theta <- as.numeric(params["theta"])
                   shape <- as.numeric(params["shape"])
                   scale <- as.numeric(params["scale"])

                   m = length(Rn)

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if ( shape <0 |scale <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## transform
                   x = Rn - theta * Ln ## Rn-theta*ln
                   location = ifelse("location"  %in% names(params), as.numeric(params["location"]) ,(min(x)-1e-6))

                   ## pdf
                   pdf=function(x,par) {VGAM::dfrechet(x=x, location = par$location, shape=par$shape, scale=par$scale, log = TRUE)}

                   ## cdf
                   cdf=function(x,par) {VGAM::pfrechet(q=x, location = par$location, shape=par$shape, scale=par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum((pdf(x=x, par=list(location = location, shape = shape, scale = scale))))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2b=0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s2b[i]=sum(  cdf( Rn[i]- theta * (Ln[i]+1):(Ln[i+1]-1), par=list(location= location, shape = shape, scale = scale) )  )
                     }
                   }
                   s2b = s2b[is.finite(s2b)]
                   s2 = sum(na.omit(s2b))

                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   ## case where last record is not last observation
                   s3=0
                   if( (Ln[m] < n) == TRUE  ) {
                     s3a=(cdf( Rn[m]-theta * (Ln[m]+1):n, par=list(location = location , shape = shape, scale = scale) )  )
                     s3 = sum(s3a[is.finite(s3a)])
                   }
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   return(s1+s2+s3)
                 }
)
## frechet explicit
register_loglik( "LDM", "records", "frechet_explicit",
                 fun <- function(data, params) {


                   ## 1. Input checks

                   required_data <- c("rec_values", "rec_times", "time")
                   if (!all(required_data %in% names(data))) {
                     stop("data must contain: rec_values, rec_times, time")
                   }

                   required_params <- c("theta", "shape", "scale")
                   if (!all(required_params %in% names(params))) {
                     stop("params must contain: theta, shape, scale")
                   }

                   ## 2. Extract inputs

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n  <- data$time[1]

                   theta <- as.numeric(params["theta"])
                   shape <- as.numeric(params["shape"])
                   scale <- as.numeric(params["scale"])

                   #theta_max <- min(Rn / Ln)
                   #if (theta >= theta_max) return(-1e12)

                   m <- length(Rn)

                   ## 3. Validate parameters

                   if (!is.finite(theta) || !is.finite(shape) || !is.finite(scale) ||
                       shape <= 0 || scale <= 0) {
                     return(-Inf)
                   }

                   if (m == 0) return(-Inf)

                   ## 4. Transform
                   x <- Rn - theta * Ln

                   location <- if ("location" %in% names(params)) {
                     as.numeric(params["location"])
                   } else {
                     #min(x) - 1e-6
                     0
                   }

                   # ---- Check validity ----
                   if (any(Rn - theta * Ln <= 0, na.rm = TRUE)) {
                     return(-Inf)   # invalid parameter region
                   }


                   # --- First term ---
                   term1 <- m * log(shape * scale^(-shape))

                   # --- Second term ---
                   term2 <- -(shape + 1) * sum(log(Rn - theta * Ln))

                   # --- Third term ---
                   term3 <- -scale^(-shape) * sum((Rn - theta * Ln)^(-shape))

                   # --- Fourth term ---
                   term4 <- numeric(m)  # better initialization

                   if (m > 1) {
                     for (j in 1:(m - 1)) {

                       # Only proceed if valid interval exists
                       if (Ln[j] + 1 <= Ln[j + 1] - 1) {

                         t_seq <- (Ln[j] + 1):(Ln[j + 1] - 1)
                         val <- Rn[j] - theta * t_seq

                         # optional safety check
                         # if (any(val <= 0)) return(-Inf)

                         term4[j] <- -sum(scale^(-shape) * val^(-shape))
                       }
                     }
                   }

                   term4 <- sum(term4, na.rm = TRUE)

                   # --- Fifth term ---
                   term5 <- 0
                   if (Ln[m] < n) {
                     t_seq <- seq(Ln[m] + 1, n)
                     val <- Rn[m] - theta * t_seq

                     if (any(val <= 0)) return(-Inf)

                     term5 <- -sum(scale^(-shape) * val^(-shape))
                   }

                   # --- Final log-likelihood ---
                   ll <- term1 + term2 + term3 + term4 + term5

                   if (!is.finite(ll)) return(-Inf)

                   return(ll)
                 }

)

##frechet inverted scale
register_loglik( "LDM", "records", "frechet_inv_scale",
                 fun <- function(data, params) {

                   ## 1. Input checks

                   required_data <- c("rec_values", "rec_times", "time")
                   if (!all(required_data %in% names(data))) {
                     stop("data must contain: rec_values, rec_times, time")
                   }

                   required_params <- c("theta", "shape", "scale")
                   if (!all(required_params %in% names(params))) {
                     stop("params must contain: theta, shape, scale")
                   }

                   ## 2. Extract inputs

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n  <- data$time[1]

                   theta <- as.numeric(params["theta"])
                   shape <- as.numeric(params["shape"])
                   scale <- as.numeric(params["scale"])

                   m <- length(Rn)

                   ## 3. Validate parameters

                   if (!is.finite(theta) || !is.finite(shape) || !is.finite(scale) ||
                       shape <= 0 || scale <= 0) {
                     return(-Inf)
                   }

                   if (m == 0) return(-Inf)


                   ## 4. Transform
                   x <- Rn - theta * Ln

                   location <- if ("location" %in% names(params)) {
                     as.numeric(params["location"])
                   } else {
                     0# min(x) - 1e-6
                   }

                   par_list <- list(location = location, shape = shape, scale = scale)

                   ## 5. PDF and CDF (log-scale)

                   pdf <- function(x) {
                     VGAM::dfrechet(x, location = par_list$location,
                                    shape = par_list$shape,
                                    scale = 1/par_list$scale, log = TRUE)
                   }

                   cdf <- function(x) {
                     VGAM::pfrechet(x, location = par_list$location,
                                    shape = par_list$shape,
                                    scale = 1/par_list$scale, log = TRUE)
                   }


                   ## 6. Likelihood components

                   ## s1: observed records
                   s1 <- sum(pdf(x)[is.finite(pdf(x))]) #sum(pdf(x), na.rm = TRUE)

                   if (!is.finite(s1)) return(-Inf)

                   ## s2: between-record intervals (vectorized)
                   s2 <- 0
                   if (m > 1) {
                     for (i in seq_len(m - 1)) {
                       gap_start <- Ln[i] + 1
                       gap_end   <- Ln[i + 1] - 1

                       if (gap_start <= gap_end) {
                         t_seq <- gap_start:gap_end
                         vals <- Rn[i] - theta * t_seq
                         contrib <- cdf(vals)
                         s2 <- s2 + sum(contrib[is.finite(contrib)])
                       }
                     }
                   }
                   if (!is.finite(s2)) return(-Inf)

                   ## s3: tail (after last record)
                   s3 <- 0
                   if (Ln[m] < n) {
                     t_seq <- (Ln[m] + 1):n
                     vals <- Rn[m] - theta * t_seq
                     contrib <- cdf(vals)
                     s3 <- sum(contrib[is.finite(contrib)])
                   }
                   if (!is.finite(s3)) return(-Inf)

                   ## 7. Return total log-likelihood

                   ll <- s1 + s2 + s3

                   if (!is.finite(ll)) return(-Inf)

                   return(ll)
                 }

)

## frechet inverted scale-explicit
register_loglik( "LDM", "records", "frechet_inv_scale_explicit",
                 fun <- function(data, params) {


                   ## 1. Input checks

                   required_data <- c("rec_values", "rec_times", "time")
                   if (!all(required_data %in% names(data))) {
                     stop("data must contain: rec_values, rec_times, time")
                   }

                   required_params <- c("theta", "shape", "scale")
                   if (!all(required_params %in% names(params))) {
                     stop("params must contain: theta, shape, scale")
                   }

                   ## 2. Extract inputs

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n  <- data$time[1]

                   theta <- as.numeric(params["theta"])
                   shape <- as.numeric(params["shape"])
                   scale <- as.numeric(params["scale"])

                   #theta_max <- min(Rn / Ln)
                   #if (theta >= theta_max) return(-1e12)

                   m <- length(Rn)

                   ## 3. Validate parameters

                   if (!is.finite(theta) || !is.finite(shape) || !is.finite(scale) ||
                       shape <= 0 || scale <= 0) {
                     return(-Inf)
                   }

                   if (m == 0) return(-Inf)

                   ## 4. Transform
                   x <- Rn - theta * Ln

                   location <- if ("location" %in% names(params)) {
                     as.numeric(params["location"])
                   } else {
                    0# min(x) - 1e-6
                   }

                    # ---- Check validity ----
                     if (any(Rn - theta * Ln <= 0, na.rm = TRUE)) {
                       return(-Inf)   # invalid parameter region
                     }


                     # --- First term ---
                     term1 <- m * log(shape * scale^(-shape))

                     # --- Second term ---
                     term2 <- -(shape + 1) * sum(log(Rn - theta * Ln))

                     # --- Third term ---
                     term3 <- -scale^(-shape) * sum((Rn - theta * Ln)^(-shape))

                     # --- Fourth term ---
                     # term4 <- 0
                     # if (m >= 1) {
                     #   for (j in 1:(m - 1)) {
                     #     t_seq <- seq(Ln[j] + 1, Ln[j + 1] - 1)
                     #     if (length(t_seq) > 0) {
                     #       val <- Rn[j] - theta * t_seq
                     #
                     #       # avoid invalid values
                     #       #if (any(val <= 0)) return(-1e12)
                     #
                     #       term4[j] <- -sum(scale^(-shape) * val^(-shape))
                     #     }
                     #   }
                     # }
                     term4 <- numeric(m)  # better initialization

                     if (m > 1) {
                       for (j in 1:(m - 1)) {

                         # Only proceed if valid interval exists
                         if (Ln[j] + 1 <= Ln[j + 1] - 1) {

                           t_seq <- (Ln[j] + 1):(Ln[j + 1] - 1)
                           val <- Rn[j] - theta * t_seq

                           # optional safety check
                           # if (any(val <= 0)) return(-Inf)

                           term4[j] <- -sum(scale^(-shape) * val^(-shape))
                         }
                       }
                     }

                     term4 <- sum(term4, na.rm = TRUE)

                     # --- Fifth term ---
                     term5 <- 0
                     if (Ln[m] < n) {
                       t_seq <- seq(Ln[m] + 1, n)
                       val <- Rn[m] - theta * t_seq

                       if (any(val <= 0)) return(-Inf)

                       term5 <- -sum(scale^(-shape) * val^(-shape))
                     }

                     # --- Final log-likelihood ---
                     ll <- term1 + term2 + term3 + term4 + term5

                   if (!is.finite(ll)) return(-Inf)

                   return(ll)
                 }

)

##Weibull
register_loglik( "LDM", "records", "weibull",
                 fun = function(data, params) {
                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n <- data$time[1]

                   if( all(c("theta", "shape", "scale") %in% names(params)) == FALSE ) stop("parameters theta, shape, scale should be present.")
                   theta <- as.numeric(params["theta"])
                   shape <- as.numeric(params["shape"])
                   scale <- as.numeric(params["scale"])

                   m = length(Rn)

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if ( shape <0 |scale <= 0) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## transform
                   x = Rn - theta * Ln ## Rn-theta*ln

                   ## pdf
                   pdf=function(x,par) {dweibull(x=x,  shape=par$shape, scale=par$scale, log = TRUE)}

                   ## cdf
                   cdf=function(x,par) {pweibull(q=x,  shape=par$shape, scale=par$scale, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum((pdf(x=x, par=list( shape = shape, scale = scale))))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2b=0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s2b[i]=sum(  cdf( Rn[i]- theta * (Ln[i]+1):(Ln[i+1]-1), par=list( shape = shape, scale = scale) )  )
                     }
                   }
                   s2b = s2b[is.finite(s2b)]
                   s2 = sum(na.omit(s2b))

                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   ## case where last record is not last observation
                   s3=0
                   if( (Ln[m] < n) == TRUE  ) {
                     s3a=(cdf( Rn[m]-theta * (Ln[m]+1):n, par=list( shape = shape, scale = scale) )  )
                     s3 = sum(s3a[is.finite(s3a)])
                   }
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   return(s1+s2+s3)
                 }
)

##exponential
register_loglik( "LDM", "records", "exp",
                 fun = function(data, params) {
                   if( all(c("rec_values", "rec_times", "time") %in% names(data)) == FALSE ) stop("a list of rec_values, rec_times, and time should be present.")

                   Rn <- data$rec_values
                   Ln <- data$rec_times
                   n <- data$time[1]

                   if( all(c("theta", "rate") %in% names(params)) == FALSE ) stop("parameters theta, rate should be present.")
                   theta <- as.numeric(params["theta"])
                   rate <- as.numeric(params["rate"])

                   m = length(Rn)

                   ## Check for invalid parameter values (e.g., if params[1] <= 0 or params[2] <= 0 or shape <= 0, we return -Inf)
                   if ( rate <0 ) { ##
                     return(-Inf)  # Invalid parameters, return a large negative value
                   }

                   ## transform
                   x = Rn - theta * Ln ## Rn-theta*ln

                   ## pdf
                   pdf=function(x,par) {dexp(x=x,  rate = par$rate, log = TRUE)}

                   ## cdf
                   cdf=function(x,par) {pexp(q=x,  rate = par$rate, log = TRUE)}

                   # Compute the terms of the likelihood function, while checking for potential issues like negative log arguments
                   s1 = sum((pdf(x=x, par=list( rate = rate))))
                   if (is.nan(s1) || !is.finite(s1)) return(-Inf)

                   s2b=0
                   for( i in 1:(m-1)  ){
                     if((Ln[i]+1 <= Ln[i+1]-1) == TRUE){ ## we have non-records in between
                       s2b[i]=sum(  cdf( Rn[i]- theta * (Ln[i]+1):(Ln[i+1]-1), par=list( rate = rate) )  )
                     }
                   }
                   s2b = s2b[is.finite(s2b)]
                   s2 = sum(na.omit(s2b))

                   if (is.nan(s2) || !is.finite(s2)) return(-Inf)

                   ## case where last record is not last observation
                   s3=0
                   if( (Ln[m] < n) == TRUE  ) {
                     s3a=(cdf( Rn[m]-theta * (Ln[m]+1):n, par=list( rate = rate) )  )
                     s3 = sum(s3a[is.finite(s3a)])
                   }
                   if (is.nan(s3) || !is.finite(s3)) return(-Inf)

                   return(s1+s2+s3)
                 }
)
