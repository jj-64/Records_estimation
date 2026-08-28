############################################################
# Robust optimization for H0 likelihood
# 1. Scale parameters
# 2. Multiple random starts
# 3. Select minimum negative log-likelihood
# 4. Inspect parameter boundaries
# 5. Verify with likelihood profile / grid
############################################################
T_val = 100

true_params <- c(
  trend =0.2,
  rate  = 1
)

series_H0 <- function(T_val, trend, par_H0){

  # Exponential Linear Drift Model
  LDM_series(
    T     = T_val,
    dist  = "exp",
    theta = trend,
    rate  = par_H0[1]
  )
}

Likelihood_under_H0 <- function(data_rec, params){

  logLik_fun_rec <-
    loglik_registry[["LDM"]][["records"]][["exp"]]

  logLik_fun_rec(
    data   = data_rec,
    params = params
  )
}

xt <- series_H0(
  T      = T_val,
  trend  = true_params["trend"],
  par_H0 = true_params["rate"]
)

R <- rec_values(xt)
L <- rec_times(xt)

## Reject degenerate samples with only one record
while(length(R) <= 1){

  xt <- series_H0(
    T      = T_val,
    trend  = true_params["trend"],
    par_H0 = true_params["rate"]
  )

  R <- rec_values(xt)
  L <- rec_times(xt)
}

data_rec <- list(
  rec_values = R,
  rec_times  = L,
  time       = T_val
)


## ---------------------------------------------------------
## 0. Parameter bounds
## ---------------------------------------------------------

lb_0 <- c(
  theta = 0.01,
  rate  = 0.01
)

ub_0 <- c(
  theta = min(R / L),
  rate  = 10
)


## ---------------------------------------------------------
## 1. Parameter scaling
##
## z = (parameter - lower) / (upper - lower)
## Therefore z is constrained to [0,1].
## ---------------------------------------------------------

scale_params <- function(z, lower, upper) {

  params <- lower + z * (upper - lower)

  names(params) <- names(lower)

  params
}


unscale_params <- function(params, lower, upper) {

  z <- (params - lower) / (upper - lower)

  names(z) <- names(params)

  z
}


## ---------------------------------------------------------
## 2. Negative log-likelihood on scaled parameter space
## ---------------------------------------------------------

negLogLik_scaled <- function(z) {

  # Transform scaled parameters back to original scale
  params <- scale_params(
    z = z,
    lower = lb_0,
    upper = ub_0
  )

  # Evaluate likelihood
  ll <- tryCatch(
    Likelihood_under_H0(
      data_rec,
      params
    ),
    error = function(e) NA_real_
  )

  # Protect optimizer against invalid likelihood values
  if (!is.finite(ll)) {
    return(1e100)
  }

  return(-ll)
}


## ---------------------------------------------------------
## 3. Multiple random starts
## ---------------------------------------------------------

set.seed(123)

n_starts <- 50

starts <- lapply(
  seq_len(n_starts),
  function(i) {

    c(
      theta = runif(1, 0, 1),
      rate  = runif(1, 0, 1)
    )

  }
)


## ---------------------------------------------------------
## 4. Run nlminb from every starting value
## ---------------------------------------------------------

fits <- vector(
  mode = "list",
  length = n_starts
)

for (i in seq_len(n_starts)) {

  cat(
    "\nOptimization start:",
    i,
    "of",
    n_starts,
    "\n"
  )

  fit_i <- tryCatch(

    nlminb(
      start = starts[[i]],
      objective = negLogLik_scaled,

      lower = c(
        theta = 0,
        rate  = 0
      ),

      upper = c(
        theta = 1,
        rate  = 1
      ),

      control = list(
        eval.max = 2000,
        iter.max = 1000,
        rel.tol  = 1e-10
      )
    ),

    error = function(e) {

      cat(
        "Optimization error:",
        conditionMessage(e),
        "\n"
      )

      NULL

    }
  )

  fits[[i]] <- fit_i
}


## ---------------------------------------------------------
## 5. Extract valid solutions
## ---------------------------------------------------------

valid <- vapply(
  fits,
  function(fit) {

    !is.null(fit) &&
      is.finite(fit$objective) &&
      all(is.finite(fit$par))

  },
  logical(1)
)


fits_valid <- fits[valid]


if (length(fits_valid) == 0) {

  stop(
    "No valid optimization solution was obtained."
  )

}


## ---------------------------------------------------------
## 6. Select solution with smallest negative log-likelihood
## ---------------------------------------------------------

objectives <- sapply(
  fits_valid,
  function(fit) fit$objective
)

best_index <- which.min(objectives)

fit_best <- fits_valid[[best_index]]


## ---------------------------------------------------------
## 7. Transform best solution back to original parameters
## ---------------------------------------------------------

params_hat <- scale_params(
  z = fit_best$par,
  lower = lb_0,
  upper = ub_0
)


theta_hat <- params_hat["theta"]
rate_hat  <- params_hat["rate"]


## ---------------------------------------------------------
## 8. Report optimization result
## ---------------------------------------------------------

cat("\n============================================\n",
    "BEST OPTIMIZATION RESULT\n",
    "============================================\n",
  "theta_hat =",
  format(theta_hat, digits = 8),
  "\n",
  "rate_hat  =",
  format(rate_hat, digits = 8),
  "\n",
  "Negative log-likelihood =",
  format(fit_best$objective, digits = 10),
  "\n",
  "Log-likelihood =",
  format(-fit_best$objective, digits = 10),
  "\n",
  "Convergence code =",
  fit_best$convergence,
  "\n",
  "Message =",
  fit_best$message,
  "\n",
  "============================================\n")


## ---------------------------------------------------------
## 9. Compare all optimization starts
## ---------------------------------------------------------

results_starts <- data.frame(
  start = seq_along(fits_valid),
  theta_scaled = sapply(
    fits_valid,
    function(fit) fit$par["theta"]
  ),
  rate_scaled = sapply(
    fits_valid,
    function(fit) fit$par["rate"]
  ),
  objective = objectives,
  convergence = sapply(
    fits_valid,
    function(fit) fit$convergence
  )
)


results_starts <- results_starts[
  order(results_starts$objective),
]


cat("\n============================================\n",
    "MULTIPLE-START RESULTS\n",
    "============================================\n")

print(
  head(
    results_starts,
    10
  )
)


## ---------------------------------------------------------
## 10. Transform all solutions to original parameter scale
## ---------------------------------------------------------

results_starts$theta <- lb_0["theta"] +
  results_starts$theta_scaled *
  (ub_0["theta"] - lb_0["theta"])

results_starts$rate <- lb_0["rate"] +
  results_starts$rate_scaled *
  (ub_0["rate"] - lb_0["rate"])


## ---------------------------------------------------------
## 11. Boundary inspection
## ---------------------------------------------------------

boundary_tol <- 1e-6

theta_at_lower <- abs(
  theta_hat - lb_0["theta"]
) <= boundary_tol

theta_at_upper <- abs(
  theta_hat - ub_0["theta"]
) <= boundary_tol

rate_at_lower <- abs(
  rate_hat - lb_0["rate"]
) <= boundary_tol

rate_at_upper <- abs(
  rate_hat - ub_0["rate"]
) <= boundary_tol


cat("\n============================================\n",
    "BOUNDARY INSPECTION\n",
    "============================================\n",
  "theta lower bound:",
  lb_0["theta"],
  "\n",
  "theta estimate:",
  theta_hat,
  "\n",
  "theta upper bound:",
  ub_0["theta"],
  "\n",
  "theta at lower boundary:",
  theta_at_lower,
  "\n",
  "theta at upper boundary:",
  theta_at_upper,
  "\n\n",
  "rate lower bound:",
  lb_0["rate"],
  "\n",
  "rate estimate:",
  rate_hat,
  "\n",
  "rate upper bound:",
  ub_0["rate"],
  "\n",
  "rate at lower boundary:",
  rate_at_lower,
  "\n",
  "rate at upper boundary:",
  rate_at_upper,
  "\n"
)


## ---------------------------------------------------------
## 12. Check whether multiple starts agree
## ---------------------------------------------------------

theta_sd <- sd(
  results_starts$theta
)

rate_sd <- sd(
  results_starts$rate
)


cat("\n============================================\n",
    "MULTIPLE-START STABILITY\n",
    "============================================\n",
  "SD of theta estimates =",
  theta_sd,
  "\n",
  "SD of rate estimates  =",
  rate_sd,
  "\n"
)


## ---------------------------------------------------------
## 13. Likelihood grid
## ---------------------------------------------------------

n_grid <- 100

theta_grid <- seq(
  lb_0["theta"],
  ub_0["theta"],
  length.out = n_grid
)

rate_grid <- seq(
  lb_0["rate"],
  ub_0["rate"],
  length.out = n_grid
)


ll_grid <- matrix(
  NA_real_,
  nrow = n_grid,
  ncol = n_grid
)


for (i in seq_along(theta_grid)) {

  for (j in seq_along(rate_grid)) {

    params <- c(
      theta = theta_grid[i],
      rate  = rate_grid[j]
    )

    ll <- tryCatch(

      Likelihood_under_H0(
        data_rec,
        params
      ),

      error = function(e) NA_real_

    )

    if (is.finite(ll)) {

      ll_grid[i, j] <- ll

    }

  }
}


## ---------------------------------------------------------
## 14. Convert to negative log-likelihood surface
## ---------------------------------------------------------

nll_grid <- -ll_grid


## ---------------------------------------------------------
## 15. Locate grid maximum
## ---------------------------------------------------------

grid_index <- which(
  ll_grid == max(
    ll_grid,
    na.rm = TRUE
  ),
  arr.ind = TRUE
)


theta_grid_max <- theta_grid[
  grid_index[1, "row"]
]

rate_grid_max <- rate_grid[
  grid_index[1, "col"]
]


cat("\n============================================\n",
    "GRID SEARCH RESULT\n",
    "============================================\n",
  "Grid theta maximum =",
  theta_grid_max,
  "\n",
  "Grid rate maximum =",
  rate_grid_max,
  "\n",
  "Optimized theta =",
  theta_hat,
  "\n",
  "Optimized rate =",
  rate_hat,
  "\n"
)


## ---------------------------------------------------------
## 16. Likelihood surface plot
## ---------------------------------------------------------

filled.contour(
  theta_grid,
  rate_grid,
  t(ll_grid),

  xlab = expression(theta),
  ylab = "rate",

  main = "Log-Likelihood Surface",

  plot.axes = {

    axis(1)
    axis(2)

    points(
      theta_hat,
      rate_hat,
      pch = 19,
      cex = 1.3
    )

    points(
      theta_grid_max,
      rate_grid_max,
      pch = 4,
      lwd = 2,
      cex = 1.5
    )

  }
)


## ---------------------------------------------------------
## 17. Likelihood profile for theta
##
## For each theta, optimize rate.
## ---------------------------------------------------------

profile_theta <- numeric(
  length(theta_grid)
)


rate_profile_theta <- numeric(
  length(theta_grid)
)


for (i in seq_along(theta_grid)) {

  theta_i <- theta_grid[i]

  objective_rate <- function(rate_scaled) {

    params <- list(
      theta = theta_i,
      rate = lb_0["rate"] +
        rate_scaled *
        (ub_0["rate"] - lb_0["rate"])
    )

    ll <- tryCatch(

      Likelihood_under_H0(
        data_rec,
        params
      ),

      error = function(e) NA_real_

    )

    if (!is.finite(ll)) {

      return(1e100)

    }

    -ll

  }


  fit_rate <- tryCatch(

    nlminb(
      start = 0.5,
      objective = objective_rate,
      lower = 0,
      upper = 1
    ),

    error = function(e) NULL

  )


  if (!is.null(fit_rate)) {

    profile_theta[i] <- -fit_rate$objective

    rate_profile_theta[i] <-
      lb_0["rate"] +
      fit_rate$par *
      (ub_0["rate"] - lb_0["rate"])

  } else {

    profile_theta[i] <- NA_real_

  }

}


## ---------------------------------------------------------
## 18. Plot profile likelihood for theta
## ---------------------------------------------------------

plot(
  theta_grid,
  profile_theta,
  type = "l",
  lwd = 2,

  xlab = expression(theta),
  ylab = "Profile log-likelihood",

  main = expression(
    "Profile log-likelihood for " ~ theta
  )
)

abline(
  v = theta_hat,
  lty = 2
)

points(
  theta_hat,
  max(
    profile_theta,
    na.rm = TRUE
  ),
  pch = 19
)


## ---------------------------------------------------------
## 19. Likelihood profile for rate
##
## For each rate, optimize theta.
## ---------------------------------------------------------

profile_rate <- numeric(
  length(rate_grid)
)


theta_profile_rate <- numeric(
  length(rate_grid)
)


for (j in seq_along(rate_grid)) {

  rate_j <- rate_grid[j]

  objective_theta <- function(theta_scaled) {

    params <- list(
      theta = lb_0["theta"] +
        theta_scaled *
        (ub_0["theta"] - lb_0["theta"]),

      rate = rate_j
    )

    ll <- tryCatch(

      Likelihood_under_H0(
        data_rec,
        params
      ),

      error = function(e) NA_real_

    )

    if (!is.finite(ll)) {

      return(1e100)

    }

    -ll

  }


  fit_theta <- tryCatch(

    nlminb(
      start = 0.5,
      objective = objective_theta,
      lower = 0,
      upper = 1
    ),

    error = function(e) NULL

  )


  if (!is.null(fit_theta)) {

    profile_rate[j] <- -fit_theta$objective

    theta_profile_rate[j] <-
      lb_0["theta"] +
      fit_theta$par *
      (ub_0["theta"] - lb_0["theta"])

  } else {

    profile_rate[j] <- NA_real_

  }

}


## ---------------------------------------------------------
## 20. Plot profile likelihood for rate
## ---------------------------------------------------------

plot(
  rate_grid,
  profile_rate,
  type = "l",
  lwd = 2,

  xlab = "rate",
  ylab = "Profile log-likelihood",

  main = "Profile log-likelihood for rate"
)

abline(
  v = rate_hat,
  lty = 2
)

points(
  rate_hat,
  max(
    profile_rate,
    na.rm = TRUE
  ),
  pch = 19
)


## ---------------------------------------------------------
## 21. Final diagnostic summary
## ---------------------------------------------------------

cat("\n============================================\n",
"FINAL DIAGNOSTIC SUMMARY\n",
"============================================\n",
  "theta_hat:",
  theta_hat,
  "\n",

  "rate_hat:",
  rate_hat,
  "\n",
  "Minimum negative log-likelihood:",
  fit_best$objective,
  "\n",
  "Convergence:",
  fit_best$convergence,
  "\n",
  "Message:",
  fit_best$message,
  "\n"
)

if (theta_at_lower || theta_at_upper) {

  cat(
    "\nWARNING: theta estimate is at/near a boundary.\n"
  )

}

if (rate_at_lower || rate_at_upper) {

  cat(
    "\nWARNING: rate estimate is at/near a boundary.\n"
  )

}

if (theta_sd > 1e-3 || rate_sd > 1e-3) {

  cat(
    "\nWARNING: estimates vary across optimization starts.\n"
  )

} else {

  cat(
    "\nMultiple starts give stable estimates.\n"
  )

}

cat(
  "\nThe likelihood surface and profile likelihoods should\n",
  "be inspected for flat regions, ridges, or boundary maxima.\n",
  "============================================\n")
