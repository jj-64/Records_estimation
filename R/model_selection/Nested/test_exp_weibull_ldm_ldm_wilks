###############################################################################
# PART A: TYPE-I ERROR OF THE WILKS LIKELIHOOD RATIO TEST
#
# H0 : Exponential LDM
# H1 : Weibull LDM
#
# Test statistic:
#
#     LR = 2 ( l_H1 - l_H0 )
#
# Under Wilks theorem and assuming H0 is nested in H1:
#
#     LR ~ Chi-square(df = p1 - p0)
#
###############################################################################

set.seed(12345)
library(dplyr)
library(purrr)
library(ggplot2)
library(scales)  # for percentage formatting if needed

###############################################################################
# DATA GENERATING MECHANISMS --------
###############################################################################

## Generate a series under H0
series_H0 <- function(T_val, par){

  # Exponential Linear Drift Model
  LDM_series(
    T     = T_val,
    dist  = "exp",
    theta = par[["trend"]],
    rate  =par[["rate"]]
  )
}

## Generate a series under H1
series_H1 <- function(T_val, par){

  # Weibull Linear Drift Model
  LDM_series(
    T      = T_val,
    dist   = "weibull",
    theta  = par[["trend"]],
    shape  = par[["shape"]],
    scale  =par[["scale"]]
  )
}

###############################################################################
# LIKELIHOOD FUNCTIONS --------
###############################################################################

Likelihood_under_H0 <- function(data_rec, params){

  logLik_fun_rec <-
    loglik_registry[["LDM"]][["records"]][["exp"]]

  logLik_fun_rec(
    data   = data_rec,
    params = params
  )
}

Likelihood_under_H1 <- function(data_rec, params){

  logLik_fun_rec <-
    loglik_registry[["LDM"]][["records"]][["weibull"]]

  logLik_fun_rec(
    data   = data_rec,
    params = params
  )
}

###############################################################################
# Part A: SIMULATION SETTINGS ------
###############################################################################
save  = TRUE

T_val      <- 75
alpha      <- 0.05
simulation <- 1000

trend_values <- seq(0.1,0.2, by = 0.1)
#Records::rec_count_mean_LDM(T=100,theta=0.1,dist="exp",rate=2)
true_params <- c(
  trend = NA,
  rate  = 2
)

H1_params = c(
  shape = 1.2
)

final_summary <- matrix(0, nrow = length(trend_values), ncol = 4)
colnames(final_summary) = c("trend", "Type_I_Error", "Power", "Nb_record")

kk <- 1

results_list = list()
results_power_list = list()

###############################################################################
# LOOP OVER TRUE DRIFT VALUES
###############################################################################

for(trend_val in trend_values){

  cat("\n====================================================\n")
  cat("True drift =", trend_val, "\n")
  cat("====================================================\n")

  true_params["trend"] = trend_val

  results <- data.frame(
    theta_H0 = rep(NA, simulation),
    rate_H0  = NA,
    logLik_H0 = NA,

    theta_H1 = NA,
    shape_H1 = NA,
    scale_H1 = NA,
    logLik_H1 = NA,

    LR       = NA,
    p_value  = NA,
    reject_H0  = NA,
    nRecords = NA,

    AIC_H0 = NA,
    AIC_H1 = NA,
    BIC_H0 = NA,
    BIC_H1 = NA

  )

  trial <- 1

  ###########################################################################
  # MONTE-CARLO LOOP
  ###########################################################################

  while(trial <= simulation){

    print(paste("Under H0, Trend = ", trend_val, ", Trial = ", trial))
    #__________________________________________________________________________
    # STEP 1: GENERATE DATA UNDER H0
    #__________________________________________________________________________

    xt <- series_H0(
      T_val      = T_val,
      par = true_params
    )

    R <- rec_values(xt)
    L <- rec_times(xt)

    ## Reject degenerate samples with only one record
    while(length(R) <= 1){

      xt <- series_H0(
        T      = T_val,
        par = true_params
      )

      R <- rec_values(xt)
      L <- rec_times(xt)
    }

    m <- length(R)

    data_rec <- list(
      rec_values = R,
      rec_times  = L,
      time       = T_val
    )

    #__________________________________________________________________________
    # STEP 2: FIT H0 (EXPONENTIAL LDM)
    #__________________________________________________________________________

    lb_0 <- c(theta = 0.01,
              rate  = 0.1)

    ub_0 <- c(theta = min(R / L),
              rate  = 10)

    x0_0 <- c(theta = 0.5 * min(R/L),# 0.01,
              rate  = 1/mean(diff(c(0,R))))


    negLogLik_H0 <- function(params){
      -Likelihood_under_H0(data_rec, params)
    }

    fit_H0 <- NULL

    for(iterations in 1:50){

      fit_H0 <- tryCatch(
        nlminb(
          start = x0_0,
          objective = negLogLik_H0,
          lower = lb_0,
          upper = ub_0
        ),
        error = function(e) NULL
      )

      if(!is.null(fit_H0) &&
         is.finite(fit_H0$objective) &&
         fit_H0$convergence == 0){

        break
      }

      x0_0 <- x0_0 + c(0.01, 0.10)
    }

    if(is.null(fit_H0) || fit_H0$convergence != 0){
      next
    }

    param_H0 <- c(
      fit_H0$par,
      logLik = -fit_H0$objective
    )

    AIC_0 = - 2 * param_H0["logLik"] + 2 * length(fit_H0$par)
    BIC_0 = - 2 * param_H0["logLik"] + log(m) * length(fit_H0$par)
    #__________________________________________________________________________
    # STEP 3: FIT H1 (WEIBULL LDM)
    #__________________________________________________________________________

    lb_1 <- c(
      theta = 0.01,
      shape = 0.1,
      scale = 0.1
    )

    ub_1 <- c(
      theta = min(R/L),
      shape = 5,
      scale = 10
    )

    x0_1 <- c(
      theta = 0.1,#param_H0[["theta"]],
      shape = 1,
      scale = 1#param_H0[["rate"]]
    )

    negLogLik_H1 <- function(params){
      -Likelihood_under_H1(data_rec, params)
    }

    fit_H1 <- NULL

    for(iterations in 1:70){

      fit_H1 <- tryCatch(
        nlminb(
          start = x0_1,
          objective = negLogLik_H1,
          lower = lb_1,
          upper = ub_1
        ),
        error = function(e) NULL
      )

      if(!is.null(fit_H1) &&
         is.finite(fit_H1$objective) &&
         fit_H1$convergence == 0){

        break
      }

      x0_1 <- x0_1 + c(0.01, 0.01, 0)
    }

    if(is.null(fit_H1) || fit_H1$convergence != 0){
      next
    }

    param_H1 <- c(
      fit_H1$par,
      logLik = -fit_H1$objective
    )

    AIC_1 = - 2 * param_H1["logLik"] + 2 * length(fit_H1$par)
    BIC_1 = - 2 * param_H1["logLik"] + log(m) * length(fit_H1$par)
    #__________________________________________________________________________
    # STEP 4: WILKS LIKELIHOOD RATIO TEST
    #__________________________________________________________________________

    logLik_H0 <- param_H0["logLik"]
    logLik_H1 <- param_H1["logLik"]

    LR <- 2 * (logLik_H1 - logLik_H0)

    ## Numerical safeguard
    LR <- max(0, LR)
    #if(LR <0) next;

    df_wilks <- 1   # Weibull adds one free parameter (shape)

    p_value <- 1 - pchisq(
      q  = LR,
      df = df_wilks
    )

    reject = ifelse(p_value < alpha , 1, 0)
    #__________________________________________________________________________
    # STEP 5: STORE RESULTS
    #__________________________________________________________________________

    results[trial, ] <- c(
      param_H0["theta"],
      param_H0["rate"],
      logLik_H0,

      param_H1["theta"],
      param_H1["shape"],
      param_H1["scale"],
      logLik_H1,

      LR,
      p_value,
      reject,
      m,

      AIC_0,
      AIC_1,
      BIC_0,
      BIC_1
    )

    trial <- trial + 1
  }

  #__________________________________________________________________________
  # ESTIMATED TYPE-I ERROR
  #__________________________________________________________________________

  # rejection_rate <- mean(
  #   results$p_value < alpha,
  #   na.rm = TRUE
  # )

  #results$rejection_H0 = ifelse( results$p_value < alpha,1,0)

  #cat("Estimated Type-I Error =", mean(results$rejection_H0), "\n")

###############################################################################
# PART B: EMPIRICAL POWER OF THE WILKS TEST
#
# Data generated under H1 (Weibull-LDM)
#
# H0 : Exponential LDM
# H1 : Weibull LDM
#
###############################################################################

results_power <- data.frame(
  theta_H0  = rep(NA, simulation),
  rate_H0   = NA,
  logLik_H0 = NA,

  theta_H1  = NA,
  shape_H1  = NA,
  scale_H1  = NA,
  logLik_H1 = NA,

  LR        = NA,
  p_value   = NA,
  reject    = NA,
  nRecords  = NA,

  AIC_H0 = NA,
  AIC_H1 = NA,
  BIC_H0 = NA,
  BIC_H1 = NA
)

trial <- 1

###############################################################################
# MONTE-CARLO LOOP
###############################################################################

while(trial <= simulation){

  print(paste("Under H1, Trend = ", trend_val, ", Trial = ", trial))

  #_________________________________________________________________________
  # STEP 1: GENERATE DATA UNDER H1
  #_________________________________________________________________________

  xt_H1 <- series_H1(
    T_val      = T_val,
    par = c(trend = true_params[["trend"]],
      shape = 2,
      scale = 1/true_params[["rate"]]
    )
  )

  R <- rec_values(xt_H1)
  L <- rec_times(xt_H1)

  while(length(R) <= 1){

    xt_H1 <- series_H1(
      T_val      = T_val,
      par = c(trend = true_params[["trend"]],
              shape = H1_params[["shape"]],
              scale = 1/true_params[["rate"]]
      )
    )

    R <- rec_values(xt_H1)
    L <- rec_times(xt_H1)
  }

   m = length(R)

  data_rec <- list(
    rec_values = R,
    rec_times  = L,
    time       = T_val
  )

  #_________________________________________________________________________
  # STEP 2: FIT H0 (EXPONENTIAL LDM)
  #_________________________________________________________________________

  lb_0 <- c(
    theta = 0.01,
    rate  = 0.01
  )

  ub_0 <- c(
    theta = min(R / L),
    rate  = 10
  )

  x0_0 <- c(
    theta = 0.5 * min(R/L), # 0.10,
    rate  = mean(diff(c(0,R)))# 0.10
  )

  negLogLik_H0 <- function(params){
    -Likelihood_under_H0(data_rec, params)
  }

  fit_H0 <- NULL

  for(iterations in 1:50){

    fit_H0 <- tryCatch(
      nlminb(
        start     = x0_0,
        objective = negLogLik_H0,
        lower     = lb_0,
        upper     = ub_0
      ),
      error = function(e) NULL
    )

    if(!is.null(fit_H0) &&
       is.finite(fit_H0$objective) &&
       fit_H0$convergence == 0){
      break
    }

    x0_0 <- x0_0 + c(0.10, 0.10)
  }

  if(is.null(fit_H0) || fit_H0$convergence != 0){
    next
  }

  param_H0 <- c(
    fit_H0$par,
    logLik = -fit_H0$objective
  )

  AIC_0 = - 2 * param_H0["logLik"] + 2 * length(fit_H0$par)
  BIC_0 = - 2 * param_H0["logLik"] + log(m) * length(fit_H0$par)
  #_________________________________________________________________________
  # STEP 3: FIT H1 (WEIBULL LDM)
  #_________________________________________________________________________

  lb_1 <- c(
    theta = 0.01,
    shape = 0.01,
    scale = 0.01
  )

  ub_1 <- c(
    theta = min(R / L),
    shape = 5,
    scale = 5
  )

  x0_1 <- c(
    theta =0.5 * min(R/L), # 0.10, param_H0[["theta"]],
    shape = 1.20,
    scale = mean(diff(c(0,R)))#param_H0[["rate"]]
  )

  negLogLik_H1 <- function(params){
    -Likelihood_under_H1(data_rec, params)
  }

  fit_H1 <- NULL

  for(iterations in 1:70){

    fit_H1 <- tryCatch(
      nlminb(
        start     = x0_1,
        objective = negLogLik_H1,
        lower     = lb_1,
        upper     = ub_1
      ),
      error = function(e) NULL
    )

    if(!is.null(fit_H1) &&
       is.finite(fit_H1$objective) &&
       fit_H1$convergence == 0){
      break
    }

    x0_1 <- x0_1 + c(0.010, 0.10, 0.10)
  }

  if(is.null(fit_H1) || fit_H1$convergence != 0){
    next
  }

  param_H1 <- c(
    fit_H1$par,
    logLik = -fit_H1$objective
  )

  AIC_1 = - 2 * param_H1["logLik"] + 2 * length(fit_H1$par)
  BIC_1 = - 2 * param_H1["logLik"] + log(m) * length(fit_H1$par)
  #_________________________________________________________________________
  # STEP 4: WILKS LIKELIHOOD RATIO TEST
  #_________________________________________________________________________

  logLik_H0 <- param_H0["logLik"]
  logLik_H1 <- param_H1["logLik"]

  LR <- 2 * (logLik_H1 - logLik_H0)

  ## Numerical safeguard
  #LR <- max(0, LR)
  if(LR < 0) next;

  critical_value <- qchisq(
    p  = 1 - alpha,
    df = 1
  )

  p_value <- 1 - pchisq(
    q  = LR,
    df = df_wilks
  )

  reject <- as.integer(
    LR >= critical_value
  )

  #_________________________________________________________________________
  # STEP 5: STORE RESULTS
  #_________________________________________________________________________

  results_power[trial, ] <- c(
    param_H0["theta"],
    param_H0["rate"],
    logLik_H0,

    param_H1["theta"],
    param_H1["shape"],
    param_H1["scale"],
    logLik_H1,

    LR,
    p_value,
    reject,
    m,

    AIC_0,
    AIC_1,
    BIC_0,
    BIC_1
  )

  trial <- trial + 1
}

###############################################################################
# TYPE I ERROR ---------
###############################################################################

type1_error <- mean(
  results$reject_H0,
  na.rm = TRUE
)

###############################################################################
# EMPIRICAL POWER
###############################################################################

empirical_power <- mean(
  results_power$reject,
  na.rm = TRUE
)

###############################################################################
# REPORT
###############################################################################

cat(
  "\nType I Error =",
  round(100 * type1_error, 2),
  "%\n"
)

cat(
  "Power =",
  round(100 * empirical_power, 2),
  "%\n"
)

###############################################################################
# SAVE FINAL RESULTS
###############################################################################

final_summary[kk, ] <- c(
  trend_val,
  100 * type1_error,
  100 * empirical_power,
  mean(results$nRecords, na.rm = TRUE)
)

results_list[[kk]] = results
results_power_list[[kk]] = results_power

kk <- kk + 1
}

##############################################################################
# Save Results
##############################################################################

# 1. Name the list items and combine into one big dataframe
final_result <- results_list %>%
  set_names(trend_values) %>%
  bind_rows(.id = "trend_value") # Creates the column at the start

# 2. (Optional) Convert the trend column to numeric if needed
final_result$trend_value <- as.numeric(final_result$trend_value)

# 1. Name the list items and combine into one big dataframe
final_result_power <- results_power_list %>%
  set_names(trend_values) %>%
  bind_rows(.id = "trend_value") # Creates the column at the start

# 2. (Optional) Convert the trend column to numeric if needed
final_result_power$trend_value <- as.numeric(final_result_power$trend_value)

library(openxlsx)

if(save){
# 1. Create a named list of your dataframes
sheets_list <- list(
  "type_I_error" = final_result,    # First sheet name = First dataframe
  "Power"   = final_result_power,     # Second sheet name = Second dataframe
  "summary"      = final_summary,          # Third sheet name = Third dataframe
  "metadata" = data.frame(T_val = T_val, rate = true_params["rate"], n_sim = simulation,
                          df = df_wilks, nb_param_0 = length(fit_H0$par),
                          nb_param_1 = length(fit_H1$par) )
)

# 2. Save everything into one workbook
write.xlsx(sheets_list, file = paste0("data/model_selection/nested/model_selection_exp_",true_params[["rate"]],"_weibull_",H1_params[["shape"]],"_ldm_ldm_T=",T_val,".xlsx"))

paste0("file written to:" , "data/model_selection/nested/model_selection_exp_",true_params[["rate"]],"_weibull_",H1_params[["shape"]],"_ldm_ldm_T=",T_val,".xlsx")
}

###############################################################################
# Plot
###############################################################################
# Ensure your data is a dataframe
final_summary <- as.data.frame(final_summary)

# Rescale Nb_record to match the scale of the left y-axis
scale_factor <- max(final_summary$Power) / max(final_summary$Nb_record)

ggplot(final_summary, aes(x = trend)) +
  geom_line(aes(y = Type_I_Error, color = "Type I Error"), size = 1) +
  geom_line(aes(y = Power, color = "Power"), size = 1) +
  geom_line(aes(y = Nb_record * scale_factor, color = "Nb_record"), linetype = "dashed", size = 1) +

  # Left y-axis
  scale_y_continuous(
    name = "Test Power | Type I Error (%)",

    # Right y-axis
    sec.axis = sec_axis(~ . / scale_factor, name = "Number of Records")
  ) +

  scale_color_manual(values = c("Type I Error" = "lightblue", "Power" = "darkblue", "Nb_record" = "#2E8B57")) +

  labs(x = "Theta", color = "Metric") +
  theme_classic() +
  theme(
    axis.title.y.right = element_text(color = "#2E8B57"),
    axis.title.y.left = element_text(color = "black"),
    legend.position = "top"
  )

