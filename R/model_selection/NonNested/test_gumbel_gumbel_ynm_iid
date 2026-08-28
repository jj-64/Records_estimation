#_______________________________________________________________
# PART A: TYPE-I ERROR OF THE WILKS LIKELIHOOD RATIO TEST
#
# H0 : Gumbel YNM
# H1 : Gumbel IID
#
# Test statistic:
#
#     LR = 2 ( l_H1 - l_H0 )
#
# Under Wilks theorem and assuming H0 is nested in H1:
#
#     LR ~ Chi-square(df = p1 - p0)
#
#_______________________________________________________________

set.seed(12345)
library(dplyr)
library(purrr)
library(ggplot2)
library(scales)  # for percentage formatting if needed
library(stats)
#devtools::load_all("~/Records")

#_______________________________________________________________
# DATA GENERATING MECHANISMS --------
#_______________________________________________________________

## Generate a series under H0
series_H0 <- function(T_val, par){

  # Exponential Linear Drift Model
  YNM_series(
    T     = T_val,
    dist  = "gumbel",
    gamma = par[["trend"]],
    location  = par[["location"]],
    scale = par[["scale"]]
  )
}

## Generate a series under H1
series_H1 <- function(T_val, par){

  VGAM::rgumbel(n=T_val, loc = par[["location"]], scale = par[["scale"]])
}

#_______________________________________________________________
# LIKELIHOOD FUNCTIONS --------
#_______________________________________________________________

Likelihood_under_H0 <- function(data_rec, params){

  logLik_fun_rec <-
    loglik_registry[["YNM"]][["records"]][["gumbel"]]

  logLik_fun_rec(
    data   = data_rec,
    params = params
  )
}

Likelihood_under_H1 <- function(data_rec, params){

  logLik_fun_rec <-
    loglik_registry[["iid"]][["records"]][["gumbel"]]

  logLik_fun_rec(
    data   = data_rec,
    params = params
  )
}

#_______________________________________________________________
# Part A: SIMULATION SETTINGS ------
#_______________________________________________________________
save  = TRUE

T_val      <- 75
alpha      <- 0.05
simulation <- 100

trend_values <- seq(1.1, 1.1, by = 0.05)
#Records::rec_count_mean_YNM(T=100,gamma=0.1,dist="exp",rate=2)
true_params <- c(
  trend = 1.2,
  location  = 1,
  scale = 1
)

H1_params = c(
  shape = 1.2
)

final_summary <- matrix(0, nrow = simulation, ncol = 16)
colnames(final_summary) = c("gamma_obs",
                            "gamma_H0","location_H0", "scale_H0","logLik_H0",
                            "gamma_H1", "location_H1", "scale_H1", "logLik_H1",
                            "AIC_0",  "AIC_1",
                            "BIC_0", "BIC_1",
                            "Lower", "Upper","Power")

#________________________________________________________________
# LOOP OVER TRUE DRIFT VALUES
#________________________________________________________________

for(trial in 1:simulation){

  trend_val = true_params[["trend"]]

  cat("\n====================================================\n",
      "True drift =", trend_val, "\n",
      "====================================================\n")

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
        T_val      = T_val,
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
    # STEP 2: FIT H0 (EXPONENTIAL YNM)
    #__________________________________________________________________________

    lb_0 <- c(gamma = 1.01,
              location  = -10,
              scale = 0.01)

    ub_0 <- c(gamma = 5,
              location = 10,
              scale = 10)

    x0_0 <- c(gamma = 1.1,
              location  = 0,
              scale = sd(R))


    negLogLik_H0 <- function(params){
      -Likelihood_under_H0(data_rec, params)
    }

    fit_H0 <- NULL

    for(iterations in 1:100){

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

      x0_0 <- x0_0 + c(0.01, 0.01, 0.01)
    }

    if(is.null(fit_H0) || fit_H0$convergence != 0){
      next
    }

    param_H0_start <- c(
      fit_H0$par,
      logLik = -fit_H0$objective
    )

    AIC_0 = - 2 * param_H0_start["logLik"] + 2 * length(fit_H0$par)
    BIC_0 = - 2 * param_H0_start["logLik"] + log(m) * length(fit_H0$par)
    #__________________________________________________________________________
    # STEP 3: FIT H1 (WEIBULL YNM)
    #__________________________________________________________________________

    lb_1 <- c(
      gamma = 1,
      location = -10,
      scale = 0.1
    )

    ub_1 <- c(
      gamma = 1,
      location = 10,
      scale = 10
    )

    x0_1 <- c(
      gamma = 1,#param_H0[["gamma"]],
      location = 0,
      scale = sd(R/L)
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

    param_H1_start <- c(
      fit_H1$par,
      logLik = -fit_H1$objective
    )

    AIC_1 = - 2 * param_H1_start["logLik"] + 2 * length(fit_H1$par)
    BIC_1 = - 2 * param_H1_start["logLik"] + log(m) * length(fit_H1$par)


##################### Go to the code of hypothesis testing and put the parameters there


    final_summary[trial,] <- c(
      trend_val,
      param_H0_start,

      param_H1_start,

      AIC_0,
      AIC_1,
      BIC_0,
      BIC_1,

      rep(0,3)

    )

#_____________________________________________
# Part A: The Critical Region of Lambda
#_____________________________________________

results <- data.frame(
  gamma_obs = rep(NA, simulation),
  gamma_H0 = NA,
  location_H0 = NA,
  scale_H0 = NA,
  logLik_H0 = NA,

  gamma_H1 = NA,
  location_H1 = NA,
  scale_H1 = NA,
  logLik_H1 = NA,

  NbRecords = NA,

  AIC_H0 = NA,
  AIC_H1 = NA,
  BIC_H0 = NA,
  BIC_H1 = NA

)

for (d in 1) {

  ## create matrixes to store records
  L_x <- matrix(NA, nrow = T_val, ncol = simulation)
  R_x <- matrix(NA, nrow = T_val, ncol = simulation)
  nbr_records <- numeric(simulation)

      #___________________________________________________________________________
      # A1- Generate simulated records from series under H0 (using the parameters estimated above)
      #___________________________________________________________________________

  sim=1
  while(sim <= simulation) {

    x_CR <- series_H0(T_val, c(trend = param_H0_start[["gamma"]], param_H0_start[c("location","scale")])) ##under H0

    ## Records and indicators: consider cases of more than 2 records
    R = rec_values(x_CR)
    L = rec_times(x_CR)
    m <- length(L)

    if(m >= 2){
      L_x[1:m, sim] <- L
      R_x[1:m, sim] <- R

      nbr_records[sim] <- m

      sim=sim+1}
  }

      #___________________________________________________
      # A2 - Perform H0/H1 likelihood maximization
      #___________________________________________________
  for (sim in 1:simulation) {

    lb_0 <- c(gamma = 1.01,
              location  = -10,
              scale = 0.01)

    ub_0 <- c(gamma = 5,
              location = 10,
              scale = 10)

    x0_0 <- c(gamma = 1.1,
              location  = 0,
              scale = sd(R/L))

    data_rec <- list(
      rec_values = R_x[1:nbr_records[sim], sim],
      rec_times  = L_x[1:nbr_records[sim], sim],
      time       = T_val
    )

    m = nbr_records[sim]

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

      x0_0 <- x0_0 + c(0.01, 0.10, 0.01)
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

  #______________________________________________
  # A3 - Perform H1 likelihood maximization
  #___________________________________________________

    lb_1 <- c(gamma = 1.0,
              location  = -10,
              scale = 0.01)

    ub_1 <- c(gamma = 1.0,
              location = 10,
              scale = 10)

    x0_1 <- c(gamma = 1.0,
              location  = 0,
              scale = sd(R/L))

    data_rec <- list(
      rec_values = R_x[1:nbr_records[sim], sim],
      rec_times  = L_x[1:nbr_records[sim], sim],
      time       = T_val
    )

    m = nbr_records[sim]

    negLogLik_H1 <- function(params){
      -Likelihood_under_H1(data_rec, params)
    }

    fit_H1 <- NULL

    for(iterations in 1:50){

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
         fit_H0$convergence == 0){

        break
      }

      x0_1 <- x0_1 + c(0.0, 0.10, 0.01)
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


    #_____________________________________________
    # Store Results
    #_____________________________________________

    results[sim + (simulation * (d - 1)), ] <- c(
      trend_val,
      param_H0,

      param_H1,

      m,

      AIC_0,
      AIC_1,
      BIC_0,
      BIC_1
    )
  }

}

      #_____________________________________________
      # A4 - Wilk's Lamda and CI
      #_____________________________________________
### Ignore simulations
# ignore_sim = c(ignore_sim_H1,ignore_sim_H0)
# if(length(ignore_sim)>0){
#   param_H1=param_H1[-ignore_sim,]
#   param_H0=param_H0[-ignore_sim,]
# }

for(d in 1:length(trend_values)){

  temp = subset(results, gamma_obs == trend_val)
  Rapport <-  temp$logLik_H1 -  temp$logLik_H0  ## Log L1 - Log L0

  Rapport_sorted <- sort(Rapport[is.finite(Rapport)])

  # 4 - confidence bounds

  final_summary[trial,"Lower"] = quantile(Rapport_sorted, alpha/2) ##LOwer
  final_summary[trial,"Upper"] = quantile(Rapport_sorted , 1-alpha/2) ##Upper
}

#_____________________________________________
# Part B:  Perform  power calculations
#_____________________________________________

results_power <- data.frame(
  gamma_obs = rep(NA, simulation),

  gamma_H0 = NA,
  location_H0 = NA,
  scale_H0 = NA,
  logLik_H0 = NA,

  gamma_H1 = NA,
  location_H1 = NA,
  scale_H1 = NA,
  logLik_H1 = NA,

  NbRecords = NA,

  AIC_H0 = NA,
  AIC_H1 = NA,
  BIC_H0 = NA,
  BIC_H1 = NA

)

for(d in 1) {

  ## create matrixes to store records
  L_x_test <- matrix(NA, nrow = T_val, ncol = simulation)
  R_x_test <- matrix(NA, nrow = T_val, ncol = simulation)
  nbr_records_test <- numeric(simulation)

  #___________________________________________________________________________
  # A1- Generate simulated records from series under H1
  #___________________________________________________________________________

  sim=1
  while (sim <= simulation) {

    x_Pow <- series_H1(T_val = T_val, param_H1_start[c("location", "scale")]  )

    ## Records and indicators: consider cases of more than 2 records
    R = rec_values(x_Pow)
    L = rec_times(x_Pow)
    m <- length(L)

    if(m >= 2){
      L_x_test[1:m, sim] <- L
      R_x_test[1:m, sim] <- R

      nbr_records_test[sim] <- m

      sim=sim+1}
  }

  #___________________________________________________
  # B2 - Perform H0/H1 likelihood maximization
  #___________________________________________________
  for (sim in 1:simulation) {

    lb_0 <- c(gamma = 1.01,
              location  = -10,
              scale = 0.01)

    ub_0 <- c(gamma = 5,
              location = 10,
              scale = 10)

    x0_0 <- c(gamma = 1.1,
              location  = 0,
              scale = sd(R/L))

    data_rec <- list(
      rec_values = R_x_test[1:nbr_records_test[sim], sim],
      rec_times  = L_x_test[1:nbr_records_test[sim], sim],
      time       = T_val
    )

    m = nbr_records[sim]

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

      x0_0 <- x0_0 + c(0.01, 0.10, 0.01)
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

    #______________________________________________
    # 3 - Perform H1 likelihood maximization
    #___________________________________________________

    lb_1 <- c(gamma = 1.0,
              location  = -10,
              scale = 0.01)

    ub_1 <- c(gamma = 1.0,
              location = 10,
              scale = 10)

    x0_1 <- c(gamma = 1.0,
              location  = 0,
              scale = sd(R/L))

    data_rec <- list(
      rec_values = R_x_test[1:nbr_records_test[sim], sim],
      rec_times  = L_x_test[1:nbr_records_test[sim], sim],
      time       = T_val
    )

    m = nbr_records_test[sim]

    negLogLik_H1 <- function(params){
      -Likelihood_under_H1(data_rec, params)
    }

    fit_H1 <- NULL

    for(iterations in 1:50){

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
         fit_H0$convergence == 0){

        break
      }

      x0_1 <- x0_1 + c(0.0, 0.10, 0.01)
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


    #_____________________________________________
    # Store Results
    #_____________________________________________

    results_power[sim + (simulation * (d - 1)), ] <- c(
      trend_val,
      param_H0,

      param_H1,

      m,

      AIC_0,
      AIC_1,
      BIC_0,
      BIC_1
    )
  }

}

#_____________________________________________
# Part C:  Calculate likelihood ratios if valid estimates are obtained
#_____________________________________________
# #ignore_sim_test = c(which(is.na(param_H0_test)),which(is.na(param_H1_test)))
# if(length(ignore_sim_test)>0){
#   param_H0_test = param_H0_test[-ignore_sim_test,]
#   param_H1_test = param_H1_test[-ignore_sim_test,]
# }

#_____________________________________________
# Part C:  Compute Rapport Lambda
#_____________________________________________

for(d in 1){

  temp = subset(results_power, gamma_obs == trend_val)
  Rapport_test <- temp$logLik_H1 - temp$logLik_H0  ## Log L1 - Log L0

  # Test power calculation based on confidence interval bounds

  final_summary[trial,"Power"] = sum(Rapport_test > final_summary[trial,"Upper"] | Rapport_test < final_summary[trial,"Lower"])*100/simulation #length(maxlogL_H0_test)
}

}
