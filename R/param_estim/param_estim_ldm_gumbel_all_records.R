library(stats)
library(dplyr)
library(GA)
library(pracma)
library(maxLik)
#library("xlsx")
library(VGAM)
library(openxlsx)
options(scipen = 6)

Plot_auto = function(variable_name, matrice = sum_matrices, file_name = "plot_bias.png"){

  # Extract values from sum_matrices
  gamma_vals <- as.numeric(names(matrice))
  gamma_bias <- sapply(matrice, function(x) x[variable_name, "theta"])
  A_bias <- sapply(matrice, function(x) x[variable_name, "scale"])
  alpha_bias <- sapply(matrice, function(x) x[variable_name, "location"])

  # Save the plot as an image with a unique name
  png(file_name, width = 800, height = 600, res = 150)

  # Create the plot
  ylim_range <- range(c(gamma_bias, A_bias, alpha_bias), na.rm = TRUE)
  plot(gamma_vals, gamma_bias, type = "o", col = "plum", xlab = "Theta", ylab = variable_name, ylim = ylim_range)
  lines(gamma_vals, A_bias, type = "o", col = "cyan4")
  lines(gamma_vals, alpha_bias, type = "o", col = "lightgreen")

  # Add a legend (uncomment if needed)
  # legend("topright", legend = c("Power", "1/Scale", "location"), col = c("plum", "cyan4", "lightgreen"), lty = 1, pch = 1, cex = 0.7)

  dev.off()  # Close the graphics device

  return(file_name)
}

Plot_auto_var = function(matrice = sum_matrices, file_name = "plot_variance.png"){

  # Extract empirical and theoretical variances for each parameter
  Emp_gamma_var <- sapply(matrice, function(x) x["AVG_Emp_var", "theta"])
  Emp_A_var <- sapply(matrice, function(x) x["AVG_Emp_var", "scale"])
  Emp_alpha_var <- sapply(matrice, function(x) x["AVG_Emp_var", "location"])

  Theo_gamma_var <- sapply(matrice, function(x) x["AVG_Theo_var", "theta"])
  Theo_A_var <- sapply(matrice, function(x) x["AVG_Theo_var", "scale"])
  Theo_alpha_var <- sapply(matrice, function(x) x["AVG_Theo_var", "location"])

  # Save the plot as an image with a unique name
  png(file_name, width = 800, height = 600, res = 150)

  # Automatically determine the y-axis limits
  ylim_range <- range(c(Emp_gamma_var, Emp_A_var, Emp_alpha_var,
                        Theo_gamma_var, Theo_A_var, Theo_alpha_var), na.rm = TRUE)
  #X-axis
  gamma_vals <- as.numeric(names(matrice))

  # Plot Empirical vs Theoretical Variance for theta
  plot(gamma_vals, Emp_gamma_var, type = "o", col = "plum", xlab = "Theta", ylab = "Std",
       main = "Empirical vs Theoretical Std", ylim = ylim_range)
  lines(gamma_vals, Theo_gamma_var, type = "o", col = "plum", lty = 2)

  # Add Empirical and Theoretical Variance for scale
  lines(gamma_vals, Emp_A_var, type = "o", col = "cyan4")
  lines(gamma_vals, Theo_A_var, type = "o", col = "cyan4", lty = 2)

  # Add Empirical and Theoretical Variance for location
  lines(gamma_vals, Emp_alpha_var, type = "o", col = "lightgreen")
  lines(gamma_vals, Theo_alpha_var, type = "o", col = "lightgreen", lty = 2)

  dev.off()  # Close the graphics device

  return(file_name)
}

list_to_df= function(sum_matrices, var="T"){
  df <- data.frame(matrix(unlist(sum_matrices), nrow=length(sum_matrices), byrow=TRUE))
  rownames(df) = paste0(var,"=",names(sum_matrices))
  colnames(df) =
    with(expand.grid(rownames(sum_matrices[[1]]), colnames(sum_matrices[[1]])), paste0(Var1,"_",Var2))
  return(df)
}

###########################   Parameters and Dataframes ###########################
simulation=10000

## Save Results in workbook
save = TRUE
save_plot = TRUE
save_details = TRUE

## gumbel parameters
scale_obs = 1
location_obs=0  ## location

## Series length
#T_val=50
T_values <- c(50, 100, 200)

# Define the sequence of theta values
trend_values <- seq(0.1, 0.50, by = 0.1)

## significance level
sign = 0.05

## Likelihood expression
info = "all"

## file name if saving
file_name = paste0("Simulation_ldm_Xt_Rn_gumbel_location=",location_obs,"_scale=", scale_obs, ".xlsx")
out_path = "data/param_est/ldm/gumbel/"

###########################  Simulation  ###########################

# Create a new workbook
wb <- createWorkbook()

results_T <- list()
results_T_record <- list()

# Create a list to save detailed results
detailed_results <- list()
detailed_results_record <- list()

for (T_val in T_values) {
  # Create a list to store sum_matrix for each trend
  sum_matrices <- list()
  sum_matrices_record <- list()

  # Loop over each trend value
  for (trend in trend_values) {

    ## Reset the matrices at the start of each iteration
    params <- matrix(0, nrow = simulation, ncol = 4)
    colnames(params) <- c("theta", "scale", "location", "N_T")

    LogL <- matrix(0, nrow = simulation, ncol = 1)

    Emp_var <- matrix(0, nrow = simulation, ncol = 3)
    colnames(Emp_var) <- c("theta", "scale", "location")

    Theo_var <- matrix(0, nrow = simulation, ncol = 3)
    colnames(Theo_var) <- c("theta", "scale", "location")

    sum_matrix =  matrix(0,nrow=5,ncol=4)
    rownames(sum_matrix)= c("AVG_param", "AVG_bias","AVG_Emp_var", "AVG_Theo_var", "Coverage_proba")
    colnames(sum_matrix) = c("theta","scale","location","N_T")

    conf_int <- matrix(0, nrow = simulation, ncol = 6)
    colnames(conf_int) <- c("g_L", "g_U", "A_L", "A_U", "a_L", "a_U")

    ## Reset the matrices at the start of each iteration
    params_record <- matrix(0, nrow = simulation, ncol = 4)
    colnames(params_record) <- c("theta", "scale", "location", "N_T")

    LogL_record <- matrix(0, nrow = simulation, ncol = 1)

    Emp_var_record <- matrix(0, nrow = simulation, ncol = 3)
    colnames(Emp_var_record) <- c("theta", "scale", "location")

    Theo_var_record <- matrix(0, nrow = simulation, ncol = 3)
    colnames(Theo_var_record) <- c("theta", "scale", "location")

    sum_matrix_record =  matrix(0,nrow=5,ncol=4)
    rownames(sum_matrix_record)= c("AVG_param", "AVG_bias","AVG_Emp_var", "AVG_Theo_var", "Coverage_proba")
    colnames(sum_matrix_record) = c("theta","scale","location","N_T")

    conf_int_record <- matrix(0, nrow = simulation, ncol = 6)
    colnames(conf_int_record) <- c("g_L", "g_U", "A_L", "A_U", "a_L", "a_U")
    ##========================== Simulation Code =======================

    no=1:simulation

    while(length(no)>=1){ #while (sim <= simulation){
      for(sim in no){

        ## generate LDM series of gumbel underlying distr
        y= LDM_series(T = T_val,theta = trend, dist = "gumbel",
                      location = location_obs, scale=scale_obs)

        R = rec_values(y)
        L = rec_times(y)

        while(length(R)<=1) {
          y= LDM_series(T = T_val,theta = trend, dist = "gumbel",
                        location = location_obs, scale=scale_obs)
          R = rec_values(y)
          L = rec_times(y)
        }  ## only one record, ignore


        ## All data
          logLik_fun <- loglik_registry[["LDM"]][["all"]][["gumbel"]]
          MLE_C = estimate_model_mle(logLik_fun, obs_type = "all", data = y,
                                     lower_bounds=c(theta=0.01, location = 0.01, scale=0.01),
                                     upper_bounds = c(theta=10, location = 100, scale=100),
                                     start_values = c(theta=0.01, location = 0.01, scale=0.01)
          )
          MLE_C = c(MLE_C$par, "objective" = - MLE_C$objective )

          ## Records Data
          logLik_fun_rec <- loglik_registry[["LDM"]][["records"]][["gumbel"]]

          data_rec = list(rec_values = R, rec_times = L, time = T_val)
          data_rec = data.frame(rec_values = R, rec_times = L, time = T_val)
          MLE_C_record = estimate_model_mle(logLik_fun_rec, obs_type = "records", data = data_rec,
                                     lower_bounds=c(theta = 0.01, location = 0.01, scale=0.01),
                                     upper_bounds = c(100, location = 100, scale=100),
                                     start_values = c(
                                       theta = 0.01,
                                       location = 0.01,
                                       scale =0.01
                                     )
          )
          MLE_C_record = c(MLE_C_record$par, "objective" = - MLE_C_record$objective ) ## return location, location, scale, objective value

        ## store estimated parameters
        params[sim,"theta"] <- MLE_C["theta"]
        params[sim,"scale"] <- MLE_C["scale"]
        params[sim,"location"] <- MLE_C["location"]
        params[sim,"N_T"] <- rec_count(y)
        LogL[sim,1] = MLE_C["objective"]


        params_record[sim,"theta"] <- MLE_C_record["theta"]
        params_record[sim,"scale"] <- MLE_C_record["scale"]
        params_record[sim,"location"] <- MLE_C_record["location"]
        params_record[sim,"N_T"] <- rec_count(y)
        LogL_record[sim,1] = MLE_C_record["objective"]
        ##____________________ Variance estimation ________________
        ### Empirical variance
          var_fun = var_logLik_records("LDM", "all","gumbel", "all")
          emp = var_fun(data = y, params=list(theta = MLE_C["theta"], scale = MLE_C["scale"], location = MLE_C["location"]))
          Emp_var[sim,"theta"] = emp[["var_theta"]]
          Emp_var[sim,"scale"] = emp[["var_scale"]]
          Emp_var[sim,"location"] = emp[["var_location"]]

          ### Theoretical variance
          theo = var_fun(data = y, params=list(theta = trend, scale = scale_obs, location = location_obs ) )
          Theo_var[sim,"theta"] = theo[["var_theta"]]
          Theo_var[sim,"scale"] = theo[["var_scale"]]
          Theo_var[sim,"location"] = theo[["var_location"]]

        ## Empirical Variance _record
          var_fun = var_logLik_records("LDM", "records","gumbel", "all")
          emp = var_fun(data = data_rec, params=list(theta = MLE_C_record["theta"], scale = MLE_C_record["scale"], location = MLE_C_record["location"]))
          #compute_var_hessian(loglik_fun = logLik_fun_rec, data_rec, par_vec = MLE_C, eps = 1e-1 )
          Emp_var_record[sim,"theta"] = emp[["var_theta"]]
          Emp_var_record[sim,"scale"] = emp[["var_scale"]]
          Emp_var_record[sim,"location"] = emp[["var_location"]]

      ## Theortical variance _record
          theo = var_fun(data = data_rec, params=list(theta = trend, scale = scale_obs, location = location_obs ) )
          Theo_var_record[sim,"theta"] = theo[["var_theta"]]
          Theo_var_record[sim,"scale"] = theo[["var_scale"]]
          Theo_var_record[sim,"location"] = theo[["var_location"]]


    }

      no <- which(
        sum( params[, "theta"] <= 0, na.rm = TRUE) > 0 |
          sum (!is.finite( params[, "theta"])) > 0 |
          rowSums(Theo_var <= 0, na.rm = TRUE) > 0 |
          rowSums(!is.finite(Theo_var)) > 0 |
          rowSums(is.na(Theo_var)) > 0 |
          rowSums(Emp_var <= 0, na.rm = TRUE) > 0 |
          rowSums(!is.finite(Emp_var)) > 0 |
          rowSums(is.na(Emp_var)) > 0 |
          sum( params_record[, "theta"] <= 0, na.rm = TRUE) > 0 |
          sum (!is.finite( params_record[, "theta"])) > 0 |
          rowSums(Theo_var_record <= 0, na.rm = TRUE) > 0 |
          rowSums(!is.finite(Theo_var_record)) > 0 |
          rowSums(is.na(Theo_var_record)) > 0 |
          rowSums(Emp_var_record <= 0, na.rm = TRUE) > 0 |
          rowSums(!is.finite(Emp_var_record)) > 0 |
          rowSums(is.na(Emp_var_record)) > 0
      )
      if(length(no)>=1){
        cat("Re-simulating for ", length(no), " simulations due to invalid estimates...\n")
      }
    }

    ##___________________________________________________________

    ## Compute the summary matrix
    sum_matrix[1, ] <- round(colMeans(params), 3)  # Average of parameters
    sum_matrix[2, ] <- round(colMeans(params) - c(trend, scale_obs, location_obs, 0), 5)  # Average bias
    sum_matrix[3, ] <- round(c((colMeans(Emp_var)), 0), 5)  # Average Empirical variance
    sum_matrix[4, ] <- round(c((colMeans(Theo_var)), 0), 5)  # Average asymptotic variance

    ##_record
    sum_matrix_record[1, ] <- round(colMeans(params_record), 3)  # Average of parameters
    sum_matrix_record[2, ] <- round(colMeans(params_record) - c(trend, scale_obs, location_obs, 0), 5)  # Average bias
    sum_matrix_record[3, ] <- round(c((colMeans(Emp_var_record)), 0), 5)  # Average Empirical variance
    sum_matrix_record[4, ] <- round(c((colMeans(Theo_var_record)), 0), 5)  # Average asymptotic variance

    ## Compute coverage probability
    z <- qnorm(1 - (sign / 2), 0, 1)

    power_g <- 0
    power_A <- 0
    power_a <- 0

    power_g_record <- 0
    power_A_record <- 0
    power_a_record <- 0

    for (sim in 1:simulation) {
      ## CI trend
      bound_g = bounds(  params[sim,"theta"] ,z=z,Theo_var[sim,"theta"] )
      if (trend < bound_g[2] && trend > bound_g[1]) {
        power_g <- power_g + 1 }

      bound_g_record = bounds(  params_record[sim,"theta"] ,z=z,Theo_var_record[sim,"theta"] )
      if (trend < bound_g_record[2] && trend > bound_g_record[1]) {
        power_g_record <- power_g_record + 1 }

      ## CI A
      bound_A = bounds(  params[sim,"scale"] ,z=z,Theo_var[sim,"scale"] )
      if (scale_obs < bound_A[2] && scale_obs > bound_A[1]) {
        power_A <- power_A + 1 }

      bound_A_record = bounds(  params_record[sim,"scale"] ,z=z,Theo_var_record[sim,"scale"] )
      if (scale_obs < bound_A_record[2] && scale_obs > bound_A_record[1]) {
        power_A_record <- power_A_record + 1 }


      ##CI a
      bound_a = bounds(  params[sim,"location"] ,z=z,Theo_var[sim,"location"] )
      if (location_obs < bound_a[2] && location_obs > bound_a[1]) {
        power_a <- power_a + 1 }

      bound_a_record = bounds(  params_record[sim,"location"] ,z=z,Theo_var_record[sim,"location"] )
      if (location_obs < bound_a_record[2] && location_obs > bound_a_record[1]) {
        power_a_record <- power_a_record + 1 }

      ## fill into matrix
      conf_int[sim,]= c(bound_g,bound_A,bound_a)
      conf_int_record[sim,]= c(bound_g_record,bound_A_record,bound_a_record)
    }

    CP_trend <- power_g / simulation
    CP_A <- power_A / simulation
    CP_alpha <- power_a / simulation
    sum_matrix[5, ] <- c(CP_trend, CP_A, CP_alpha, 0)

    CP_trend_record <- power_g_record / simulation
    CP_A_record <- power_A_record / simulation
    CP_alpha_record <- power_a_record / simulation
    sum_matrix_record[5, ] <- c(CP_trend_record, CP_A_record, CP_alpha_record, 0)

    ## Store the sum_matrix for this trend
    sum_matrices[[as.character(trend)]] <- sum_matrix
    sum_matrices_record[[as.character(trend)]] <- sum_matrix_record

    ## Saving
    if(save_details){
      # 1. Convert to data frames and add distinct prefixes to column names
      df_theo   <- as.data.frame(Theo_var) %>% rename_with(~paste0("Theo_", .))
      df_emp    <- as.data.frame(Emp_var)  %>% rename_with(~paste0("Emp_", .))
      df_params <- as.data.frame(params)   %>% rename_with(~paste0("Param_", .))
      df_logl   <- as.data.frame(LogL)     %>% rename_with(~"LogL")

      # 2. Combine them all horizontally
      combined_data <- cbind(df_theo, df_emp, df_params, df_logl) %>%
        mutate(T_value = T_val, theta = trend) %>% # 4. Add the T_val column at the beginning (or end) of the dataframe
        relocate(T_value) # Moves T_value to the first column

      # 5. Save this iteration's data frame into the list using T_val as the key
      detailed_results[[paste0("T_", T_val, "theta_", trend)]] <- combined_data


      df_theo_record   <- as.data.frame(Theo_var_record) %>% rename_with(~paste0("Theo_", .))
      df_emp_record    <- as.data.frame(Emp_var_record)  %>% rename_with(~paste0("Emp_", .))
      df_params_record <- as.data.frame(params_record)   %>% rename_with(~paste0("Param_", .))
      df_logl_record   <- as.data.frame(LogL_record)     %>% rename_with(~"LogL")

      # 2. Combine them all horizontally
      combined_data_record <- cbind(df_theo_record, df_emp_record, df_params_record, df_logl_record) %>%
        mutate(T_value = T_val, theta = trend) %>% # 4. Add the T_val column at the beginning (or end) of the dataframe
        relocate(T_value) # Moves T_value to the first column

      # 5. Save this iteration's data frame into the list using T_val as the key
      detailed_results_record[[paste0("T_", T_val, "theta_", trend)]] <- combined_data_record

    }

  }

  results_T[[paste0("T=", T_val)]] <- sum_matrices
  results_T_record[[paste0("T=", T_val)]] <- sum_matrices_record

}

if(save_details){
  # All
  final_combined_data <- bind_rows(detailed_results)
  addWorksheet(wb, "All_Results")
  writeData(wb, sheet = "All_Results", final_combined_data, rowNames = FALSE, colNames = TRUE)

  # Record
  final_combined_data_record <- bind_rows(detailed_results_record)
  addWorksheet(wb, "All_Results_record")
  writeData(wb, sheet = "All_Results_record", final_combined_data_record, rowNames = FALSE, colNames = TRUE)

  ## summary table by Number of records
  summary_table_NT <- final_combined_data %>%
    group_by(Param_N_T) %>%
    dplyr::summarise(
      count_sim = n(),
      Proba_N_T = n() / nrow(final_combined_data),
      across(everything(), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::relocate(count_sim, Proba_N_T, .after = last_col())
  addWorksheet(wb, "summary_NT")
  writeData(wb, sheet = "summary_NT", summary_table_NT, rowNames = FALSE, colNames = TRUE)

  summary_table_NT_record <- final_combined_data_record %>%
    dplyr::group_by(Param_N_T) %>%
    dplyr::summarise(
      count_sim = n(),
      Proba_N_T = n() / nrow(final_combined_data_record),
      across(everything(), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::relocate(count_sim, Proba_N_T, .after = last_col())
  addWorksheet(wb, "summary_NT_record")
  writeData(wb, sheet = "summary_NT_record", summary_table_NT_record, rowNames = FALSE, colNames = TRUE)


}

## All
if(save){
  # Initialize empty list
  all_data <- list()

  for (T_name in names(results_T)) {
    for (theta_name in names(results_T[[T_name]])) {

      mat <- results_T[[T_name]][[theta_name]]

      df <- as.data.frame(mat)

      # Add row names as a column
      df$Metric <- rownames(df)

      # Add identifiers
      df$T <- T_name
      df$theta_obs <- theta_name
      df$theta_obs = as.numeric(df$theta_obs)

      all_data[[length(all_data) + 1]] <- df
    }
  }

  # Combine everything
  final_df <- do.call(rbind, all_data)

  # Reorder columns (optional)
  final_df <- final_df[, c("T", "theta_obs", "Metric", colnames(mat))]
  #final_df$theta = as.numeric(final_df$theta)

  # Write to CSV
  addWorksheet(wb, c("All par"))
  writeData(wb,  c("All par") , final_df, rowNames = FALSE , colNames = TRUE)

}

## Record
if(save){
  # Initialize empty list
  all_data <- list()

  for (T_name in names(results_T_record)) {
    for (theta_name in names(results_T_record[[T_name]])) {

      mat <- results_T_record[[T_name]][[theta_name]]

      df <- as.data.frame(mat)

      # Add row names as a column
      df$Metric <- rownames(df)

      # Add identifiers
      df$T <- T_name
      df$theta_obs <- theta_name
      df$theta_obs = as.numeric(df$theta_obs)

      all_data[[length(all_data) + 1]] <- df
    }
  }

  # Combine everything
  final_df_record <- do.call(rbind, all_data)

  # Reorder columns (optional)
  final_df_record <- final_df_record[, c("T", "theta_obs", "Metric", colnames(mat))]
  #final_df$theta = as.numeric(final_df$theta)

  # Write to CSV
  addWorksheet(wb, c("All par_record"))
  writeData(wb,  c("All par_record") , final_df_record, rowNames = FALSE , colNames = TRUE)

}

# Plots -------------------------------------------------------------------

# Bias plot vs trend
p = Plot_auto(variable_name = "AVG_bias", file_name = "plot_bias.png")
addWorksheet(wb, "Bias vs. trend")
insertImage(wb, "Bias vs. trend", p, startRow = 2, startCol = 2, width = 6, height = 4)

# Coverage Probability plot vs trend
p2 = Plot_auto(variable_name = "Coverage_proba", file_name = "plot_coverage.png")
addWorksheet(wb, "Cov vs. trend")
insertImage(wb, "Cov vs. trend", p2, startRow = 2, startCol = 2, width = 6, height = 4)

# Variance plot
p3 = Plot_auto_var(sum_matrices, file_name = "plot_variance.png")
addWorksheet(wb, "Var vs. trend")
insertImage(wb, "Var vs. trend", p3, startRow = 2, startCol = 2, width = 6, height = 4)

# ##1. plot theta vs true theta, grouped by T_________________
#   # Convert types if needed
# final_df$T <- factor(final_df$T)
#   # Keep only AVG_param row (estimates)
# df_param <- final_df %>%
#   filter(Metric == "AVG_param")
# ggplot(df_param, aes(x = theta, y = theta, color = T, group = T)) +
#   geom_line() +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   labs(
#     #title = "Estimated vs True Theta",
#     x = expression(theta),
#     y = expression(hat(theta)),
#   ) +
#   theme_minimal()
# ggplot(df_param, aes(x = theta, y = theta,
#                        linetype = T, location = T, group = T)) +
#     geom_line(color = "black", size = 0.6) +
#     geom_point(color = "black", size = 2) +
#     geom_abline(slope = 1, intercept = 0,
#                 linetype = "dashed", color = "gray40") +
#     scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
#     labs(
#       x = expression(theta),
#       y = expression(hat(theta)),
#       title = NULL
#     ) +
#     theme_extremes()
#
#   ##2.  Bias plot_______
# df_bias <- final_df %>%
#     filter(Metric == "AVG_bias")
# ggplot2::ggplot(df_bias, aes(x = theta, y = theta, color = T)) +
#     geom_line() +
#     geom_point() +
#     labs(
#       title = "Bias of Theta Estimate",
#       x = "True Theta",
#       y = "Bias"
#     ) +
#     theme_minimal()
# ggplot(df_bias, aes(x = theta, y = theta,
#                     linetype = T, group = T)) +
#   geom_line(color = "black", size = 0.6) +
#   geom_point(color = "black", size = 2) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
#   scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
#   labs(
#     x = expression(theta),
#     y = "Bias",
#     title = NULL
#   ) +
#   theme_extremes()
#
# # 3. Standard deviation comparison_____________
# df_sd <- final_df %>%
#   filter(Metric %in% c("AVG_Emp_var", "AVG_Theo_var"))
# ggplot(df_sd, aes(x = theta, y = theta, color = Metric)) +
#   geom_line() +
#   facet_wrap(~ T) +
#   labs(
#     title = "Empirical vs Theoretical Std Dev",
#     x = "Theta",
#     y = "Std Dev"
#   ) +
#   theme_minimal()
# ggplot(df_sd, aes(x = theta, y = theta,
#                   linetype = Metric, group = Metric)) +
#   geom_line(color = "black", size = 0.6) +
#   facet_wrap(~ T) +
#   scale_linetype_manual(values = c("solid", "dashed")) +
#   labs(
#     x = expression(theta),
#     y = "Standard deviation"
#   ) +
#   theme_extremes()
#
# ## Multi -parameter plot
# df_long <- df_param %>%
#   tidyr::pivot_longer(cols = c(theta, scale, location),
#                names_to = "Parameter",
#                values_to = "Value")
#
# ggplot(df_long, aes(x = theta, y = Value, color = T)) +
#   geom_line() +
#   facet_wrap(~ Parameter, scales = "free_y") +
#   labs(
#     title = "Parameter Estimates Across Theta",
#     x = "True Theta",
#     y = "Estimate"
#   ) +
#   theme_minimal()

plots <- plot_extremes_color(
  final_df,
  obs_col = "theta_obs",
  param_cols = c("theta", "scale", "location")
)

# Display
plots$bias$theta
plots$est_vs_true$theta
plots$std_dev$theta
plots$multi_param
plots$coverage

if(save_plot){
  # addWorksheet(wb, "Theta Estimates")
  # insertImage(wb, "Theta Estimates", plots_theta$est_vs_true, startRow = 2, startCol = 2, width = 6, height = 4)
  #
  # addWorksheet(wb, "Theta Bias")
  # insertImage(wb, "Theta Bias", plots_theta$bias, startRow = 2, startCol = 2, width = 6, height = 4)
  #
  # addWorksheet(wb, "Theta Std Dev")
  # insertImage(wb, "Theta Std Dev", plots_theta$std_dev, startRow = 2, startCol = 2, width = 6, height = 4)
  #
  # addWorksheet(wb, "Multi-Param Plot")
  # insertImage(wb, "Multi-Param Plot", plots_theta$multi_param, startRow = 2, startCol = 2, width = 6, height = 4)

  ggsave(paste0(out_path,"theta_est_vs_true.png"), plots$est_vs_true$theta, width = 6, height = 4)
  ggsave(paste0(out_path,"theta_bias.png"), plots$bias$theta, width = 6, height = 4)
  ggsave(paste0(out_path,"theta_var.png"), plots$std_dev$theta, width = 6, height = 4)

  ggsave(paste0(out_path,"location_est_vs_true.png"), plots$est_vs_true$location, width = 6, height = 4)
  ggsave(paste0(out_path,"location_bias.png"), plots$bias$location, width = 6, height = 4)
  ggsave(paste0(out_path,"location_var.png"), plots$std_dev$location, width = 6, height = 4)

  ggsave(paste0(out_path,"scale_est_vs_true.png"), plots$est_vs_true$scale, width = 6, height = 4)
  ggsave(paste0(out_path,"scale_bias.png"), plots$bias$scale, width = 6, height = 4)
  ggsave(paste0(out_path,"scale_var.png"), plots$std_dev$scale, width = 6, height = 4)

  ggsave(paste0(out_path,"multi_param.png"), plots$multi_param, width = 6, height = 4)
  ggsave(paste0(out_path,"coverage_probability.png"), plots$coverage, width = 6, height = 4)

}

## Save the workbook--------------------
if(save){
  saveWorkbook(wb, paste0(out_path,file_name), overwrite = TRUE)
  print(paste0("file saved to ",out_path,file_name ))
}
