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
  gamma_bias <- sapply(matrice, function(x) x[variable_name, "gamma"])
  A_bias <- sapply(matrice, function(x) x[variable_name, "scale"])
  alpha_bias <- sapply(matrice, function(x) x[variable_name, "shape"])

  # Save the plot as an image with a unique name
  png(file_name, width = 800, height = 600, res = 150)

  # Create the plot
  ylim_range <- range(c(gamma_bias, A_bias, alpha_bias), na.rm = TRUE)
  plot(gamma_vals, gamma_bias, type = "o", col = "plum", xlab = "gamma", ylab = variable_name, ylim = ylim_range)
  lines(gamma_vals, A_bias, type = "o", col = "cyan4")
  lines(gamma_vals, alpha_bias, type = "o", col = "lightgreen")

  # Add a legend (uncomment if needed)
  # legend("topright", legend = c("Power", "1/Scale", "shape"), col = c("plum", "cyan4", "lightgreen"), lty = 1, pch = 1, cex = 0.7)

  dev.off()  # Close the graphics device

  return(file_name)
}

Plot_auto_var = function(matrice = sum_matrices, file_name = "plot_variance.png"){

  # Extract empirical and theoretical variances for each parameter
  Emp_gamma_var <- sapply(matrice, function(x) x["AVG_Emp_var", "gamma"])
  Emp_A_var <- sapply(matrice, function(x) x["AVG_Emp_var", "scale"])
  Emp_alpha_var <- sapply(matrice, function(x) x["AVG_Emp_var", "shape"])

  Theo_gamma_var <- sapply(matrice, function(x) x["AVG_Theo_var", "gamma"])
  Theo_A_var <- sapply(matrice, function(x) x["AVG_Theo_var", "scale"])
  Theo_alpha_var <- sapply(matrice, function(x) x["AVG_Theo_var", "shape"])

  # Save the plot as an image with a unique name
  png(file_name, width = 800, height = 600, res = 150)

  # Automatically determine the y-axis limits
  ylim_range <- range(c(Emp_gamma_var, Emp_A_var, Emp_alpha_var,
                        Theo_gamma_var, Theo_A_var, Theo_alpha_var), na.rm = TRUE)
  #X-axis
  gamma_vals <- as.numeric(names(matrice))

  # Plot Empirical vs Theoretical Variance for gamma
  plot(gamma_vals, Emp_gamma_var, type = "o", col = "plum", xlab = "gamma", ylab = "Std",
       main = "Empirical vs Theoretical Std", ylim = ylim_range)
  lines(gamma_vals, Theo_gamma_var, type = "o", col = "plum", lty = 2)

  # Add Empirical and Theoretical Variance for scale
  lines(gamma_vals, Emp_A_var, type = "o", col = "cyan4")
  lines(gamma_vals, Theo_A_var, type = "o", col = "cyan4", lty = 2)

  # Add Empirical and Theoretical Variance for shape
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
simulation=5000

## Save Results in workbook
save = TRUE
save_plot = FALSE
save_details = TRUE

## frechet parameters
scale_obs = 1
scale_obs_inv = 1/scale_obs
shape_obs = 2  ## shape

## Series length
#T_val=50
T_values <- c(50,100 , 200)

# Define the sequence of gamma values
trend_values <- seq(1.1, 1.3, by = 0.05)

## significance level
sign = 0.05

## Likelihood expression
info = "all"

## file name if saving
file_name = paste0("Simulation_ynm_Xt_Rn_frechet_shape=",shape_obs,"_scale=", scale_obs_inv, "_simulation2_", simulation, ".xlsx")
out_path = "data/param_est/ynm/frechet_inv_scale/"

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
    colnames(params) <- c("gamma", "scale", "shape", "N_T")

    LogL <- matrix(0, nrow = simulation, ncol = 1)

    Emp_var <- matrix(0, nrow = simulation, ncol = 3)
    colnames(Emp_var) <- c("gamma", "scale", "shape")

    Theo_var <- matrix(0, nrow = simulation, ncol = 3)
    colnames(Theo_var) <- c("gamma", "scale", "shape")

    sum_matrix =  matrix(0,nrow=5,ncol=4)
    rownames(sum_matrix)= c("AVG_param", "AVG_bias","AVG_Emp_var", "AVG_Theo_var", "Coverage_proba")
    colnames(sum_matrix) = c("gamma","scale","shape","N_T")

    conf_int <- matrix(0, nrow = simulation, ncol = 6)
    colnames(conf_int) <- c("g_L", "g_U", "A_L", "A_U", "a_L", "a_U")

    ## Reset the matrices at the start of each iteration
    params_record <- matrix(0, nrow = simulation, ncol = 4)
    colnames(params_record) <- c("gamma", "scale", "shape", "N_T")

    LogL_record <- matrix(0, nrow = simulation, ncol = 1)

    Emp_var_record <- matrix(0, nrow = simulation, ncol = 3)
    colnames(Emp_var_record) <- c("gamma", "scale", "shape")

    Theo_var_record <- matrix(0, nrow = simulation, ncol = 3)
    colnames(Theo_var_record) <- c("gamma", "scale", "shape")

    sum_matrix_record =  matrix(0,nrow=5,ncol=4)
    rownames(sum_matrix_record)= c("AVG_param", "AVG_bias","AVG_Emp_var", "AVG_Theo_var", "Coverage_proba")
    colnames(sum_matrix_record) = c("gamma","scale","shape","N_T")

    conf_int_record <- matrix(0, nrow = simulation, ncol = 6)
    colnames(conf_int_record) <- c("g_L", "g_U", "A_L", "A_U", "a_L", "a_U")
    ##========================== Simulation Code =======================

    no=1:simulation

    while(length(no)>=1){ #while (sim <= simulation){
      for(sim in no){

        ## generate ynm series of frechet underlying distr
        y= YNM_series(T = T_val,gamma = trend, dist = "frechet",
                      shape = shape_obs, scale=scale_obs)

        R = rec_values(y)
        L = rec_times(y)

        while(length(R)<=1) {
          y= YNM_series(T = T_val,gamma = trend, dist = "frechet",
                        shape = shape_obs, scale=scale_obs)
          R = rec_values(y)
          L = rec_times(y)
        }  ## only one record, ignore


        ## All data
        gamma0 = exp(coef(lm(log(y) ~ seq(1,T_val)))[2])[[1]]
        logLik_fun <- loglik_registry[["YNM"]][["all"]][["frechet_inv_scale"]]
        MLE_C = estimate_model_mle(logLik_fun, obs_type = "all", data = y,
                                   lower_bounds=c(gamma=1.01, shape = 0.01, scale=0.01, location = 0),
                                   upper_bounds = c(gamma=2, shape = 10, scale=10, location = 0),
                                   start_values = c(gamma=1.2, shape = 0.01, scale = 0.5, location = 0)
        )
        MLE_C = c(MLE_C$par, "objective" = - MLE_C$objective )

        ## Records Data
        logLik_fun_rec <- loglik_registry[["YNM"]][["records"]][["frechet_inv_scale"]]

        data_rec = list(rec_values = R, rec_times = L, time = T_val)
        data_rec = data.frame(rec_values = R, rec_times = L, time = T_val)
        MLE_C_record = estimate_model_mle(logLik_fun_rec, obs_type = "records", data = data_rec,
                                         lower_bounds=c(gamma=1.001, shape = 0.01, scale=0.01, location = 0),
                                         upper_bounds =c( 100, shape = 100, scale=100, location = 0),
                                         start_values = c(
                                           gamma = 1.01,
                                           shape = 0.01,
                                           scale =0.01,
                                           location = 0
                                         )
        )

        MLE_C_record = c(MLE_C_record$par, "objective" = - MLE_C_record$objective ) ## return shape, shape, scale, objective value

        ## store estimated parameters
        params[sim,"gamma"] <- MLE_C["gamma"]
        params[sim,"scale"] <- MLE_C["scale"]
        params[sim,"shape"] <- MLE_C["shape"]
        params[sim,"N_T"] <- rec_count(y)
        LogL[sim,1] = MLE_C["objective"]


        params_record[sim,"gamma"] <- MLE_C_record["gamma"]
        params_record[sim,"scale"] <- MLE_C_record["scale"]
        params_record[sim,"shape"] <- MLE_C_record["shape"]
        params_record[sim,"N_T"] <- rec_count(y)
        LogL_record[sim,1] = MLE_C_record["objective"]
        ##____________________ Variance estimation ________________
        ### Empirical variance
        var_fun = var_logLik_records("YNM", "all","frechet_inv_scale", "all")
        emp = var_fun(data = y, params=list(gamma = MLE_C["gamma"], scale = MLE_C["scale"], shape = MLE_C["shape"]))
        #H= emp$Hessian
        #cov2cor(solve(-H))
        Emp_var[sim,"gamma"] = emp[["var_gamma"]]
        Emp_var[sim,"scale"] = emp[["var_scale"]]
        Emp_var[sim,"shape"] = emp[["var_shape"]]

        ### Theoretical variance
        theo = var_fun(data = y, params=list(gamma = trend, scale = scale_obs_inv, shape = shape_obs ) )
        Theo_var[sim,"gamma"] = theo[["var_gamma"]]
        Theo_var[sim,"scale"] = theo[["var_scale"]]
        Theo_var[sim,"shape"] = theo[["var_shape"]]

        ## Empirical Variance _record
        var_fun = var_logLik_records("YNM", "records","frechet_inv_scale", "all")
        emp = var_fun(data = data_rec, params=list(gamma = MLE_C_record["gamma"], scale = MLE_C_record["scale"], shape = MLE_C_record["shape"]))
        #emp = compute_hessian(loglik_fun = logLik_fun_rec, data = data_rec, par_vec = MLE_C_record[1:3], eps = 1e-1 )
        Emp_var_record[sim,"gamma"] = emp[["var_gamma"]]
        Emp_var_record[sim,"scale"] = emp[["var_scale"]]
        Emp_var_record[sim,"shape"] = emp[["var_shape"]]

        ## Theortical variance _record
        theo = var_fun(data = data_rec, params=list(gamma = trend, scale = scale_obs_inv, shape = shape_obs ) )
        Theo_var_record[sim,"gamma"] = theo[["var_gamma"]]
        Theo_var_record[sim,"scale"] = theo[["var_scale"]]
        Theo_var_record[sim,"shape"] = theo[["var_shape"]]


      }

      no <- which(
        sum( params[, "gamma"] <= 1, na.rm = TRUE) > 0 |
          sum (!is.finite( params[, "gamma"])) > 0 |
          rowSums(Theo_var <= 0, na.rm = TRUE) > 0 |
          rowSums(Theo_var >= 10, na.rm = TRUE) > 0 |
          rowSums(!is.finite(Theo_var)) > 0 |
          rowSums(is.na(Theo_var)) > 0 |
          rowSums(Emp_var <= 0, na.rm = TRUE) > 0 |
          rowSums(Emp_var >= 10, na.rm = TRUE) > 0 |
          rowSums(!is.finite(Emp_var)) > 0 |
          rowSums(is.na(Emp_var)) > 0 |
          sum( params_record[, "gamma"] <= 1, na.rm = TRUE) > 0 |
          sum (!is.finite( params_record[, "gamma"])) > 0 |
          rowSums(Theo_var_record <= 0, na.rm = TRUE) > 0 |
          rowSums(Theo_var_record >= 10, na.rm = TRUE) > 0 |
          rowSums(!is.finite(Theo_var_record)) > 0 |
          rowSums(is.na(Theo_var_record)) > 0 |
          rowSums(Emp_var_record <= 0, na.rm = TRUE) > 0 |
          rowSums(Emp_var_record >= 10, na.rm = TRUE) > 0 |
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
    sum_matrix[2, ] <- round(colMeans(params) - c(trend, scale_obs_inv, shape_obs, 0), 5)  # Average bias
    sum_matrix[3, ] <- round(c(colMeans(Emp_var), 0), 5)  # Average Empirical variance
    sum_matrix[4, ] <- round(c(colMeans(Theo_var), 0), 5)  # Average asymptotic variance

    ##_record
    sum_matrix_record[1, ] <- round(colMeans(params_record), 3)  # Average of parameters
    sum_matrix_record[2, ] <- round(colMeans(params_record) - c(trend, scale_obs_inv, shape_obs, 0), 5)  # Average bias
    sum_matrix_record[3, ] <- round(c(colMeans(Emp_var_record), 0), 5)  # Average Empirical variance
    sum_matrix_record[4, ] <- round(c(colMeans(Theo_var_record), 0), 5)  # Average asymptotic variance

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
      bound_g = bounds(  params[sim,"gamma"] ,z=z,Theo_var[sim,"gamma"] )
      if (trend < bound_g[2] && trend > bound_g[1]) {
        power_g <- power_g + 1 }

      bound_g_record = bounds(  params_record[sim,"gamma"] ,z=z,Theo_var_record[sim,"gamma"] )
      if (trend < bound_g_record[2] && trend > bound_g_record[1]) {
        power_g_record <- power_g_record + 1 }

      ## CI A
      bound_A = bounds(  params[sim,"scale"] ,z=z,Theo_var[sim,"scale"] )
      if (scale_obs_inv < bound_A[2] && scale_obs_inv > bound_A[1]) {
        power_A <- power_A + 1 }

      bound_A_record = bounds(  params_record[sim,"scale"] ,z=z,Theo_var_record[sim,"scale"] )
      if (scale_obs_inv < bound_A_record[2] && scale_obs_inv > bound_A_record[1]) {
        power_A_record <- power_A_record + 1 }


      ##CI a
      bound_a = bounds(  params[sim,"shape"] ,z=z,Theo_var[sim,"shape"] )
      if (shape_obs < bound_a[2] && shape_obs > bound_a[1]) {
        power_a <- power_a + 1 }

      bound_a_record = bounds(  params_record[sim,"shape"] ,z=z,Theo_var_record[sim,"shape"] )
      if (shape_obs < bound_a_record[2] && shape_obs > bound_a_record[1]) {
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
        mutate(T_value = T_val, gamma = trend) %>% # 4. Add the T_val column at the beginning (or end) of the dataframe
        relocate(T_value) # Moves T_value to the first column

      # 5. Save this iteration's data frame into the list using T_val as the key
      detailed_results[[paste0("T_", T_val, "gamma_", trend)]] <- combined_data


      df_theo_record   <- as.data.frame(Theo_var_record) %>% rename_with(~paste0("Theo_", .))
      df_emp_record    <- as.data.frame(Emp_var_record)  %>% rename_with(~paste0("Emp_", .))
      df_params_record <- as.data.frame(params_record)   %>% rename_with(~paste0("Param_", .))
      df_logl_record   <- as.data.frame(LogL_record)     %>% rename_with(~"LogL")

      # 2. Combine them all horizontally
      combined_data_record <- cbind(df_theo_record, df_emp_record, df_params_record, df_logl_record) %>%
        mutate(T_value = T_val, gamma = trend) %>% # 4. Add the T_val column at the beginning (or end) of the dataframe
        relocate(T_value) # Moves T_value to the first column

      # 5. Save this iteration's data frame into the list using T_val as the key
      detailed_results_record[[paste0("T_", T_val, "gamma_", trend)]] <- combined_data_record

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
# Check the results
#print(results_T)


## All
if(save){
  # Initialize empty list
  all_data <- list()

  for (T_name in names(results_T)) {
    for (gamma_name in names(results_T[[T_name]])) {

      mat <- results_T[[T_name]][[gamma_name]]

      df <- as.data.frame(mat)

      # Add row names as a column
      df$Metric <- rownames(df)

      # Add identifiers
      df$T <- T_name
      df$gamma_obs <- gamma_name
      df$gamma_obs = as.numeric(df$gamma_obs)

      all_data[[length(all_data) + 1]] <- df
    }
  }

  # Combine everything
  final_df <- do.call(rbind, all_data)

  # Reorder columns (optional)
  final_df <- final_df[, c("T", "gamma_obs", "Metric", colnames(mat))]
  #final_df$gamma = as.numeric(final_df$gamma)

  # Write to CSV
  addWorksheet(wb, c("All par"))
  writeData(wb,  c("All par") , final_df, rowNames = FALSE , colNames = TRUE)

}

## Record
if(save){
  # Initialize empty list
  all_data <- list()

  for (T_name in names(results_T_record)) {
    for (gamma_name in names(results_T_record[[T_name]])) {

      mat <- results_T_record[[T_name]][[gamma_name]]

      df <- as.data.frame(mat)

      # Add row names as a column
      df$Metric <- rownames(df)

      # Add identifiers
      df$T <- T_name
      df$gamma_obs <- gamma_name
      df$gamma_obs = as.numeric(df$gamma_obs)

      all_data[[length(all_data) + 1]] <- df
    }
  }

  # Combine everything
  final_df_record <- do.call(rbind, all_data)

  # Reorder columns (optional)
  final_df_record <- final_df_record[, c("T", "gamma_obs", "Metric", colnames(mat))]
  #final_df$gamma = as.numeric(final_df$gamma)

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

# ##1. plot gamma vs true gamma, grouped by T_________________
#   # Convert types if needed
# final_df$T <- factor(final_df$T)
#   # Keep only AVG_param row (estimates)
# df_param <- final_df %>%
#   filter(Metric == "AVG_param")
# ggplot(df_param, aes(x = gamma, y = gamma, color = T, group = T)) +
#   geom_line() +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   labs(
#     #title = "Estimated vs True gamma",
#     x = expression(gamma),
#     y = expression(hat(gamma)),
#   ) +
#   theme_minimal()
# ggplot(df_param, aes(x = gamma, y = gamma,
#                        linetype = T, shape = T, group = T)) +
#     geom_line(color = "black", size = 0.6) +
#     geom_point(color = "black", size = 2) +
#     geom_abline(slope = 1, intercept = 0,
#                 linetype = "dashed", color = "gray40") +
#     scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
#     labs(
#       x = expression(gamma),
#       y = expression(hat(gamma)),
#       title = NULL
#     ) +
#     theme_extremes()
#
#   ##2.  Bias plot_______
# df_bias <- final_df %>%
#     filter(Metric == "AVG_bias")
# ggplot2::ggplot(df_bias, aes(x = gamma, y = gamma, color = T)) +
#     geom_line() +
#     geom_point() +
#     labs(
#       title = "Bias of gamma Estimate",
#       x = "True gamma",
#       y = "Bias"
#     ) +
#     theme_minimal()
# ggplot(df_bias, aes(x = gamma, y = gamma,
#                     linetype = T, group = T)) +
#   geom_line(color = "black", size = 0.6) +
#   geom_point(color = "black", size = 2) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
#   scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
#   labs(
#     x = expression(gamma),
#     y = "Bias",
#     title = NULL
#   ) +
#   theme_extremes()
#
# # 3. Standard deviation comparison_____________
# df_sd <- final_df %>%
#   filter(Metric %in% c("AVG_Emp_var", "AVG_Theo_var"))
# ggplot(df_sd, aes(x = gamma, y = gamma, color = Metric)) +
#   geom_line() +
#   facet_wrap(~ T) +
#   labs(
#     title = "Empirical vs Theoretical Std Dev",
#     x = "gamma",
#     y = "Std Dev"
#   ) +
#   theme_minimal()
# ggplot(df_sd, aes(x = gamma, y = gamma,
#                   linetype = Metric, group = Metric)) +
#   geom_line(color = "black", size = 0.6) +
#   facet_wrap(~ T) +
#   scale_linetype_manual(values = c("solid", "dashed")) +
#   labs(
#     x = expression(gamma),
#     y = "Standard deviation"
#   ) +
#   theme_extremes()
#
# ## Multi -parameter plot
# df_long <- df_param %>%
#   tidyr::pivot_longer(cols = c(gamma, scale, shape),
#                names_to = "Parameter",
#                values_to = "Value")
#
# ggplot(df_long, aes(x = gamma, y = Value, color = T)) +
#   geom_line() +
#   facet_wrap(~ Parameter, scales = "free_y") +
#   labs(
#     title = "Parameter Estimates Across gamma",
#     x = "True gamma",
#     y = "Estimate"
#   ) +
#   theme_minimal()

plots <- plot_extremes_color(
  final_df,
  obs_col = "gamma_obs",
  param_cols = c("gamma", "scale", "shape")
)

# Display
plots$bias$gamma
plots$est_vs_true$gamma
plots$std_dev$gamma
plots$multi_param
plots$coverage

if(save_plot){
  # addWorksheet(wb, "gamma Estimates")
  # insertImage(wb, "gamma Estimates", plots_gamma$est_vs_true, startRow = 2, startCol = 2, width = 6, height = 4)
  #
  # addWorksheet(wb, "gamma Bias")
  # insertImage(wb, "gamma Bias", plots_gamma$bias, startRow = 2, startCol = 2, width = 6, height = 4)
  #
  # addWorksheet(wb, "gamma Std Dev")
  # insertImage(wb, "gamma Std Dev", plots_gamma$std_dev, startRow = 2, startCol = 2, width = 6, height = 4)
  #
  # addWorksheet(wb, "Multi-Param Plot")
  # insertImage(wb, "Multi-Param Plot", plots_gamma$multi_param, startRow = 2, startCol = 2, width = 6, height = 4)

  ggsave(paste0(out_path,"gamma_est_vs_true.png"), plots$est_vs_true$gamma, width = 6, height = 4)
  ggsave(paste0(out_path,"gamma_bias.png"), plots$bias$gamma, width = 6, height = 4)
  ggsave(paste0(out_path,"gamma_var.png"), plots$std_dev$gamma, width = 6, height = 4)

  ggsave(paste0(out_path,"shape_est_vs_true.png"), plots$est_vs_true$shape, width = 6, height = 4)
  ggsave(paste0(out_path,"shape_bias.png"), plots$bias$shape, width = 6, height = 4)
  ggsave(paste0(out_path,"shape_var.png"), plots$std_dev$shape, width = 6, height = 4)

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
