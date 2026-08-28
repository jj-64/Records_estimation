library(stats)
library(GA)
library(VGAM)
library(pracma)
library(maxLik)
#library("xlsx")
library(openxlsx)
options(scipen = 6)

Plot_auto = function(variable_name, matrice = sum_matrices, file_name = "plot_bias.png", xlab = "Series length (T)"){

  # Extract values from sum_matrices
  gamma_vals <- as.numeric(names(matrice))
  gamma_bias <- sapply(matrice, function(x) x[variable_name, "gamma"])
  A_bias <- sapply(matrice, function(x) x[variable_name, "min"])
  alpha_bias <- sapply(matrice, function(x) x[variable_name, "max"])

  # Save the plot as an image with a unique name
  png(file_name, width = 800, height = 600, res = 150)

  # Create the plot
  ylim_range <- range(c(gamma_bias, A_bias, alpha_bias), na.rm = TRUE)
  plot(gamma_vals, alpha_bias, type = "o", col = "plum", xlab = xlab , ylab = variable_name, ylim = ylim_range)
  #lines(gamma_vals, A_bias, type = "o", col = "cyan4")
  #lines(gamma_vals, alpha_bias, type = "o", col = "lightgreen")

  # Add a legend (uncomment if needed)
  # legend("topright", legend = c("Power", "1/Scale", "Shape"), col = c("plum", "cyan4", "lightgreen"), lty = 1, pch = 1, cex = 0.7)

  dev.off()  # Close the graphics device

  return(file_name)
}

Plot_auto_var = function(matrice = sum_matrices, file_name = "plot_variance.png"){

  # Extract empirical and theoretical variances for each parameter
  Emp_gamma_var <- sapply(matrice, function(x) x["AVG_Emp_var", "gamma"])
  Emp_A_var <- sapply(matrice, function(x) x["AVG_Emp_var", "min"])
  Emp_alpha_var <- sapply(matrice, function(x) x["AVG_Emp_var", "max"])

  Theo_gamma_var <- sapply(matrice, function(x) x["AVG_Theo_var", "gamma"])
  Theo_A_var <- sapply(matrice, function(x) x["AVG_Theo_var", "min"])
  Theo_alpha_var <- sapply(matrice, function(x) x["AVG_Theo_var", "max"])

  # Save the plot as an image with a unique name
  png(file_name, width = 800, height = 600, res = 150)

  # Automatically determine the y-axis limits
  ylim_range <- range(c(Emp_gamma_var, Emp_A_var, Emp_alpha_var,
                        Theo_gamma_var, Theo_A_var, Theo_alpha_var), na.rm = TRUE)
  #X-axis
  gamma_vals <- as.numeric(names(matrice))

  # Plot Empirical vs Theoretical Variance for gamma
  plot(gamma_vals, Emp_gamma_var, type = "o", col = "plum", xlab = "Gamma", ylab = "var",
       main = "Empirical vs Theoretical var", ylim = ylim_range)
  lines(gamma_vals, Theo_gamma_var, type = "o", col = "plum", lty = 2)

  # Add Empirical and Theoretical Variance for min
  lines(gamma_vals, Emp_A_var, type = "o", col = "cyan4")
  lines(gamma_vals, Theo_A_var, type = "o", col = "cyan4", lty = 2)

  # Add Empirical and Theoretical Variance for max
  lines(gamma_vals, Emp_alpha_var, type = "o", col = "lightgreen")
  lines(gamma_vals, Theo_alpha_var, type = "o", col = "lightgreen", lty = 2)

  dev.off()  # Close the graphics device

  return(file_name)
}

list_to_df = function(sum_matrices, var="T"){
  df <- data.frame(matrix(unlist(sum_matrices), nrow=length(sum_matrices), byrow=TRUE))
  rownames(df) = paste0(var,"=",names(sum_matrices))
  colnames(df) =
    with(expand.grid(rownames(sum_matrices[[1]]), colnames(sum_matrices[[1]])), paste0(Var1,"_",Var2))
  return(df)
}

###########################   Parameters ###########################
simulation=5000

## Save Results in workbook
save = TRUE
save_plot = TRUE
save_details = TRUE

## Series length
T_values=seq(50,300,by =10 )

## maximum (parameter of unfiorm)
max_obs = 1
min_obs = - max_obs

# Define the sequence of theta values
gamma =1

## significance level
sign = 0.05

## Likelihood expression
info = "all"

## biased?
biased = FALSE ## TRUE: no bias correction; FALSE: bias correction

## file name if saving
file_name = paste0("Simulation_dtrw_Xt_Rn_overT_unif_max=",max_obs,"_biasCorrected_simulation_",simulation,".xlsx")
out_path = "data/param_est/dtrw/unif/"
###########################  Simulation  ###########################
# Create a new workbook
wb <- openxlsx::createWorkbook()

# Create a list to store sum_matrix
sum_matrices <- list()
sum_matrices_record <- list()

# Create a list to save detailed results
detailed_results <- list()
detailed_results_record <- list()

# Loop over each gamma value
for (T_val in T_values) {

  ## Reset the matrices at the start of each iteration
  params <- matrix(0, nrow = simulation, ncol = 4)
  colnames(params) <- c("gamma", "min", "max", "N_T")

  LogL <- matrix(0, nrow = simulation, ncol = 1)

  Emp_var <- matrix(0, nrow = simulation, ncol = 3)
  colnames(Emp_var) <- c("gamma", "min", "max")

  Theo_var <- matrix(0, nrow = simulation, ncol = 3)
  colnames(Theo_var) <- c("gamma", "min", "max")

  sum_matrix = matrix(0,nrow=5,ncol=4)
  rownames(sum_matrix)= c("AVG_param", "AVG_bias","AVG_Emp_var", "AVG_Theo_var", "Coverage_proba")
  colnames(sum_matrix) = c("gamma","min","max","N_T")

  conf_int <- matrix(0, nrow = simulation, ncol = 6)
  colnames(conf_int) <- c("g_L", "g_U", "min_L", "min_U", "max_L", "max_U")

  ## Reset the matrices at the start of each iteration
  params_record <- matrix(0, nrow = simulation, ncol = 4)
  colnames(params_record) <- c("gamma", "min", "max", "N_T")

  LogL_record <- matrix(0, nrow = simulation, ncol = 1)

  Emp_var_record <- matrix(0, nrow = simulation, ncol = 3)
  colnames(Emp_var_record) <- c("gamma", "min", "max")

  Theo_var_record <- matrix(0, nrow = simulation, ncol = 3)
  colnames(Theo_var_record) <- c("gamma", "min", "max")

  sum_matrix_record = matrix(0,nrow=5,ncol=4)
  rownames(sum_matrix_record)= c("AVG_param", "AVG_bias","AVG_Emp_var", "AVG_Theo_var", "Coverage_proba")
  colnames(sum_matrix_record) = c("gamma","min","max","N_T")

  conf_int_record <- matrix(0, nrow = simulation, ncol = 6)
  colnames(conf_int_record) <- c("g_L", "g_U", "A_L", "A_U", "a_L", "a_U")

  ############################ Simulation Code #############################

  no=1:simulation

  while(length(no)>=1){ #while (sim <= simulation){
    for(sim in no){

      ## generate yang series of Frechet underlying distr
      y=DTRW_series(T =  T_val, dist = "unif", min = - max_obs, max = max_obs )
      R=rec_values(y)
      L = rec_times(y)

      while(length(R)<=1) {
        y=DTRW_series(T=T_val, dist = "unif", min = -max_obs, max = max_obs )
        R=rec_values(y)
        L = rec_times(y)
      }  ## only one record, ignore

      ## Likelihood optimizer
       ## All
        logLik_fun <- loglik_registry[["DTRW"]][["all"]][["unif"]]
        # MLE_C = estimate_model_mle(logLik_fun_rec, obs_type = "all", data = y,
        #                            lower_bounds=c(max = 0.01),
        #                            upper_bounds = c(max = 50),
        #                            start_values = c(max = sd((y)) ) )
        # MLE_C = c(-MLE_C$par, MLE_C$par, - MLE_C$objective ) ## return mean, max, objective value
        max_est = max(abs(diff(y)))
        if(biased == FALSE) {max_est = (T_val)/(T_val-1) * max_est}
        MLE_C = c( min = -max_est, max = max_est, objective = - logLik_fun(data = y, params = c(max = max_est)) )

        ## Records
        logLik_fun_rec <- loglik_registry[["DTRW"]][["records"]][["unif"]]
        ## Records Data
        data_rec = list(rec_values = R, rec_times = L, time = T_val)
        data_rec = data.frame(rec_values = R, rec_times = L, time = T_val)
        # MLE_C = estimate_model_mle(logLik_fun_rec, obs_type = "records", data = data_rec,
        #                            lower_bounds=c(max = 0.01),
        #                            upper_bounds = c(max = 50),
        #                            start_values = c(max = sd((y)) ) )
        # MLE_C = c(-MLE_C$par, MLE_C$par, - MLE_C$objective ) ## return mean, max, objective value
        max_est_record = max(abs(diff(R)))
        if(biased == FALSE) {max_est_record = max_est_record * (T_val)/(T_val-1)} #bias correction
        MLE_C_record = c(min = -max_est_record, max= max_est_record, objective = - logLik_fun_rec(data = data_rec, params = c(max = max_est)) )


      ## store estimated parameters
      params[sim,1] <- 1  ## Gamma always 1 in iid
      params[sim,"min"] <- MLE_C[["min"]]
      params[sim,"max"] <- MLE_C["max"]
      params[sim,"N_T"] <- rec_count(y)
      LogL[sim,1] = MLE_C[["objective"]]

      params_record[sim,1] <- 1  ## Gamma always 1 in iid
      params_record[sim, "min"] <- MLE_C_record[["min"]]
      params_record[sim, "max"] <- MLE_C_record[["max"]]
      params_record[sim,"N_T"] <- rec_count(y)
      LogL_record[sim,1] = MLE_C_record[["objective"]]

      ########################## Variance estimation ##################
      ## All data
      var_fun = var_logLik_records("dtrw", "all","unif", "max")

      emp_var = var_fun(data = y, params=list(max = MLE_C["max"], min = MLE_C["min"] ),biased = biased)
      Emp_var[sim,1] = 1
      Emp_var[sim, "min"] = emp_var[["max"]]
      Emp_var[sim, "max"] = emp_var[["max"]]

      theo_var = var_fun(data = y, params=list(max = max_obs, min = min_obs ), biased = biased)
      Theo_var[sim,1] =1
      Theo_var[sim,"min"] = theo_var
      Theo_var[sim,"max"] = theo_var

      ## Record data
      var_fun = var_logLik_records("dtrw", "records","unif", "max")

      emp_var = var_fun(data = data_rec, params=list(max = MLE_C["max"], min = MLE_C["min"] ),biased = biased)
      Emp_var_record[sim,1] = 1
      Emp_var_record[sim, "min"] = emp_var[["max"]]
      Emp_var_record[sim, "max"] = emp_var[["max"]]

      theo_var = var_fun(data = data_rec, params=list(max = max_obs, min = min_obs ), biased = biased)
      Theo_var_record[sim,1] =1
      Theo_var_record[sim,"min"] = theo_var
      Theo_var_record[sim,"max"] = theo_var
    }

    no <- which(
      rowSums(Theo_var <= 0, na.rm = TRUE) > 0 |
        rowSums(!is.finite(Theo_var)) > 0 |
        rowSums(is.na(Theo_var)) > 0 |
        rowSums(Emp_var <= 0, na.rm = TRUE) > 0 |
        rowSums(!is.finite(Emp_var)) > 0 |
        rowSums(is.na(Emp_var)) > 0 |
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

    ##########################################################################

    ## Compute the summary matrix
    sum_matrix[1, ] <- round(colMeans(params), 3)  # Average of parameters
    sum_matrix[2, ] <- round(colMeans(params) - c(1,-max_obs, max_obs, 0), 5)  # Average bias
    sum_matrix[3, ] <- round(c(colMeans(Emp_var), 0), 5)  # Average Empirical variance
    sum_matrix[4, ] <- round(c(colMeans(Theo_var), 0), 5)  # Average asymptotic variance

    ##_record
    sum_matrix_record[1, ] <- round(colMeans(params_record), 3)  # Average of parameters
    sum_matrix_record[2, ] <- round(colMeans(params_record) - c( 1, -max_obs, max_obs, 0), 5)  # Average bias
    sum_matrix_record[3, ] <- round(c(colMeans(Emp_var_record), 0), 5)  # Average Empirical variance
    sum_matrix_record[4, ] <- round(c(colMeans(Theo_var_record), 0), 5)  # Average asymptotic variance


    ## Compute coverage probability
    z <- qnorm(1 - (sign / 2), 0, 1)

    power_max <- 0
    power_max_record <- 0

    for (sim in 1:simulation) {

      ##CI a
      bound_max = bounds(  params[sim,"max"] ,z=z,Theo_var[sim,"max"] )
      if (max_obs < bound_max[2] && max_obs > bound_max[1]) {
        power_max <- power_max + 1 }

      bound_max_record = bounds(  params_record[sim,"max"] ,z=z,Theo_var_record[sim,"max"] )
      if (max_obs < bound_max_record[2] && max_obs > bound_max_record[1]) {
        power_max_record <- power_max_record + 1 }

      ## fill into matrix
      conf_int[sim,]= c(rep(0,4),bound_max)
      conf_int_record[sim,]= c(rep(0,4),bound_max_record)

    }

    CP_max <- power_max / simulation
    sum_matrix[5, ] <- c(0,0, CP_max, 0)

    CP_max_record <- power_max_record / simulation
    sum_matrix_record[5, ] <- c(0,0, CP_max_record, 0)

    ## Store the sum_matrix for this gamma
    sum_matrices[[as.character(T_val)]] <- sum_matrix
    sum_matrices_record[[as.character(T_val)]] <- sum_matrix_record

    ## Saving
    if(save_details){
      # Add each matrix to a new sheet
      # 1. Convert to data frames and add distinct prefixes to column names
      df_theo   <- as.data.frame(Theo_var) %>% rename_with(~paste0("Theo_", .))
      df_emp    <- as.data.frame(Emp_var)  %>% rename_with(~paste0("Emp_", .))
      df_params <- as.data.frame(params)   %>% rename_with(~paste0("Param_", .))
      df_logl   <- as.data.frame(LogL)     %>% rename_with(~"LogL")

      # 2. Combine them all horizontally
      combined_data <- cbind(df_theo, df_emp, df_params, df_logl)
      # Remove any column containing "gamma"
      combined_data <- combined_data %>%
        dplyr::select(-contains(c("gamma", "min"))) %>%
        mutate(T_value = T_val) %>% # 4. Add the T_val column at the beginning (or end) of the dataframe
        relocate(T_value) # Moves T_value to the first column
      # 5. Save this iteration's data frame into the list using T_val as the key
      detailed_results[[paste0("T_", T_val)]] <- combined_data

      df_theo_record   <- as.data.frame(Theo_var_record) %>% rename_with(~paste0("Theo_", .))
      df_emp_record    <- as.data.frame(Emp_var_record)  %>% rename_with(~paste0("Emp_", .))
      df_params_record <- as.data.frame(params_record)   %>% rename_with(~paste0("Param_", .))
      df_logl_record   <- as.data.frame(LogL_record)     %>% rename_with(~"LogL")

      # 2. Combine them all horizontally
      combined_data_record <- cbind(df_theo_record, df_emp_record, df_params_record, df_logl_record)
      # Remove any column containing "gamma"
      combined_data_record <- combined_data_record %>%
        dplyr::select(-contains(c("gamma","min")) )%>%
        mutate(T_value = T_val) %>% # 4. Add the T_val column at the beginning (or end) of the dataframe
        relocate(T_value) # Moves T_value to the first column
      # 5. Save this iteration's data frame into the list using T_val as the key
      detailed_results[[paste0("T_", T_val)]] <- combined_data
      detailed_results_record[[paste0("T_", T_val)]] <- combined_data_record

    }

}

## Check the results
# print(sum_matrices)

if(save_details){
  final_combined_data <- bind_rows(detailed_results)
  addWorksheet(wb, "All_Results") # 3. Save to a single sheet in Excel
  writeData(wb, sheet = "All_Results", final_combined_data, rowNames = FALSE, colNames = TRUE)

  final_combined_data_record <- bind_rows(detailed_results_record)
  addWorksheet(wb, "All_Results_record") # 3. Save to a single sheet in Excel
  writeData(wb, sheet = "All_Results_record", final_combined_data_record, rowNames = FALSE, colNames = TRUE)

  ## summary table by Number of records
  summary_table_NT <- final_combined_data %>%
    dplyr::group_by(Param_N_T) %>%
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

## Add to EXCEL
if(save) {
  # Remove any column containing "gamma"
  All_Par <- list_to_df(sum_matrices, var="T")
  All_Par <- All_Par %>%
    dplyr::select(-contains(c("gamma", "min")))
  addWorksheet(wb, c("All Par"))
  writeData(wb, c("All Par"), All_Par, rowNames = TRUE , colNames = TRUE)

  All_Par_record <- list_to_df(sum_matrices_record, var="T")
  All_Par_record <- All_Par_record %>%
    dplyr::select(-contains(c("gamma", "min")))
  addWorksheet(wb, c("All Par_record"))
  writeData(wb, c("All Par_record"), All_Par_record, rowNames = TRUE , colNames = TRUE)
}

# Plots -------------------------------------------------------------------

if(save_plot){
  # Bias plot vs Gamma
  p = Plot_auto(variable_name = "AVG_bias", file_name = "plot_bias.png", matrice = sum_matrices)
  addWorksheet(wb, "Bias vs. gamma")
  insertImage(wb, "Bias vs. gamma", p, startRow = 2, startCol = 2, width = 6, height = 4)

  # Coverage Probability plot vs Gamma
  p2 = Plot_auto(variable_name = "Coverage_proba", file_name = "plot_coverage.png")
  addWorksheet(wb, "Cov vs. gamma")
  insertImage(wb, "Cov vs. gamma", p2, startRow = 2, startCol = 2, width = 6, height = 4)

  # Variance plot
  p3 = Plot_auto_var(sum_matrices, file_name = "plot_variance.png")
  addWorksheet(wb, "Var vs. gamma")
  insertImage(wb, "Var vs. gamma", p3, startRow = 2, startCol = 2, width = 6, height = 4)
}

# Save the workbook -------------------------------------------------------
if(save){
  saveWorkbook(wb, paste0(out_path,file_name), overwrite = TRUE)
  print(paste0("file saved to ",out_path,file_name ))
}



