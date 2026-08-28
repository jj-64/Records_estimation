compute_coverage <- function(est, theo_std, true_value, sign_level = 0.05) {

  z <- qnorm(1 - sign_level / 2)

  lower <- est - z * theo_std
  upper <- est + z * theo_std

  covered <- (true_value >= lower) & (true_value <= upper)

  return(mean(covered, na.rm = TRUE))
}


library(dplyr)


## by number of records
aggregate_by_NT <- function(data,
                           true_values,
                           param_names = c("sd"),
                           sign_level = 0.05) {

  total_sim <- nrow(data)
  z <- qnorm(1 - sign_level / 2)

  results_list <- list()

  for (param in param_names) {

    # Dynamically build column names
    est_col  <- paste0("Param_", param)
    emp_col  <- paste0("Emp_", param)
    theo_col <- paste0("Theo_", param)

    # Extract true value for this parameter
    true_value <- true_values[[param]]

    # Safety check
    if (!all(c(est_col, emp_col, theo_col) %in% colnames(data))) {
      warning(paste("Missing columns for parameter:", param))
      next
    }

    # Aggregate
    res <- data %>%
      group_by(Param_N_T) %>%
      summarise(
        Parameter = param,

        Mean_Est = mean(.data[[est_col]], na.rm = TRUE),

        true_value = true_value,

        Bias = mean(.data[[est_col]] - true_value, na.rm = TRUE),

        Emp_std = mean(.data[[emp_col]]^0.5, na.rm = TRUE),
        Theo_std = mean(.data[[theo_col]]^0.5, na.rm = TRUE),

        EmpVar = mean((.data[[emp_col]]), na.rm = TRUE),
        TheoVar = mean((.data[[theo_col]]), na.rm = TRUE),

        Count = n(),
        Percent = n() / total_sim,

        # Coverage (vectorized inside summarise)
        Coverage = mean(
          (true_value >= (.data[[est_col]] - z * .data[[theo_col]]^0.5)) &
            (true_value <= (.data[[est_col]] + z * .data[[theo_col]]^0.5)),
          na.rm = TRUE
        ),

        .groups = "drop"
      )

    results_list[[param]] <- res
  }

  # Combine all parameters
  final_df <- bind_rows(results_list)

  return(final_df)
}


# agg_T <- aggregate_by_T(
#   data = sim_data,
#   param_names = c("theta", "scale", "location"),
#   true_values = list(
#     scale = 1,
#     location = 0
#   )
# )

## Case 2: only theta
# agg_theta <- aggregate_by_T(
#   data = sim_data,
#   param_names = "theta"
# )
aggregate_by_T <- function(data,
                           true_values = list(),     # named list for constant params
                           param_names = c("sd"),
                           sign_level = 0.05) {

  total_sim <- nrow(data)
  z <- qnorm(1 - sign_level / 2)

  results_list <- list()

  for (param in param_names) {

    # Column names
    est_col  <- paste0("Param_", param)
    emp_col  <- paste0("Emp_", param)
    theo_col <- paste0("Theo_", param)

    # Detect TRUE value source
    if (param %in% c("theta", "gamma") && param %in% colnames(data)) {
      true_col_mode <- TRUE
    } else {
      true_col_mode <- FALSE
      if (is.null(true_values[[param]])) {
        stop(paste("Missing true value for parameter:", param))
      }
      true_value_const <- true_values[[param]]
    }

    # Check columns exist
    required_cols <- c(est_col, emp_col, theo_col)
    if (!all(required_cols %in% colnames(data))) {
      warning(paste("Missing columns for parameter:", param))
      next
    }

    # Aggregate
    res <- data %>%
      group_by(T_value) %>%
      summarise(
        Parameter = param,

        Mean_Est = mean(.data[[est_col]], na.rm = TRUE),

        # True value
        true_value = if (true_col_mode) {
          .data[[param]][1]
        } else {
          true_value_const[1]
        },

        # Bias
        Bias = if (true_col_mode) {
          mean(.data[[est_col]] - .data[[param]], na.rm = TRUE)
        } else {
          mean(.data[[est_col]] - true_value_const, na.rm = TRUE)
        },

        Emp_std = mean(.data[[emp_col]]^0.5, na.rm = TRUE),
        Theo_std = mean(.data[[theo_col]]^0.5, na.rm = TRUE),

        EmpVar = mean((.data[[emp_col]]), na.rm = TRUE),
        TheoVar = mean((.data[[theo_col]]), na.rm = TRUE),

        Count = n(),
        Percent = n() / total_sim,

        # Coverage probability
        Coverage = if (true_col_mode) {
          mean(
            (.data[[param]] >= (.data[[est_col]] - z * .data[[theo_col]]^0.5)) &
              (.data[[param]] <= (.data[[est_col]] + z * .data[[theo_col]]^0.5)),
            na.rm = TRUE
          )
        } else {
          mean(
            (true_value_const >= (.data[[est_col]] - z * .data[[theo_col]])) &
              (true_value_const <= (.data[[est_col]] + z * .data[[theo_col]])),
            na.rm = TRUE
          )
        },

        .groups = "drop"
      )

    results_list[[param]] <- res
  }

  return(bind_rows(results_list))
}
