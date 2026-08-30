library(ggplot2)
library(dplyr)
library(readxl)

## Functions----------
plot_estimation_results <- function(data,
                                    x_var = "NT",
                                    facet_by = c("Model", "Distribution"),
                                    param_filter = NULL) {

  df <- data

  # Optional filter on parameters
  if (!is.null(param_filter)) {
    df <- df %>% filter(Parameter %in% param_filter)
  }

  # _______________________________
  # 1. Bias plot
  # _______________________________
  p_bias <- ggplot(df, aes_string(x = x_var, y = "Bias",
                                  color = "Method", linetype = "Method")) +
    geom_line(size = 1) +
    geom_point() +
    facet_grid(as.formula(paste(paste(facet_by, collapse = "~"), "~ Parameter"))) +
    labs(title = "Bias Comparison",
         y = "Bias",
         x = x_var) +
    theme_minimal()


  # _______________________________
  # 2. Variance comparison
  # _______________________________
  df_var <- df %>%
    mutate(VarRatio = EmpVar / lag(EmpVar[Method == "Full"], default = NA))

  p_var <- ggplot(df, aes_string(x = x_var, y = "EmpVar",
                                 color = "Method", linetype = "Method")) +
    geom_line(size = 1) +
    geom_point() +
    facet_grid(as.formula(paste(paste(facet_by, collapse = "~"), "~ Parameter"))) +
    labs(title = "Empirical Variance Comparison",
         y = "Empirical Variance",
         x = x_var) +
    theme_minimal()


  # _______________________________
  # 3. Empirical vs Theoretical variance
  # _______________________________
  p_var_diag <- ggplot(df, aes(x = TheoVar, y = EmpVar, color = Method)) +
    geom_point(alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    facet_grid(as.formula(paste(paste(facet_by, collapse = "~"), "~ Parameter"))) +
    labs(title = "Empirical vs Theoretical Variance",
         x = "Theoretical Variance",
         y = "Empirical Variance") +
    theme_minimal()


  # _______________________________
  # 4. Coverage probability
  # _______________________________
  p_cov <- ggplot(df, aes_string(x = x_var, y = "Coverage",
                                 color = "Method", linetype = "Method")) +
    geom_line(size = 1) +
    geom_point() +
    geom_hline(yintercept = 0.95, linetype = "dotted", color = "black") +
    facet_grid(as.formula(paste(paste(facet_by, collapse = "~"), "~ Parameter"))) +
    labs(title = "Coverage Probability",
         y = "Coverage",
         x = x_var) +
    theme_minimal()


  # _______________________________
  # 5. Efficiency ratio plot
  # _______________________________
  df_ratio <- df %>%
    dplyr::select(Model, Distribution, Parameter, Method, !!sym(x_var), EmpVar) %>%
    tidyr::pivot_wider(names_from = Method, values_from = EmpVar) %>%
    mutate(EfficiencyRatio = Records / Full)

  p_eff <- ggplot(df_ratio,
                  aes_string(x = x_var, y = "EfficiencyRatio")) +
    geom_line(color = "red", size = 1) +
    geom_point(color = "red") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_grid(as.formula(paste(paste(facet_by, collapse = "~"), "~ Parameter"))) +
    labs(title = "Variance Ratio (Records / Full Data)",
         y = "Efficiency Ratio",
         x = x_var) +
    theme_minimal()


  return(list(
    BiasPlot = p_bias,
    VariancePlot = p_var,
    VarDiagnostic = p_var_diag,
    CoveragePlot = p_cov,
    EfficiencyPlot = p_eff
  ))
}

# function to save plot
save_plots <- function(plot_list, prefix = "results") {
  for (name in names(plot_list)) {
    ggsave(paste0(prefix, "_", name, ".png"),
           plot_list[[name]],
           width = 8, height = 6)
  }
}


## Data ------------
## Parameters
model = "DTRW"
dist = "norm"
parameter = "sd"
parameter_obs = 1

## Read Data
df_records_raw <- read_excel("data/param_est/dtrw/norm/Simulation_dtrw_Xt_Rn_overT_norm_sd=1.xlsx",
                                                    sheet = "summary_NT_record") %>% as.data.frame()
df_all_raw <- read_excel("data/param_est/dtrw/norm/Simulation_dtrw_Xt_Rn_overT_norm_sd=1.xlsx",
                             sheet = "summary_NT") %>% as.data.frame()


# Convert to required format
df_records <- df_records_raw %>%
  transmute(
    Model = model,
    Distribution = dist,

    Parameter = parameter,               # parameter name
    Method = "Records",
    T = T_val,
    NT = Param_N_T,

    Bias = Param_sd - parameter_obs,            # <-- replace 5 with true value

    EmpVar = Emp_sd^2,
    TheoVar = Theo_sd^2,

    Coverage = Proba_N_T             # or replace if true CP
  )

df_full <- df_all_raw %>%
  transmute(
    Model = model,
    Distribution = dist,

    Parameter = parameter,               # parameter name
    Method = "Full",
    T = T_val,
    NT = Param_N_T,

    Bias = Param_sd - parameter_obs,

    EmpVar = Emp_sd^2,
    TheoVar = Theo_sd^2,

    Coverage = Proba_N_T             # or replace if true CP
  )

## combine both
df <- bind_rows(df_records, df_full)

## How to use it
plot_estimation_results(df, param_filter = "sd")

# Display plots
plots$BiasPlot
plots$VariancePlot
plots$VarDiagnostic
plots$CoveragePlot
plots$EfficiencyPlot

## Summary Table -----------

create_summary_table_B <- function(data, avg_over = "T") {

  # Step 1: average over T (or number of records if different)
  df_avg <- data %>%
    group_by(Model, Distribution, Parameter, Method) %>%
    summarise(
      Bias = mean(Bias, na.rm = TRUE),
      EmpVar = mean(EmpVar, na.rm = TRUE),
      Coverage = mean(Coverage, na.rm = TRUE),
      .groups = "drop"
    )

  # Step 2: reshape to wide format (Full vs Records)
  df_wide <- df_avg %>%
    pivot_wider(
      names_from = Method,
      values_from = c(Bias, EmpVar, Coverage)
    )

  # Step 3: compute ratios and differences
  summary_table <- df_wide %>%
    mutate(
      VarRatio = EmpVar_Records / EmpVar_Full,
      BiasRatio = abs(Bias_Records) / abs(Bias_Full),
      CPDiff = Coverage_Records - Coverage_Full
    ) %>%
    dplyr::select(Model, Distribution, Parameter,
           VarRatio, BiasRatio, CPDiff)

  return(summary_table)
}

## call the function
table_B <- create_summary_table_B(df)
print(table_B)

# nice table
table_B %>%
  mutate(
    VarRatio = round(VarRatio, 3),
    BiasRatio = round(BiasRatio, 3),
    CPDiff = round(CPDiff, 3)
  )

## Sort to highlight worst information loss
table_B %>%
  arrange(desc(VarRatio))

## optional
# ✅ Interpretation guide (important for your paper)
#
# VarRatio > 1 → records lose efficiency
# BiasRatio > 1 → records introduce more bias
# CPDiff < 0 → undercoverage (bad inference)
#
# 👉 This table is your main quantitative evidence of information loss.
create_aggregated_table_B <- function(table_B) {
  table_B %>%
    group_by(Model, Distribution) %>%
    summarise(
      AvgVarRatio = mean(VarRatio, na.rm = TRUE),
      AvgBiasRatio = mean(BiasRatio, na.rm = TRUE),
      AvgCPDiff = mean(CPDiff, na.rm = TRUE),
      .groups = "drop"
    )
}
