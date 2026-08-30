library(dplyr)
library(tidyr)

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
    select(Model, Distribution, Parameter,
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
