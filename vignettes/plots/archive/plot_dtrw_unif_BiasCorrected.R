library(readxl)
library(tidyverse)
library(patchwork)
library(scales)
library(ggplot2)
library(dplyr)

#___________________________________
# Read data ---------
#___________________________________

df <- read_excel("data/param_est/dtrw/unif/Simulation_dtrw_Xt_Rn_overT_unif_max=1_biasCorrected.xlsx", sheet = "All Par")
df_record <- read_excel("data/param_est/dtrw/unif/Simulation_dtrw_Xt_Rn_overT_unif_max=1_biasCorrected.xlsx", sheet = "All Par_record")

colnames(df)[1] = "T"
colnames(df_record)[1] = "T"

#___________________________________
# Extract numeric T ----------
#___________________________________

df_plot <- bind_rows(
  df %>% mutate(Type = "All observations"),
  df_record %>% mutate(Type = "Records")
) %>%
  mutate(T_num = as.numeric(gsub("T=", "", T)))

#___________________________________
# 1. Bias plot --------
#___________________________________

# Bias plot
bias_plot <- ggplot(df_plot,
                    aes(x = T_num,
                        y = AVG_bias_max,
                        color = Type,
                        linetype = Type,
                        shape = Type,
                        group = Type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "black") +
  scale_color_manual(values = c(
    "All observations" = "#1f77b4",
    "Records" = "#d62728"
  )) +
  scale_linetype_manual(values = c(
    "All observations" = "solid",
    "Records" = "dashed"
  )) +
  scale_shape_manual(values = c(
    "All observations" = 16,
    "Records" = 17
  )) +
  theme_bw() +
  labs(
    title = "Average Bias across T",
    x = "series length (T)",
    y = "Bias",
    color = "",
    linetype = "",
    shape = ""
  )

bias_plot

#___________________________________
# 2. Coverage Probability ----------
#___________________________________

coverage_plot <- ggplot(df_plot,
                        aes(x = T_num,
                            y = Coverage_proba_max,
                            color = Type,
                            linetype = Type,
                            shape = Type,
                            group = Type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.95,
             linetype = "dotted",
             color = "black") +
  scale_color_manual(values = c(
    "All observations" = "#1f77b4",
    "Records" = "#d62728"
  )) +
  scale_linetype_manual(values = c(
    "All observations" = "solid",
    "Records" = "dashed"
  )) +
  scale_shape_manual(values = c(
    "All observations" = 16,
    "Records" = 17
  )) +
  theme_bw() +
  labs(
    title = "Coverage Probability across T",
    x = "series length (T)",
    y = "Coverage Probability",
    color = "",
    linetype = "",
    shape = ""
  )

coverage_plot

#___________________________________
# 3. Empirical vs Theoretical Variance ---------
#___________________________________


df_plot2 <- bind_rows(
  df %>% mutate(Scenario = "All observations"),
  df_record %>% mutate(Scenario = "Records")
) %>%
  mutate(T_num = as.numeric(gsub("T=", "", T))) %>%
  pivot_longer(
    cols = c(AVG_Emp_var_max, AVG_Theo_var_max),
    names_to = "Variance_Type",
    values_to = "Variance"
  ) %>%
  mutate(
    Variance_Type = recode(
      Variance_Type,
      AVG_Emp_var_max = "Empirical variance",
      AVG_Theo_var_max = "Theoretical variance"
    ),
    Line = paste(Scenario, Variance_Type, sep = " - ")
  )

variance_plot <- ggplot(
  df_plot2,
  aes(
    x = T_num,
    y = Variance,
    color = Scenario,
    linetype = Variance_Type,
    shape = Variance_Type,
    group = Line
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "All observations" = "#1f77b4",
    "Records" = "#d62728"
  )) +
  scale_linetype_manual(values = c(
    "Empirical variance" = "solid",
    "Theoretical variance" = "dashed"
  )) +
  scale_shape_manual(values = c(
    "Empirical variance" = 16,
    "Theoretical variance" = 17
  )) +
  theme_bw() +
  labs(
    title = "Empirical and Theoretical Variance across T",
    x = "series length (T)",
    y = "Variance",
    color = "Scenario",
    linetype = "Variance type",
    shape = "Variance type"
  )

variance_plot

#___________________________________
# 4. Ratio Empirical/Theoretical Variance ---------
#___________________________________
#
# ratio_df <- df_long %>%
#   filter(Metric %in% c("AVG_Emp_var","AVG_Theo_var")) %>%
#   pivot_wider(
#     names_from = Metric,
#     values_from = Value
#   ) %>%
#   mutate(
#     Ratio = AVG_Emp_var / AVG_Theo_var
#   )
#
# ratio_plot <- ggplot(
#   ratio_df,
#   aes(x = gamma_obs,
#       y = Ratio,
#       colour = T,
#       group = T)
# ) +
#   geom_line(size = 1) +
#   geom_point(size = 2) +
#   geom_hline(yintercept = 1,
#              linetype = "dashed") +
#   facet_wrap(~Parameter,
#              scales = "free_y") +
#   theme_bw() +
#   labs(
#     title = "Empirical / Theoretical Variance",
#     x = expression(gamma),
#     y = "Ratio"
#   )

#___________________________________
# 5. Average Parameter Estimates-----
#___________________________________

# avg_param_plot <- df_long %>%
#   filter(Metric == "AVG_param") %>%
#   ggplot(aes(x = gamma_obs,
#              y = Value,
#              color = T,
#              group = T)) +
#   geom_line(size = 1) +
#   geom_point(size = 2) +
#   facet_wrap(~Parameter,
#              scales = "free_y") +
#   theme_bw() +
#   labs(
#     title = "Average Estimated Parameters",
#     x = expression(gamma),
#     y = "Estimate"
#   )

#___________________________________
# 6. Average Number of Records
#___________________________________

NT_plot <- ggplot(df_plot %>% filter(Type == "Records"),
    aes(x = T_num,
       y = AVG_param_N_T
       )
    ) +
    geom_line(linewidth = 1, linetype = "dashed", color = "#d62728") +
    geom_point(size = 3, color = "#d62728") +
    scale_y_continuous(breaks = seq(10,20,2)) +
    theme_bw() +
    labs(
      title = "Average Number of Observed Records",
      x = "series length (T)",
      y = "Average Number of Observed Records"
    )
NT_plot

#___________________________________
# Save all figures
#___________________________________

ggsave("vignettes/plots_output/dtrw/unif/BiasCorrected/Bias_plot.png",
       bias_plot,
       width = 10,
       height = 6,
       dpi = 300)

ggsave("vignettes/plots_output/dtrw/unif/BiasCorrected/Coverage_plot.png",
       coverage_plot,
       width = 10,
       height = 6,
       dpi = 300)

ggsave("vignettes/plots_output/dtrw/unif/BiasCorrected/Variance_plot.png",
       variance_plot,
       width = 12,
       height = 8,
       dpi = 300)


ggsave("vignettes/plots_output/dtrw/unif/BiasCorrected/Avg_NT_plot.png",
       NT_plot,
       width = 8,
       height = 5,
       dpi = 300)
