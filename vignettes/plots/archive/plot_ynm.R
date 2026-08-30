library(readxl)
library(tidyverse)
library(patchwork)
library(scales)

#___________________________________
# Read data ---------
#___________________________________

df <- read_excel("data/param_est/ynm/gumbel/Simulation_ynm_Xt_Rn_gumbel_location=0_scale=1.xlsx", sheet = "All par")

df$T <- factor(df$T, levels = c("T=50","T=100","T=200"))

#___________________________________
# Long format ----------
#___________________________________

df_long <- df %>%
  pivot_longer(
    cols = c(gamma, scale, location),
    names_to = "Parameter",
    values_to = "Value"
  )

#___________________________________
# 1. Bias plot --------
#___________________________________

bias_plot <- df_long %>%
  filter(Metric == "AVG_bias") %>%
  ggplot(aes(x = gamma_obs,
             y = Value,
             color = T,
             group = T)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~Parameter, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Average Bias",
    x = expression(gamma),
    y = "Bias"
  )

bias_plot
#___________________________________
# 2. Coverage Probability ----------
#___________________________________

coverage_plot <- df_long %>%
  filter(Metric == "Coverage_proba") %>%
  ggplot(aes(x = gamma_obs,
             y = Value,
             color = T,
             group = T)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.95,
             linetype = "dashed",
             colour = "red") +
  facet_wrap(~Parameter) +
  theme_bw() +
  labs(
    title = "Coverage Probability",
    x = expression(gamma),
    y = "Coverage"
  )

#___________________________________
# 3. Empirical vs Theoretical Variance ---------
#___________________________________

var_df <- df_long %>%
  filter(Metric %in% c("AVG_Emp_var","AVG_Theo_var"))

variance_plot <- ggplot(
  var_df,
  aes(x = gamma_obs,
      y = Value,
      color = Metric,
      linetype = Metric)
) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_grid(Parameter ~ T,
             scales = "free_y") +
  theme_bw() +
  labs(
    title = "Empirical vs Theoretical Variance",
    x = expression(gamma),
    y = "Variance"
  )

#___________________________________
# 4. Ratio Empirical/Theoretical Variance ---------
#___________________________________

ratio_df <- df_long %>%
  filter(Metric %in% c("AVG_Emp_var","AVG_Theo_var")) %>%
  pivot_wider(
    names_from = Metric,
    values_from = Value
  ) %>%
  mutate(
    Ratio = AVG_Emp_var / AVG_Theo_var
  )

ratio_plot <- ggplot(
  ratio_df,
  aes(x = gamma_obs,
      y = Ratio,
      colour = T,
      group = T)
) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1,
             linetype = "dashed") +
  facet_wrap(~Parameter,
             scales = "free_y") +
  theme_bw() +
  labs(
    title = "Empirical / Theoretical Variance",
    x = expression(gamma),
    y = "Ratio"
  )

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

NT_plot <- ggplot(
  df %>%
    distinct(T, gamma_obs, N_T),
  aes(x = gamma_obs,
      y = N_T,
      color = T,
      group = T)
) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    title = "Average Number of Observed Records",
    x = expression(gamma),
    y = expression(bar(N)[T])
  )

#___________________________________
# Save all figures
#___________________________________

ggsave("Bias_plot.png",
       bias_plot,
       width = 10,
       height = 6,
       dpi = 300)

ggsave("Coverage_plot.png",
       coverage_plot,
       width = 10,
       height = 6,
       dpi = 300)

ggsave("Variance_plot.png",
       variance_plot,
       width = 12,
       height = 8,
       dpi = 300)

ggsave("Variance_ratio_plot.png",
       ratio_plot,
       width = 10,
       height = 6,
       dpi = 300)

ggsave("Parameter_plot.png",
       avg_param_plot,
       width = 10,
       height = 6,
       dpi = 300)

ggsave("Avg_NT_plot.png",
       NT_plot,
       width = 8,
       height = 5,
       dpi = 300)
