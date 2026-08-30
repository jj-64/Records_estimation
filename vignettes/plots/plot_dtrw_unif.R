library(readxl)
library(tidyverse)
library(patchwork)
library(scales)
library(ggplot2)
library(dplyr)
path_out = "vignettes/plots_output/dtrw/unif/"

#___________________________________
# Read data ---------
#___________________________________
true_param = 1

df_biased <- read_excel(paste0("data/param_est/dtrw/unif/Simulation_dtrw_Xt_Rn_overT_unif_max=",true_param,"_biased.xlsx"), sheet = "All Par")
df_record_biased <- read_excel(paste0("data/param_est/dtrw/unif/Simulation_dtrw_Xt_Rn_overT_unif_max=",true_param,"_biased.xlsx"), sheet = "All Par_record")
df <- read_excel(paste0("data/param_est/dtrw/unif/Simulation_dtrw_Xt_Rn_overT_unif_max=",true_param,"_BiasCorrected_simulation_5000.xlsx"), sheet = "All Par")
df_record <- read_excel(paste0("data/param_est/dtrw/unif/Simulation_dtrw_Xt_Rn_overT_unif_max=",true_param,"_BiasCorrected_simulation_5000.xlsx"), sheet = "All Par_record")

colnames(df_biased)[1] = "T"
colnames(df_record_biased)[1] = "T"
colnames(df)[1] = "T"
colnames(df_record)[1] = "T"

#___________________________________
# Extract numeric T ----------
#___________________________________

df_plot <- bind_rows(
  df %>% mutate(Scenario = "All observations"),
  df_record %>% mutate(Scenario = "Records"),
  df_biased %>% mutate(Scenario = "All observations (biased)"),
  df_record_biased %>% mutate(Scenario = "Records (biased)")
) %>%
  mutate(T_num = as.numeric(gsub("T=", "", T)))

#___________________________________
# 1. Bias plot --------
#___________________________________

# Bias plot
bias_plot <- ggplot(df_plot,
                    aes(x = T_num,
                        y = AVG_bias_max,
                        color = Scenario,
                        linetype = Scenario,
                        shape = Scenario,
                        group = Scenario)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "black") +
  scale_color_manual(values = c(
    "All observations" = "#1f77b4",
    "Records" = "#009E73",
    "All observations (biased)" = "#001f77",
    "Records (biased)" = "#006400"
  )) +
  scale_linetype_manual(values = c(
    "All observations" = "solid",
    "Records" = "dashed",
    "All observations (biased)" = "solid",
    "Records (biased)" = "dashed"
  )) +
  scale_shape_manual(values = c(
    "All observations" = 16,
    "Records" = 17,
    "All observations (biased)" = 18,
    "Records (biased)" = 19
  )) +
  theme_bw() +
  labs(
    title = "", #Average Bias across T",
    x = "series length (T)",
    y = "Bias",
    color = "",
    linetype = "",
    shape = ""
  ) +  theme(legend.position = "none")

bias_plot

#___________________________________
# 2. Coverage Probability ----------
#___________________________________

coverage_plot <- ggplot(df_plot,
                        aes(x = T_num,
                            y = Coverage_proba_max,
                            color = Scenario,
                            linetype = Scenario,
                            shape = Scenario,
                            group = Scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.95,
             linetype = "dotted",
             color = "black") +
  scale_color_manual(values = c(
    "All observations" = "#1f77b4",
    "Records" = "#009E73",
    "All observations (biased)" = "#001f77",
    "Records (biased)" = "#006400"
  )) +
  scale_linetype_manual(values = c(
    "All observations" = "solid",
    "Records" = "dashed",
    "All observations (biased)" = "solid",
    "Records (biased)" = "dashed"
  )) +
  scale_shape_manual(values = c(
    "All observations" = 16,
    "Records" = 17,
    "All observations (biased)" = 18,
    "Records (biased)" = 19
  ))  +
  theme_bw() +
  labs(
    title = "", #"Coverage Probability across T",
    x = "series length (T)",
    y = "Coverage Probability",
    color = "",
    linetype = "",
    shape = ""
  )+  theme(legend.position = "none")

coverage_plot

#___________________________________
# 3. Empirical vs Theoretical Variance ---------
#___________________________________

df_plot2 <- bind_rows(
  df %>% mutate(Scenario = "All observations"),
  df_record %>% mutate(Scenario = "Records"),
  df_biased %>% mutate(Scenario = "All observations (biased)"),
  df_record_biased %>% mutate(Scenario = "Records (biased)")
) %>%
  mutate(T_num = as.numeric(gsub("T=", "", T))) %>%
  pivot_longer(
    cols = c(AVG_Emp_var_max, AVG_Theo_var_max),
    names_to = "Variance_Scenario",
    values_to = "Variance"
  ) %>%
  mutate(
    Variance_Scenario = recode(
      Variance_Scenario,
      AVG_Emp_var_max = "Empirical variance",
      AVG_Theo_var_max = "Theoretical variance"
    ),
    Line = paste(Scenario, Variance_Scenario, sep = " - ")
  )

variance_plot <- ggplot(
  df_plot2,
  aes(
    x = T_num,
    y = Variance,
    color = Scenario,
    linetype = Variance_Scenario,
    shape = Variance_Scenario,
    group = Line
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "All observations" = "#1f77b4",
    "Records" = "#009E73",
    "All observations (biased)" = "#001f77",
    "Records (biased)" = "#006400"
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
    title = "", #"Empirical and Theoretical Variance across T",
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

# _________________________________
# 8. RMSE -------
# ___________________________________

rmse_data <-
  df_plot %>% mutate(
    RMSE_max = AVG_bias_max^2 + AVG_Emp_var_max) %>%
  mutate(T_num = as.numeric(gsub("T=", "", T))) %>%
  pivot_longer(
    cols = c(RMSE_max),
    names_to = "Parameter",
    values_to = "RMSE"
  ) %>%
  mutate(
    Parameter = recode(
      Parameter,
      RMSE_max = "Max"
    )
  )

rmse_plot <- ggplot(
  rmse_data,
  aes(
    x = T_num,
    y = RMSE,
    color = Scenario,
    linetype = Scenario,
    shape = Scenario,
    group = Scenario
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  # geom_hline(yintercept = 0.95,
  #            linetype = "dashed",
  #            color = "black") +
  # facet_wrap(~Parameter, scales = "free_y") +
  theme_bw() +
  labs(
    title = "", #"Coverage Probability",
    x = "series length (T)",
    y = "RMSE"
  ) +
  scale_color_manual(values = c(
    "All observations" = "#1f77b4",
    "Records" = "#009E73",
    "All observations (biased)" = "#001f77",
    "Records (biased)" = "#006400"
  )) +
  scale_linetype_manual(values = c(
    "All observations" = "solid",
    "Records" = "dashed",
    "All observations (biased)" = "solid",
    "Records (biased)" = "dashed"
  )) +
  scale_shape_manual(values = c(
    "All observations" = 16,
    "Records" = 17,
    "All observations (biased)" = 18,
    "Records (biased)" = 19
  ))+
  theme(legend.position = "none")

rmse_plot


#___________________________________
# 6. Average Number of Records  -----------
#___________________________________

NT_plot <- ggplot(df_plot %>% filter(Scenario == "Records"),
    aes(x = T_num,
       y = AVG_param_N_T
       )
    ) +
    geom_line(linewidth = 1, linetype = "dashed", color = "#E69F00") +
    geom_point(size = 3, color = "#E69F00") +
    scale_y_continuous(breaks = seq(10,20,2)) +
    theme_bw() +
    labs(
      title = "", #"Average Number of Observed Records",
      x = "series length (T)",
      y = expression(bar(N)[T])
    )
NT_plot

#___________________________________
# Save all figures -----------
#___________________________________

ggsave(paste0(path_out,"Bias_plot.png"),
       bias_plot,
       width = 10,
       height = 6,
       dpi = 300)

ggsave(paste0(path_out,"Coverage_plot.png"),
       coverage_plot,
       width = 10,
       height = 6,
       dpi = 300)

ggsave(paste0(path_out,"Variance_plot.png"),
       variance_plot,
       width = 12,
       height = 8,
       dpi = 300)

ggsave(paste0(path_out,"RMSE_plot.png"),
       rmse_plot,
       width = 8,
       height = 5,
       dpi = 300)

ggsave(paste0(path_out,"Avg_NT_plot.png"),
       NT_plot,
       width = 8,
       height = 5,
       dpi = 300)

# Common Springer theme --------
theme_extremes <- theme_bw(base_family = "Times New Roman") +
  theme(
    text = element_text(size = 11),
    plot.title = element_text(
      hjust = 0.5,
      size = 11,
      family = "Times New Roman"
    ),
    plot.margin = margin(10, 10, 10, 10)
  )


bias_plot <- bias_plot +
  labs(
    title = "Bias",
    tag = "(a)"
  ) +
  theme_extremes +
  theme(
    legend.position = "none",
    plot.tag.position = c(0.02, 0.97),
    plot.tag = element_text(
      family = "Times New Roman",
      size = 11,
      face = "plain"
    )
  )
coverage_plot <- coverage_plot +
  labs(
    title = "Coverage Probability",
    tag = "(b)"
  ) +
  theme_extremes +
  theme(
    legend.position = "right",
    plot.tag.position = c(0.02, 0.97),
    plot.tag = element_text(
      family = "Times New Roman",
      size = 11,
      face = "plain"
    )
  )

variance_plot <- variance_plot +
  labs(
    title = "Empirical and Theoretical Variance",
    tag = "(c)"
  ) +
  theme_extremes +
  theme(
    # Keep the legend here so patchwork can collect it
    legend.position = "right",
    plot.tag.position = c(0.02, 0.97),
    plot.title = element_text(
      hjust = 0.5,
      size = 11,
      family = "Times New Roman",
      margin = margin(b = 8)
    )
  )


rmse_plot <- rmse_plot +
  labs(
    title = "Root Mean Square Error",
    tag = "(c)"
  ) +
  theme_extremes +
  theme(
    legend.position = "none",
    plot.tag.position = c(0.02, 0.97),
    plot.tag = element_text(
      family = "Times New Roman",
      size = 11,
      face = "plain"
    )
  )

NT_plot <- NT_plot +
  labs(
    title = "Average Number of Observed Records",
    tag = "(d)"
  ) +
  theme_extremes +
  theme(
    legend.position = "none",
    plot.tag.position = c(0.02, 0.97),
    plot.tag = element_text(
      family = "Times New Roman",
      size = 11,
      face = "plain"
    )
  )

# Combine with one common legend
combined_plot <- (bias_plot | coverage_plot) /
  (rmse_plot |  NT_plot)  + # variance_plot
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  )

combined_plot

ggsave(paste0(path_out,"dtrw_unif.pdf"),
       combined_plot,
       width = 9,
       height = 9,
       device = cairo_pdf)

ggsave(paste0(path_out,"dtrw_unif.png"),
       combined_plot,
       width = 12,
       height = 9,
       dpi = 300)
