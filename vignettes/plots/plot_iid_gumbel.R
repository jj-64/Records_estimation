library(readxl)
library(tidyverse)
library(patchwork)
library(scales)
library(dplyr)
library(tidyr)
library(ggplot2)
path_out = "vignettes/plots_output/iid/gumbel/"
#___________________________________
# Read data ---------
#___________________________________

df <- read_excel("data/param_est/iid/gumbel/Simulation_iid_Xt_Rn_overT_gumbel_scale=1_location=0_simulation_5000.xlsx", sheet = "All Par") %>% mutate(Scenario = "All observations")
df_record <- read_excel("data/param_est/iid/gumbel/Simulation_iid_Xt_Rn_overT_gumbel_scale=1_location=0_simulation_5000.xlsx", sheet = "All Par_record")%>% mutate(Scenario = "Records")

colnames(df)[1] = "T"
colnames(df_record)[1] = "T"

df_long =  bind_rows(
  df,
  df_record
)

param_labs <- c(
  scale    = expression(sigma),
  location = expression(mu)
)
#___________________________________
# Long format ----------
#___________________________________

bias_df <- df_long %>%
  mutate(T_num = as.numeric(gsub("T=", "", T))) %>%
  pivot_longer(
    cols = c(AVG_bias_location, AVG_bias_scale),
    names_to = "Parameter",
    values_to = "Bias"
  ) %>%
  mutate(
    Parameter = recode(
      Parameter,
      AVG_bias_location = "location",
      AVG_bias_scale = "Scale"
    )
  )

#___________________________________
# 1. Bias plot --------
#___________________________________

bias_plot <- ggplot(
  bias_df,
  aes(
    x = T_num,
    y = Bias,
    color = Scenario,
    linetype = Scenario,
    shape = Scenario,
    group = Scenario
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0,
             linetype = "dotted",
             color = "black") +
  facet_wrap(~Parameter, scales = "free_y",
             labeller = labeller(
               Parameter = as_labeller(
                 c(Scale="sigma",
                   location="mu"),
                 label_parsed
               )
             )) +
  theme_bw() +
  labs(
    title = "", #"Average Bias",
    x = "series length (T)",
    y = "Bias"
  ) +
  scale_color_manual(values = c(
    "All observations" = "#0072B2",
    "Records" = "#009E73"
  )) +
  scale_linetype_manual(values = c(
    "All observations" = "solid",
    "Records" = "dashed"
  )) +
  scale_shape_manual(values = c(
    "All observations" = 16,
    "Records" = 17
  ))+
  theme(legend.position = "none")

bias_plot

#___________________________________
# 2. Coverage Probability ----------
#___________________________________

coverage_df <- df_long %>%
  mutate(T_num = as.numeric(gsub("T=", "", T))) %>%
  pivot_longer(
    cols = c(Coverage_proba_location,
             Coverage_proba_scale),
    names_to = "Parameter",
    values_to = "Coverage"
  ) %>%
  mutate(
    Parameter = recode(
      Parameter,
      Coverage_proba_location = "location",
      Coverage_proba_scale = "Scale"
    )
  )

coverage_plot <- ggplot(
  coverage_df,
  aes(
    x = T_num,
    y = Coverage,
    color = Scenario,
    linetype = Scenario,
    shape = Scenario,
    group = Scenario
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.95,
             linetype = "dashed",
             color = "black") +
  facet_wrap(~Parameter,
             labeller = labeller(
               Parameter = as_labeller(
                 c(
                   Scale="sigma",
                   location="mu"),
                 label_parsed
               )
             )) +
  theme_bw() +
  labs(
    title = "", #"Coverage Probability",
    x = "series length (T)",
    y = "Coverage Probability"
  ) +
  scale_color_manual(values = c(
    "All observations" = "#0072B2",
    "Records" = "#009E73"
  )) +
  scale_linetype_manual(values = c(
    "All observations" = "solid",
    "Records" = "dashed"
  )) +
  scale_shape_manual(values = c(
    "All observations" = 16,
    "Records" = 17
  ))+
  theme(legend.position = "none")

coverage_plot

#___________________________________
# 3. Empirical vs Theoretical Variance ---------
#___________________________________

var_df <- df_long%>%
  mutate(T_num = as.numeric(gsub("T=", "", T))) %>%
  pivot_longer(
    cols = c(
      AVG_Emp_var_location,
      AVG_Theo_var_location,
      AVG_Emp_var_scale,
      AVG_Theo_var_scale
    ),
    names_to = "Variable",
    values_to = "Variance"
  ) %>%
  mutate(
    Parameter = ifelse(grepl("location", Variable),
                       "Location",
                       "Scale"),
    Variance_Type = ifelse(grepl("Emp", Variable),
                           "Empirical",
                           "Theoretical")
  )

variance_plot <- ggplot(
  var_df,
  aes(
    x = T_num,
    y = Variance,
    color = Scenario,
    linetype = Variance_Type,
    shape = Variance_Type,
    group = interaction(Scenario, Variance_Type)
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~Parameter, scales = "free_y",
             labeller = labeller(
               Parameter = as_labeller(
                 c(Scale="sigma",
                   Location = "mu"),
                 label_parsed
               )
             )) +
  scale_color_manual(values = c(
    "All observations" = "#0072B2",
    "Records" = "#009E73"
  )) +
  scale_linetype_manual(values = c(
    "Empirical" = "solid",
    "Theoretical" = "dashed"
  )) +
  scale_shape_manual(values = c(
    "Empirical" = 16,
    "Theoretical" = 17
  )) +
  theme_bw() +
  labs(
    title = "Empirical and Theoretical Variance",
    x = "series length (T)",
    y = "Variance",
    color = "Scenario",
    linetype = "Variance Type",
    shape = "Variance Type"
  )

variance_plot

# #___________________________________
# # 4. Ratio Empirical/Theoretical Variance ---------
# #___________________________________
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

rmse_data <- bind_rows(
  df %>% mutate(
    RMSE_scale = AVG_bias_scale^2 + AVG_Emp_var_scale,
    RMSE_location = AVG_bias_location^2 + AVG_Emp_var_location),
  df_record %>% mutate(
    RMSE_scale = AVG_bias_scale^2 + AVG_Emp_var_scale,
    RMSE_location = AVG_bias_location^2 + AVG_Emp_var_location)
) %>%
  mutate(T_num = as.numeric(gsub("T=", "", T))) %>%
  pivot_longer(
    cols = c(RMSE_location,
             RMSE_scale),
    names_to = "Parameter",
    values_to = "RMSE"
  ) %>%
  mutate(
    Parameter = recode(
      Parameter,
      RMSE_scale = "Scale",
      RMSE_location = "location"
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
  geom_hline(yintercept = 0.95,
             linetype = "dashed",
             color = "black") +
  facet_wrap(~Parameter,
             labeller = labeller(
               Parameter = as_labeller(
                 c(Scale="sigma",
                   location = "mu"),
                 label_parsed
               )
             )) +
  theme_bw() +
  labs(
    title = "", #"Coverage Probability",
    x = "series length (T)",
    y = "RMSE"
  ) +
  scale_color_manual(values = c(
    "All observations" = "#0072B2",
    "Records" = "#009E73"
  )) +
  scale_linetype_manual(values = c(
    "All observations" = "solid",
    "Records" = "dashed"
  )) +
  scale_shape_manual(values = c(
    "All observations" = 16,
    "Records" = 17
  ))+
  theme(legend.position = "none")

rmse_plot

#___________________________________
# 6. Average Number of Records
#___________________________________

NT_plot <- ggplot(
  bias_df %>%
    filter(Scenario == "Records" & Parameter == "Scale"),
  aes(x = T_num,
      y = AVG_param_N_T)
) +
  geom_line(size = 1,
            color = "#E69F00") +
  geom_point(size = 3,
             color = "#E69F00") +
  theme_bw() +
  labs(
    title = "", #"Average Number of Observed Records",
    x = "series length (T)",
    y = expression(bar(N)[T])
  )

NT_plot

#___________________________________
# Save all figures ----
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
# theme_extremes <- theme_bw(base_family = "Times New Roman") +
#   theme(
#     text = element_text(size = 11),
#     plot.tag = element_text(size = 11, family = "Times New Roman"),
#     legend.position = "bottom"
#   )

theme_extremes <- theme_bw(base_family = "Times New Roman") +
  theme(
    text = element_text(size = 11),
    plot.title = element_text(
      hjust = 0.5,
      size = 11,
      family = "Times New Roman"
    ),
    #legend.position = "bottom",
    plot.margin = margin(10, 10, 10, 10)
  )

bias_plot <- bias_plot +
  labs( title = "Bias",
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
    legend.position = "none",
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
    legend.position = "bottom",
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
    legend.position = "bottom",
    plot.tag.position = c(0.02, 0.97),
    plot.tag = element_text(
      family = "Times New Roman",
      size = 11,
      face = "plain"
    )
  )

# Combine with one common legend
combined_plot <- (bias_plot | coverage_plot) /
  (rmse_plot | NT_plot) + #
  #plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  )

combined_plot

ggsave(paste0(path_out,"iid_gumbel.pdf"),
       combined_plot,
       width = 9,
       height = 9,
       device = cairo_pdf)

ggsave(paste0(path_out,"iid_gumbel.png"),
       combined_plot,
       width = 12,
       height = 9,
       dpi = 300)
