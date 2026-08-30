library(readxl)
library(tidyverse)
library(patchwork)
library(scales)
library(ggplot2)
path_out = "vignettes/plots_output/ynm/ynm_frechet_inv_scale/"

#___________________________________
# Read data ---------
#___________________________________
simulation = 5000

df_full <- read_excel(paste0("data/param_est/ynm/frechet_inv_scale/Simulation_ynm_Xt_Rn_frechet_shape=2_scale=1_simulation_", simulation, ".xlsx"), sheet = "All par") %>% mutate(Scenario = "All observations")
df_records <- read_excel(paste0("data/param_est/ynm/frechet_inv_scale/Simulation_ynm_Xt_Rn_frechet_shape=2_scale=1_simulation_", simulation, ".xlsx"), sheet = "All par_record")%>% mutate(Scenario = "Records")

df = rbind(df_full, df_records)

df$T <- factor(df$T, levels = c("T=50","T=100","T=200"))

#___________________________________
# Long format ----------
#___________________________________

df_long <- df %>%
  pivot_longer(
    cols = c(gamma, scale, shape),
    names_to = "Parameter",
    values_to = "Value"
  )

# df_long <- df_long %>%
#   mutate(
#     Parameter = factor(
#       Parameter,
#       levels = c("gamma","scale","shape"),
#       labels = c(
#         expression(gamma),
#         expression(sigma),
#         expression(alpha)
#       )
#     )
#   )

param_labs <- c(
  gamma    = expression(gamma),
  scale    = expression(sigma),
  shape = expression(alpha)
)

# facet_wrap(
#   ~Parameter,
#   labeller = labeller(
#     Parameter = as_labeller(param_labs,
#                             label_parsed)
#   )
# )

#___________________________________
# 1. Bias plot --------
#___________________________________

# bias_plot <- df_long %>%
#   filter(Metric == "AVG_bias") %>%
#   ggplot(
#     aes(
#       x = gamma_obs,
#       y = Value,
#       colour = Scenario,
#       linetype = Scenario
#     )
#   ) +
#   geom_line(size = 1) +
#   geom_point(size = 2) +
#   facet_grid(Parameter ~ T,
#              scales = "free_y") +
#   theme_bw() +
#   labs(
#     title = "Bias Comparison",
#     x = expression(gamma),
#     y = "Bias"
#   )
#
# bias_plot


bias_df <- df_long %>%
  filter(Metric == "AVG_bias")

bias_plot <- ggplot( bias_df,
  aes(
    x = gamma_obs,
    y = Value,
    color = Scenario,
    linetype = Scenario,
    shape = Scenario,
    group = Scenario
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(
    Parameter ~ T,
    labeller = labeller(
      Parameter = as_labeller(
        c(gamma="gamma",
          scale="sigma",
          shape="alpha"),
        label_parsed
      )
    )
  ) +
  scale_colour_manual(
    values = c(
      "All observations" = "#0072B2",
      "Records" = "#009E73"
    )
    )+
  scale_linetype_manual(
    values = c(
    "All observations" = "solid",
    "Records" = "dashed"
  )) +
  scale_shape_manual(values = c(
    "All observations" = 16,
    "Records" = 17
  )) +
  theme_bw() +
  theme(
    legend.position = "none")+
  labs(
    title = "", #"Empirical and Theoretical Variance",
    x = expression(gamma),
    y = "Bias") +
  scale_x_continuous(n.breaks=3)

bias_plot

#___________________________________
# 2. Coverage Probability ----------
#___________________________________

coverage_plot = ggplot( df_long %>%
          filter(Metric == "Coverage_proba"),
        aes(
          x = gamma_obs,
          y = Value,
          color = Scenario,
          linetype = Scenario,
          shape = Scenario,
          group = Scenario
        )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(
    Parameter ~ T,
    labeller = labeller(
      Parameter = as_labeller(
        c(gamma="gamma",
          scale="sigma",
          shape="alpha"),
        label_parsed
      )
    )
  ) +
  scale_colour_manual(
    values = c(
      "All observations" = "#0072B2",
      "Records" = "#009E73"
    )
  )+
  scale_linetype_manual(
    values = c(
      "All observations" = "solid",
      "Records" = "dashed"
    )) +
  scale_shape_manual(values = c(
    "All observations" = 16,
    "Records" = 17
  )) +
  theme_bw() +
  labs(
    title = "", #"Empirical and Theoretical Variance",
    x = expression(gamma),
    y = "Coverage Probability") +
  theme(
    legend.position = "none")+
  scale_x_continuous(n.breaks=3)

coverage_plot

#___________________________________
# 3. Empirical vs Theoretical Variance ---------
#___________________________________

var_df <- df_long %>%
  filter(Metric %in% c("AVG_Emp_var","AVG_Theo_var"))

variance_plot <- var_df %>%
  ggplot(
    aes(
      x = gamma_obs,
      y = Value,
      colour = Scenario,
      linetype = Metric,
      shape = Metric,
      group = interaction(Scenario, Metric)
    )
  ) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_grid(Parameter ~ T,
             scales = "free_y",
             labeller = labeller(
               Parameter = as_labeller(
                 c(gamma="gamma",
                   scale="sigma",
                   shape="alpha"),
                 label_parsed
               )
             )) +
  labs(
    title = "", #"Empirical and Theoretical Variance",
    x = expression(gamma),
    y = "Variance",
    color = "Scenario",
    linetype = "Variance Type",
    shape = "Variance Type"
  ) +
  scale_linetype_manual(
    values = c(
      "AVG_Emp_var"  = "solid",
      "AVG_Theo_var" =  "dashed"
    ),
    ,
    labels = c(
      "Empirical Variance",
      "Theoretical Variance"
    )
  )+ scale_color_manual(values = c(
    "All observations" = "#0072B2",
    "Records" = "#009E73"
  )) +
  scale_shape_manual(
    values = c(
      "AVG_Emp_var"  = 16,
      "AVG_Theo_var" = 17
    ),
    labels = c(
      "Empirical Variance",
      "Theoretical Variance"
    )
    ) +
  theme_bw() +
  scale_x_continuous(n.breaks = 3)

variance_plot


#___________________________________
# 4. Relative Efficiency ---------
#___________________________________
#
# efficiency_df <- df_long %>%
#   filter(Metric == "AVG_Emp_var") %>%
#   dplyr::select(
#     T,
#     gamma_obs,
#     Parameter,
#     Scenario,
#     Value
#   ) %>%
#   pivot_wider(
#     names_from = Scenario,
#     values_from = Value
#   ) %>%
#   mutate(
#     RelativeEfficiency = Records / `All observations`
#   )
#
# eff_plot <- ggplot(
#   efficiency_df,
#   aes(
#     x = gamma_obs,
#     y = RelativeEfficiency,
#     colour = T,
#     group = T
#   )
# ) +
#   geom_line(size = 1) +
#   geom_point(size = 2) +
#   facet_wrap(~Parameter,
#              scales = "free_y") +
#   theme_bw() +
#   labs(
#     title = "Efficiency Loss Using Record Data",
#     x = expression(gamma),
#     y = expression(
#       frac(Var[Record], Var[`All observations`])
#     )
#   )
#
# eff_plot

# _________________________________
# 8. RMSE -------
# ___________________________________

rmse_data <- df %>% dplyr::select(-N_T) %>%
  filter(Metric %in% c("AVG_bias", "AVG_Emp_var")) %>%
  pivot_wider(
    names_from = Metric,
    values_from = c(gamma, scale, shape)
  ) %>%
  mutate(
    RMSE_gamma = sqrt(gamma_AVG_bias^2 + gamma_AVG_Emp_var),
    RMSE_scale = sqrt(scale_AVG_bias^2 + scale_AVG_Emp_var),
    RMSE_shape = sqrt(shape_AVG_bias^2 + shape_AVG_Emp_var)
  )%>%
  dplyr::select(T, gamma_obs,
                RMSE_gamma,
                RMSE_scale,
                RMSE_shape,
                Scenario) %>%
  pivot_longer(
    starts_with("RMSE_"),
    names_to = "Parameter",
    values_to = "RMSE"
  )


rmse_plot = ggplot( rmse_data ,
                    aes(
                      x = gamma_obs,
                      y = RMSE,
                      color = Scenario,
                      linetype = Scenario,
                      shape = Scenario,
                      group = Scenario
                    )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(
    Parameter ~ T,
    labeller = labeller(
      Parameter = as_labeller(
        c(RMSE_gamma="gamma",
          RMSE_scale="sigma",
          RMSE_shape="alpha"),
        label_parsed
      )
    )
  ) +
  scale_colour_manual(
    values = c(
      "All observations" = "#0072B2",
      "Records" = "#009E73"
    )
  )+
  scale_linetype_manual(
    values = c(
      "All observations" = "solid",
      "Records" = "dashed"
    )) +
  scale_shape_manual(values = c(
    "All observations" = 16,
    "Records" = 17
  )) +
  theme_bw() +
  theme(legend.position = "none")+
  labs(x= expression(gamma))

rmse_plot


#___________________________________
# 6. Average Number of Records -------
#___________________________________

NT_plot <- df %>%
  filter(N_T >0) %>%
  distinct(T,
           gamma_obs, N_T) %>%
  ggplot(
    aes(
      x = gamma_obs,
      y = N_T,
      colour = T,
      group = T,
      linetype = T
    )
  ) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    title = "", # "Average Number of Records",
    x = expression(gamma),
    y = expression(bar(N)[T])
  ) +
  scale_color_manual(values = c(
    "T=50" = "#0072B2",
    "T=100" = "#009E73",
    "T=200" = "#E69F00"
  )) +
  scale_x_continuous(n.breaks = 3)

NT_plot

#___________________________________
# 7. Relative efficiency -------
#___________________________________
#
# eff_NT <- efficiency_df %>%
#   left_join(
#     df %>% filter(N_T > 0 ) %>%
#       distinct(
#         T,
#         gamma_obs,
#         N_T
#       ),
#     by = c("T","gamma_obs")
#   )
#
# eff_NT = eff_NT[seq(1,nrow(eff_NT),b=2),]
#
# NT_eff_plot <- ggplot(
#   eff_NT,
#   aes(
#     x = N_T,
#     y = RelativeEfficiency,
#     colour = T
#   )
# ) +
#   geom_point(size = 3) +
#   geom_smooth(
#     Scenario = "lm",
#     se = FALSE
#   ) +
#   facet_wrap(~Parameter,
#              scales = "free_y") +
#   theme_bw() +
#   labs(
#     title = "Efficiency Loss versus Number of Records",
#     x = expression(bar(N)[T]),
#     y = expression(
#       frac(Var[Record], Var[Full])
#     )
#   ) +scale_colour_manual(
#     values = c(
#       "T=50"  = "#E69F00",
#       "T=100" = "#0072B2",
#       "T=200" = "#009E73"
#     )
#   )
#
# NT_eff_plot

#___________________________________
# Save all figures -----
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

# ggsave(paste0(path_out,"Variance_ratio_plot.png"),
#        ratio_plot,
#        width = 10,
#        height = 6,
#        dpi = 300)

# ggsave(paste0(path_out,"Parameter_plot.png"),
#        avg_param_plot,
#        width = 10,
#        height = 6,
#        dpi = 300)

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

# theme_extremes <- theme_bw(base_family = "Times New Roman") +
#   theme(
#     text = element_text(size = 11),
#     plot.title = element_text(
#       hjust = 0.5,
#       size = 11,
#       family = "Times New Roman"
#     ),
#     #legend.position = "bottom",
#     plot.margin = margin(10, 10, 10, 10)
#   )
#
# bias_plot <- bias_plot +
#   labs( title = "Bias",
#     tag = "(a)"
#     ) +
#   theme_extremes +
#   theme(
#   legend.position = "none",
#   plot.tag.position = c(0.02, 0.97),
#   plot.tag = element_text(
#     family = "Times New Roman",
#     size = 11,
#     face = "plain"
#   )
# )
#
# coverage_plot <- coverage_plot +
#   labs(
#     title = "Coverage Probability",
#     tag = "(b)"
#   ) +
#   theme_extremes +
#   theme(
#     legend.position = "none",
#     plot.tag.position = c(0.02, 0.97),
#     plot.tag = element_text(
#       family = "Times New Roman",
#       size = 11,
#       face = "plain"
#     )
#   )
#
# variance_plot <- variance_plot +
#   labs(
#     title = "Empirical and Theoretical Variance",
#     tag = "(c)"
#   ) +
#   theme_extremes +
#   theme(
#     legend.position = "bottom",
#     plot.tag.position = c(0.02, 0.97),
#     plot.title = element_text(
#       hjust = 0.5,
#       size = 11,
#       family = "Times New Roman",
#       margin = margin(b = 8)
#     )
#   )
#
#
# NT_plot <- NT_plot +
#   labs(
#     title = "Average Number of Observed Records",
#     tag = "(d)"
#   ) +
#   theme_extremes +
#   theme(
#     legend.position = "bottom",
#     plot.tag.position = c(0.02, 0.97),
#     plot.tag = element_text(
#       family = "Times New Roman",
#       size = 11,
#       face = "plain"
#     )
#   )
#
# # Combine with one common legend
# combined_plot <- (bias_plot | coverage_plot) /
#   (variance_plot | NT_plot) +
#   #plot_layout(guides = "collect") &
#   theme(
#     legend.position = "bottom",
#     legend.justification = "center"
#   )
#
# combined_plot

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

# Combine plots and collect the legend globally
combined_plot <- (bias_plot | coverage_plot) /
  (rmse_plot |  NT_plot)  + # variance_plot
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  )

combined_plot

ggsave(paste0(path_out,"ynm_frechet.pdf"),
       combined_plot,
       width = 9,
       height = 9,
       device = cairo_pdf)

ggsave(paste0(path_out,"ynm_frechet.png"),
       combined_plot,
       width = 12,
       height = 9,
       dpi = 300)
