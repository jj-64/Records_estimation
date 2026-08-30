# ==============================================================================
# Script: Consolidated Model Performance & Record Count Plots
# ==============================================================================

# 1. Load Required Packages ---------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# 2. Read and Prepare Data ----------------------------------------------------
t_values <- c(75, 100, 200)

read_summary_sheet <- function(t_val) {
  file_path <- sprintf("data/model_selection/nested/model_selection_exp_2_weibull_1.2_ldm_ldm_T=%d.xlsx", t_val)
  if (!file.exists(file_path)) return(NULL)

  read_excel(path = file_path, sheet = "summary") %>%
    mutate(T_size = factor(t_val, levels = t_values)) %>%
    filter(trend >= 0.1)
}

df_summary <- map_dfr(t_values, read_summary_sheet) %>%
  mutate(
    trend = as.numeric(trend),
    Type_I_Error = as.numeric(Type_I_Error),
    Power = as.numeric(Power),
    Nb_record = as.numeric(Nb_record)
  )

# Reshape data to long format for combining Power & Type I Error
df_rates <- df_summary %>%
  pivot_longer(
    cols = c(Type_I_Error, Power),
    names_to = "Metric",
    values_to = "Rate"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = c("Power", "Type_I_Error"), labels = c("Power", "Type I Error"))
  )

# 3. Styling Theme & Palettes --------------------------------------------------
pub_theme <- theme_bw(base_family = "Times New Roman") +
  theme(
    text = element_text(size = 11),
    panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title = element_text(),
    axis.text = element_text(color = "black"),
    #legend.position = "top",
    #legend.box = "vertical",
    #legend.title = element_text(size = 12),
    plot.title = element_text(
      hjust = 0.5,
      size = 11,
      family = "Times New Roman"
    )
  )

palette_t <- c("75" = "#0072B2", "100" = "#E69F00", "200" = "#009E73")

# 4. Generate Figures ----------------------------------------------------------

# Main Plot: Power and Type I Error together (distinguished by line type)
p_rates <- ggplot(df_rates, aes(x = trend, y = Rate, color = T_size, linetype = Metric)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = Metric), size = 2) +
  scale_color_manual(values = palette_t, name = "Sample Size (T)") +
  scale_linetype_manual(values = c("Power" = "solid", "Type I Error" = "dashed"), name = "Metric") +
  scale_shape_manual(values = c("Power" = 16, "Type I Error" = 17), name = "Metric") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = "Hypothesis Testing Performance",
    tag = "(a)",
    x = expression(gamma),
    y = "Type I Error | Power (%)"
  ) +
  pub_theme +
  theme(
    legend.position = "none",
    plot.tag.position = c(0.02, 0.97),
    plot.tag = element_text(
      family = "Times New Roman",
      size = 11,
      face = "plain"
    )
  )

# Side Plot: Number of Records
p_records <- ggplot(df_summary, aes(x = trend, y = Nb_record, color = T_size)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2, shape = 15) +
  scale_color_manual(values = palette_t, name = "Sample Size (T)") +
  labs(
    title = "Average Number of Observed Records",
    tag = "(b)",
    x = expression(gamma),
    y = "Number of Records"
  ) +
  pub_theme +
  theme(
    legend.position = "right",
    plot.tag.position = c(0.02, 0.97),
    plot.tag = element_text(
      family = "Times New Roman",
      size = 11,
      face = "plain"
    )
  )


# p <- p +
#   labs(
#     title = "Distribution of Likelihood-Ratio Test",
#     tag = "(c)",
#     x = expression(gamma),
#     y = "Number of Records"
#   ) +
#   pub_theme +
#   theme(
#     legend.position = "right",
#     plot.tag.position = c(0.02, 0.97),
#     plot.tag = element_text(
#       family = "Times New Roman",
#       size = 11,
#       face = "plain"
#     )
#   )
# 5. Combine and Save ---------------------------------------------------------

# Combine side-by-side with grouped legend alignment
final_figure <- (p_rates | p_records | p) +
  plot_layout(widths = c(1, 1, 1), guides = "collect") &
  theme(legend.position = "bottom",
legend.justification = "center")

print(final_figure)

#ggsave("vignettes/plots_power_error/plot_Power_TypeI_ldm.pdf", plot = final_figure, width = 11, height = 5, units = "in", dpi = 300)
ggsave("vignettes/plots_power_error/plot_Power_TypeI_ldm.png", plot = final_figure, width = 11, height = 5, units = "in", dpi = 300)
