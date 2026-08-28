library(ggplot2)
library(dplyr)
library(tidyr)

# ---------------------------
# 1. Journal theme
# ---------------------------
theme_extremes <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      strip.background = element_blank(),
      strip.text = element_text(size = 11, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_blank(),
      legend.position = "top"
    )
}

# ---------------------------
# 2. Main plotting function
# ---------------------------
plot_extremes <- function(data, param_col) {

  data$theta <- as.numeric(data$theta)
  data$T <- factor(data$T)

  # ---------------------------
  # Plot 1: Estimated vs True
  # ---------------------------
  p1 <- data %>%
    filter(Metric == "AVG_param") %>%
    ggplot(aes(x = theta, y = .data[[param_col]],
               linetype = T, shape = T, group = T)) +
    geom_line(color = "black", size = 0.6) +
    geom_point(color = "black", size = 2) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "gray40") +
    labs(
      x = expression(theta),
      y = bquote(hat(.(param_col)))
    ) +
    theme_extremes()

  # ---------------------------
  # Plot 2: Bias
  # ---------------------------
  p2 <- data %>%
    filter(Metric == "AVG_bias") %>%
    ggplot(aes(x = theta, y = .data[[param_col]],
               linetype = T, group = T)) +
    geom_line(color = "black", size = 0.6) +
    geom_point(color = "black", size = 2) +
    geom_hline(yintercept = 0,
               linetype = "dashed", color = "gray40") +
    labs(
      x = expression(theta),
      y = "Bias"
    ) +
    theme_extremes()

  # ---------------------------
  # Plot 3: Std deviation
  # ---------------------------
  p3 <- data %>%
    filter(Metric %in% c("AVG_Emp_std", "AVG_Theo_std")) %>%
    ggplot(aes(x = theta, y = .data[[param_col]],
               linetype = Metric, group = Metric)) +
    geom_line(color = "black", size = 0.6) +
    facet_wrap(~ T) +
    labs(
      x = expression(theta),
      y = "Standard deviation"
    ) +
    theme_extremes()

  # ---------------------------
  # Plot 4: Multi-parameter (optional)
  # ---------------------------
  p4 <- data %>%
    filter(Metric == "AVG_param") %>%
    pivot_longer(cols = c(thetaHat, scaleHat, shapeHat),
                 names_to = "Parameter",
                 values_to = "Value") %>%
    ggplot(aes(x = theta, y = Value,
               linetype = T, group = T)) +
    geom_line(color = "black", size = 0.6) +
    facet_wrap(~ Parameter, scales = "free_y") +
    labs(
      x = expression(theta),
      y = "Estimate"
    ) +
    theme_extremes()

  return(list(
    est_vs_true = p1,
    bias = p2,
    std_dev = p3,
    multi_param = p4
  ))
}
