library(ggplot2)
library(dplyr)
library(tidyr)

# ---------------------------
# 1. Palette + Theme
# ---------------------------
extreme_colors <- c("#2C3E50", "#1F7A8C", "#BF8B2E", "#9E2A2B")

theme_extremes_color <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "gray90", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      strip.background = element_rect(fill = "gray95", color = NA),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "top",
      legend.title = element_blank()
    )
}

# ---------------------------
# 2. Generic plotting function
# ---------------------------
plot_extremes_color <- function(data, obs_col, param_cols = NULL) {

  data[[obs_col]] <- as.numeric(data[[obs_col]])
  data$T <- factor(data$T)

  # Detect parameter columns automatically
  if (is.null(param_cols)) {
    param_cols <- setdiff(
      colnames(data),
      c("T", obs_col, "Metric")
    )
  }

  # ---------------------------
  # Plot 1: Estimated vs True
  # ---------------------------
  p1_list <- lapply(param_cols, function(param) {

    p <- ggplot(
      data %>% filter(Metric == "AVG_param"),
      aes(x = .data[[obs_col]],
          y = .data[[param]],
          color = T, linetype = T, group = T)
    ) +
      geom_line(size = 0.8) +
      geom_point(size = 2) +
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "gray40") +
      scale_color_manual(values = extreme_colors) +
      labs(
        x = obs_col,
        y = paste("Estimate of", param)
      ) +
      ggtitle(param) +
      theme_extremes_color()

    return(p)
  })

  names(p1_list) <- param_cols

  # ---------------------------
  # Plot 2: Bias
  # ---------------------------
  p2_list <- lapply(param_cols, function(param) {

    ggplot(
      data %>% filter(Metric == "AVG_bias"),
      aes(x = .data[[obs_col]],
          y = .data[[param]],
          color = T, linetype = T, group = T)
    ) +
      geom_line(size = 0.8) +
      geom_point(size = 2) +
      geom_hline(yintercept = 0,
                 linetype = "dashed", color = "gray40") +
      scale_color_manual(values = extreme_colors) +
      labs(
        x = obs_col,
        y = paste("Bias of", param)
      ) +
      ggtitle(param) +
      theme_extremes_color()
  })

  names(p2_list) <- param_cols

  # ---------------------------
  # Plot 3: Std deviation
  # ---------------------------
  p3_list <- lapply(param_cols, function(param) {

    ggplot(
      data %>% filter(Metric %in% c("AVG_Emp_std", "AVG_Theo_std")),
      aes(x = .data[[obs_col]],
          y = .data[[param]],
          color = Metric,
          linetype = Metric,
          group = Metric)
    ) +
      geom_line(size = 0.8) +
      facet_wrap(~ T) +
      scale_color_manual(values = c("#2C3E50", "#BF8B2E")) +
      labs(
        x = obs_col,
        y = paste("Std Dev of", param)
      ) +
      ggtitle(param) +
      theme_extremes_color()
  })

  names(p3_list) <- param_cols

  # ---------------------------
  # Plot 4: Multi-parameter
  # ---------------------------
  p4 <- data %>%
    filter(Metric == "AVG_param") %>%
    pivot_longer(cols = all_of(param_cols),
                 names_to = "Parameter",
                 values_to = "Value") %>%
    ggplot(aes(x = .data[[obs_col]],
               y = Value,
               color = T, group = T)) +
    geom_line(size = 0.8) +
    facet_wrap(~ Parameter, scales = "free_y") +
    scale_color_manual(values = extreme_colors) +
    labs(
      x = obs_col,
      y = "Estimate"
    ) +
    theme_extremes_color()

  # ---------------------------
  # Plot 5: Coverage probability (NEW 👍)
  # ---------------------------
  if ("Coverage_proba" %in% data$Metric) {

    p5 <- data %>%
      filter(Metric == "Coverage_proba") %>%
      pivot_longer(cols = all_of(param_cols),
                   names_to = "Parameter",
                   values_to = "Coverage") %>%
      ggplot(aes(x = .data[[obs_col]],
                 y = Coverage,
                 color = T, group = T)) +
      geom_line(size = 0.8) +
      facet_wrap(~ Parameter) +
      geom_hline(yintercept = 0.95,
                 linetype = "dashed",
                 color = "gray40") +
      scale_color_manual(values = extreme_colors) +
      labs(
        x = obs_col,
        y = "Coverage probability"
      ) +
      theme_extremes_color()

  } else {
    p5 <- NULL
  }

  return(list(
    est_vs_true = p1_list,
    bias = p2_list,
    std_dev = p3_list,
    multi_param = p4,
    coverage = p5
  ))
}

# library(ggplot2)
# library(dplyr)
# library(tidyr)
#
# # ---------------------------
# # 1. Color palette
# # ---------------------------
# extreme_colors <- c(
#   "#2C3E50",  # deep blue-gray
#   "#1F7A8C",  # muted teal
#   "#BF8B2E",  # soft amber
#   "#9E2A2B"   # muted red
# )
#
# # ---------------------------
# # 2. Theme (color version)
# # ---------------------------
# theme_extremes_color <- function() {
#   theme_bw(base_size = 12) +
#     theme(
#       panel.grid.major = element_line(color = "gray90", size = 0.3),
#       panel.grid.minor = element_blank(),
#       panel.border = element_rect(color = "black", fill = NA),
#       strip.background = element_rect(fill = "gray95", color = NA),
#       strip.text = element_text(size = 11, face = "bold"),
#       axis.title = element_text(size = 12),
#       axis.text = element_text(size = 10),
#       legend.position = "top",
#       legend.title = element_blank()
#     )
# }
#
# # ---------------------------
# # 3. Main function
# # ---------------------------
# plot_extremes_color <- function(data, param_col) {
#
#   data$theta <- as.numeric(data$theta)
#   data$T <- factor(data$T)
#
#   # ---------------------------
#   # Plot 1: Estimated vs True
#   # ---------------------------
#   p1 <- data %>%
#     filter(Metric == "AVG_param") %>%
#     ggplot(aes(x = theta, y = .data[[param_col]],
#                color = T, linetype = T, group = T)) +
#     geom_line(size = 0.8) +
#     geom_point(size = 2) +
#     geom_abline(slope = 1, intercept = 0,
#                 linetype = "dashed", color = "gray40") +
#     scale_color_manual(values = extreme_colors) +
#     labs(
#       x = expression(theta),
#       y = bquote(hat(.(param_col)))
#     ) +
#     theme_extremes_color()
#
#   # ---------------------------
#   # Plot 2: Bias
#   # ---------------------------
#   p2 <- data %>%
#     filter(Metric == "AVG_bias") %>%
#     ggplot(aes(x = theta, y = .data[[param_col]],
#                color = T, linetype = T, group = T)) +
#     geom_line(size = 0.8) +
#     geom_point(size = 2) +
#     geom_hline(yintercept = 0,
#                linetype = "dashed", color = "gray40") +
#     scale_color_manual(values = extreme_colors) +
#     labs(
#       x = expression(theta),
#       y = "Bias"
#     ) +
#     theme_extremes_color()
#
#   # ---------------------------
#   # Plot 3: Std deviation
#   # ---------------------------
#   p3 <- data %>%
#     filter(Metric %in% c("AVG_Emp_std", "AVG_Theo_std")) %>%
#     ggplot(aes(x = theta, y = .data[[param_col]],
#                color = Metric, linetype = Metric,
#                group = Metric)) +
#     geom_line(size = 0.8) +
#     facet_wrap(~ T) +
#     scale_color_manual(values = c("#2C3E50", "#BF8B2E")) +
#     labs(
#       x = expression(theta),
#       y = "Standard deviation"
#     ) +
#     theme_extremes_color()
#
#   # ---------------------------
#   # Plot 4: Multi-parameter
#   # ---------------------------
#   p4 <- data %>%
#     filter(Metric == "AVG_param") %>%
#     pivot_longer(cols = c(thetaHat, scaleHat, shapeHat),
#                  names_to = "Parameter",
#                  values_to = "Value") %>%
#     ggplot(aes(x = theta, y = Value,
#                color = T, group = T)) +
#     geom_line(size = 0.8) +
#     facet_wrap(~ Parameter, scales = "free_y") +
#     scale_color_manual(values = extreme_colors) +
#     labs(
#       x = expression(theta),
#       y = "Estimate"
#     ) +
#     theme_extremes_color()
#
#   return(list(
#     est_vs_true = p1,
#     bias = p2,
#     std_dev = p3,
#     multi_param = p4
#   ))
# }
