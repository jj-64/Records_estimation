# ####################################################
# # Assess validity of Wilks' chi-square approximation
# # Exponential (H0) vs Weibull (H1)
# ####################################################

library(ggplot2)
library(grid)
library(gridExtra)
library(readxl)

###################################################
# INPUT
###################################################

data <- read_excel("data/model_selection/nested/model_selection_exp_2_weibull_1.2_ynm_ynm_T=50.xlsx",
                   sheet = "type_I_error")
Lambda <-data %>% filter (trend_value == 1.2) %>% pull(LR)

df_chi <- 1
alpha  <- 0.05

###################################################
# DIAGNOSTICS
###################################################

ks_res <- ks.test(Lambda, "pchisq", df=df_chi)

crit <- qchisq(1-alpha, df_chi)

emp_size <- mean(Lambda > crit)

diag_text <- paste0(
  "n = ", length(Lambda),
  "\nMean = ", round(mean(Lambda),3),
  "\nVar = ", round(var(Lambda),3),
  "\nKS p = ", signif(ks_res$p.value,3),
  "\nEmpirical size = ", round(emp_size,3),
  "\nNominal size = ", alpha
)

###################################################
# DATA FOR CHI-SQUARE CURVE
###################################################

x <- seq(0, max(max(Lambda), crit)*1.2, length.out = 1000)

chi_df <- data.frame(
  x = x,
  density = dchisq(x, df=df_chi)
)

###################################################
# EXTREME THEME
###################################################

theme_extreme <- theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.grid.major = element_line(
      color = "grey86",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    axis.text = element_text(
      color = "black"
    ),
    axis.title = element_text(
      color = "black",
      face = "bold"
    ),
    plot.title = element_text(
      color = "black",
      face = "bold",
      size = 12,
      hjust = .5
    ),
    plot.subtitle = element_text(
      color = "grey80",
      hjust = .5
    ),
    legend.background = element_rect(
      fill = "white"
    ),
    legend.text = element_text(
      color = "black"
    ),
    legend.title = element_text(
      color = "black"
    )
  )

###################################################
# MAIN PLOT
###################################################

p <- ggplot(data.frame(Lambda), aes(Lambda)) +

  geom_histogram(
    aes(y = after_stat(density),
        fill = "Empirical"),
    bins = 35,
    alpha = .75,
    color = "black"
  ) +

  scale_x_continuous(limits = c(0,10)) +
  scale_y_continuous(limits = c(0,1)) +
  geom_line(
    data = chi_df,
    aes(x, density,
        color = "χ²(1)"),
    linewidth = 1.2
  ) +

  geom_vline(
    xintercept = crit,
    color = "#009E73",
    linewidth = 1.1,
    linetype = "dashed"
  ) +

  annotate(
    "label",
    x = Inf,
    y = Inf,
    label = diag_text,
    hjust = 1.05,
    vjust = 1.05,
    size = 5.2,
    color = "black",
    fill = "white",
    label.size = 0.2
  ) +

  scale_fill_manual(
    values = c("Empirical" = "#0072B2")
  ) +

  scale_color_manual(
    values = c("χ²(1)" = "#E69F00")
  ) +

  labs(
    # title =
    #   expression(
    #     paste(
    #       "Likelihood Ratio Distribution under ",
    #       H[0]
    #     )
    #   ),
    # subtitle =
    #   expression(
    #     paste(
    #       "Empirical Distribution vs. Theoretical ",
    #       chi[1]^2
    #     )
    #   ),
    x = expression(Lambda),
    y = "Density",
    fill = "",
    color = ""
  ) +

  theme_extreme

p
###################################################
# SAVE HIGH-RES FIGURE
###################################################

ggsave(
  "vignettes/plots_wilks_lambda/wilks_exp_weibull_ynm_ynm.png",
  p,
  width = 10,
  height = 7,
  dpi = 600,
  bg = "#0F1117"
)

print(p)
