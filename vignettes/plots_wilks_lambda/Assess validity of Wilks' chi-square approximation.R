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

data <- read_excel("data/model_selection/nested/model_selection_exp_2_weibull_1.2_ldm_ldm_T=100.xlsx",
                                                                  sheet = "type_I_error")
Lambda <-data %>% filter (trend_value == 0.1) %>% pull(LR)

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
  "vignettes/plots_wilks_lambda/Wilks_Diagnostic_Extreme.png",
  p,
  width = 10,
  height = 7,
  dpi = 600,
  bg = "#0F1117"
)

print(p)


# # Read LR statistics
# dat <- read.csv("Lambda_H0.csv")
#
# Lambda <- dat$LR
#
# # Degrees of freedom
# df <- 1
#
# ####################################################
# # 1. Visual comparison
# ####################################################
#
# hist(Lambda,
#      probability = TRUE,
#      breaks = 30,
#      col = "lightgray",
#      border = "white",
#      main = expression(paste("Empirical distribution of ", Lambda)),
#      xlab = expression(Lambda),
#      xlim = c(0,20))
#
# curve(dchisq(x, df = df),
#       add = TRUE,
#       col = "red",
#       lwd = 3)
#
# legend("topright",
#        legend = c("Chi-square(1)"),
#        col = "red",
#        lwd = 3)
#
# ####################################################
# # 2. QQ plot
# ####################################################
#
# n <- length(Lambda)
#
# theoretical <- qchisq(ppoints(n), df = df)
# empirical   <- sort(Lambda)
#
# qqplot(theoretical,
#        empirical,
#        main = "QQ Plot: Empirical vs Chi-square",
#        xlab = expression(chi[1]^2~quantiles),
#        ylab = expression(Lambda~quantiles))
#
# abline(0, 1, col = "red", lwd = 2)
#
# ####################################################
# # 3. Kolmogorov-Smirnov test
# ####################################################
#
# ks <- ks.test(Lambda,
#               "pchisq",
#               df = df)
#
# print(ks)
#
# ####################################################
# # 4. Type-I error control
# ####################################################
#
# alpha_levels <- c(0.10, 0.05, 0.01)
#
# result <- data.frame(
#   alpha = alpha_levels,
#   theoretical_critical = NA,
#   empirical_rejection_rate = NA
# )
#
# for(i in seq_along(alpha_levels)) {
#
#   alpha <- alpha_levels[i]
#
#   crit <- qchisq(1 - alpha, df = df)
#
#   rej <- mean(Lambda > crit)
#
#   result$theoretical_critical[i] <- crit
#   result$empirical_rejection_rate[i] <- rej
# }
#
# print(result)
#
# ####################################################
# # 5. Confidence intervals for rejection rates
# ####################################################
#
# for(i in 1:nrow(result)) {
#
#   alpha <- result$alpha[i]
#
#   crit <- result$theoretical_critical[i]
#
#   x <- sum(Lambda > crit)
#
#   n <- length(Lambda)
#
#   ci <- prop.test(x, n)$conf.int
#
#   cat("\n")
#   cat("Alpha =", alpha, "\n")
#   cat("Observed rejection rate =", round(x/n,4), "\n")
#   cat("95% CI = [",
#       round(ci[1],4), ",",
#       round(ci[2],4), "]\n")
# }
#
