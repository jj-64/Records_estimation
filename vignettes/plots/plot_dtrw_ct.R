library(readxl)
library(tidyverse)
library(patchwork)
library(scales)
library(ggplot2)
library(dplyr)

#___________________________________
# Read data ---------
#___________________________________
true_param = 1

df_unif <- read_excel(paste0("data/param_est/dtrw/unif/Simulation_dtrw_Xt_Rn_overT_unif_max=",true_param,"_BiasCorrected.xlsx"), sheet = "All Par_record")
df_unif_biased <- read_excel(paste0("data/param_est/dtrw/unif/Simulation_dtrw_Xt_Rn_overT_unif_max=",true_param,"_biased.xlsx"), sheet = "All Par_record")
df_norm <- read_excel(paste0("data/param_est/dtrw/norm/Simulation_dtrw_Xt_Rn_overT_norm_sd=",true_param,".xlsx"), sheet = "All Par_record")

colnames(df_unif)[1:2] = c("T", "AVG_param")
colnames(df_unif_biased)[1:2] = c("T", "AVG_param")
colnames(df_norm)[1:2] = c("T", "AVG_param")

#___________________________________
# Extract numeric T ----------
#___________________________________

df_plot <- bind_rows(
  df_unif[,1:2] %>% mutate(Type = "Uniform"),
  df_unif_biased[,1:2] %>% mutate(Type = "Uniform (biased)"),
  df_norm[,1:2] %>% mutate(Type = "Normal")
) %>%
  mutate(T_num = as.numeric(gsub("T=", "", T)),
         cT = AVG_param/true_param)

#___________________________________
# 1. Bias plot --------
#___________________________________

# Bias plot
cT_plot <- ggplot(df_plot,
                    aes(x = T_num,
                        y = cT,
                        color = Type,
                        linetype = Type,
                        shape = Type,
                        group = Type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Normal" = "#0072B2",
    "Uniform" = "#009E73",
    "Uniform (biased)" = "#E69F00"
  )) +
  scale_linetype_manual(values = c(
    "Normal" = "solid",
    "Uniform" = "dashed",
    "Uniform (biased)" = "longdash"
  )) +
  scale_shape_manual(values = c(
    "Normal" = 16,
    "Uniform" = 17,
    "Uniform (biased)" = 18
  )) +
  scale_y_continuous(n.breaks=10) +
  theme_bw() +
  labs(
    title = "Average Bias across T",
    x = "series length (T)",
    y = "Bias",
    color = "",
    linetype = "",
    shape = ""
  )

cT_plot


#___________________________________
# Save all figures --------
#___________________________________

ggsave(paste0("vignettes/plots_output/dtrw/cT_plot_param=",true_param,".png"),
       cT_plot,
       width = 10,
       height = 6,
       dpi = 300)
