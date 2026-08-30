
library(readxl)
sim_data <- read_excel("data/param_est/dtrw/norm/Simulation_dtrw_Xt_Rn_overT_norm_sd=1.xlsx",
                                                    sheet = "All_Results_record")

sim_data <- read_excel("data/param_est/ldm/gumbel/Simulation_ldm_Xt_Rn_gumbel_location=0_scale=1.xlsx",
                                                    sheet = "All_Results_record")
head(sim_data)

agg_T <- aggregate_by_T(
  data = sim_data,
  true_values = list(location = 0, scale = 1),
  param_names = c("location", "scale", "theta")
)

agg_NT <- aggregate_by_NT(
  data = sim_data,
  true_values = list(location = 0, scale = 1),
  param_names = c("location", "scale")
)

head(agg_T)
head(agg_NT)
