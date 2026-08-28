library(stats)

logLik_fun_rec <- loglik_registry[["LDM"]][["records"]][["frechet_inv_scale"]]

# parameters
trend =0.1
T_val = 50
shape_obs = 2
scale_obs = 1

## generate LDM series of frechet underlying distr
y= LDM_series(T = T_val,theta = trend, dist = "frechet",
              shape = shape_obs, scale=scale_obs)
# y= c(0.924219871299229,2.20326092425526,0.870623759490709,1.04964455250712,
# 5.12134923643906,2.03810714442534,1.62031496550323,2.62095242113047,
# 3.11796315373072,1.59611497057797,2.53617457658537,5.60308193575109,3.08680613719124,
# 3.14959630464053,2.58470422042793,2.21190347307693,4.2031266612816,7.02506088795092,
# 4.3362425360437,3.24210815322174,4.1290599255864,3.594343133097,3.11128381017565,
# 3.03490749813946,4.73299643271331,3.66051446529391,3.71563557421047,4.42844240689788,
# 3.74193440302009,3.92227415545899,3.64998185409224,4.99899091160745,4.84432238563662,
# 5.03761890155542,4.41567963370727,4.63502271851145,5.26355712021795,4.68760883620642,
# 5.63546684460518,4.85013197421985,6.33503150590016,4.93213862832236,5.37663598946662,
# 5.76630223954827,4.96849122527482,8.02052681621,7.4724379976162,7.08886487466141,
# 6.0144882701532,6.33477597734309)

R = rec_values(y)
L = rec_times(y)

# R = c(0.924219871299229,2.20326092425526,5.12134923643906,5.60308193575109,7.02506088795092,8.02052681621)
# L = c(1,2,5,12,18,46)

data_rec = list(rec_values = R, rec_times = L, time = T_val)
data_rec = data.frame(rec_values = R, rec_times = L, time = T_val)

# data_rec
# rec_values rec_times time
# 1  0.9242199         1   50
# 2  2.2032609         2   50
# 3  5.1213492         5   50
# 4  5.6030819        12   50
# 5  7.0250609        18   50
# 6  8.0205268        46   50

## estimate Max LogLik
MLE_C_record = estimate_model_mle(logLik_fun_rec, obs_type = "records", data = data_rec,
                                  lower_bounds=c(theta=0.001, shape = 0.01, scale=0.01),
                                  upper_bounds = c(theta= min(R/L)-1e-10, shape = 10, scale=10),
                                  start_values = c(
                                    theta = 0.01,
                                    shape = 0.01,
                                    scale =0.01
                                  )
)

MLE_C_record = c(MLE_C_record$par, "objective" = - MLE_C_record$objective ) ## return shape, shape, scale, objective value
MLE_C_record

# MLE_C_record
# theta        shape        scale    objective
# 0.09808864   2.34482291   0.86796437 -20.43898085

H = compute_hessian(loglik_fun = logLik_fun_rec, data_rec, par_vec = MLE_C_record[1:3], eps = 1e-2 )

# H
# $hessian
#        theta     shape      scale
# theta -2109.79540 58.760070 203.596795
# shape    58.76007 -5.420841  -8.192845
# scale   203.59679 -8.192845 -43.789278
#
# $vcov
#          theta        shape        scale
# theta 0.0009493848  0.005046731  0.003469905
# shape 0.0050467305  0.284030063 -0.029676584
# scale 0.0034699047 -0.029676584  0.044522248
#
# $corr
#         theta      shape      scale
# theta 1.0000000  0.3073313  0.5337131
# shape 0.3073313  1.0000000 -0.2639023
# scale 0.5337131 -0.2639023  1.0000000

theta.hat = MLE_C_record["theta"]
scale.hat =  MLE_C_record["scale"]
shape.hat =  MLE_C_record["shape"]

sd.theta = sqrt(diag(H$vcov)[["theta"]])
sd.scale = sqrt(diag(H$vcov)[["scale"]])
sd.shape = sqrt(diag(H$vcov)[["shape"]])

## Log lik of theta vs scale and shape --------
shape_grid <- seq(0.1, 5, length.out = 100)
scale_grid <- seq(0.1, 5, length.out = 100)

## Matrix to store log-likelihood values
LL <- matrix(NA, nrow = length(shape_grid),
             ncol = length(scale_grid))


## Evaluate log-likelihood
for(i in seq_along(shape_grid)) {
  for(j in seq_along(scale_grid)) {

    params <- c(
      theta    = trend,
      shape    = shape_grid[i],
      scale    = scale_grid[j],
      location = 0
    )

    LL[i, j] <- logLik_fun_rec(data = data_rec, params = params)
  }
}

colnames(LL) = scale_grid
rownames(LL) = shape_grid

## Remove infinities if present
LL[!is.finite(LL)] <- NA

## plot Relative contour
LL_rel <- LL - max(LL, na.rm = TRUE)

contour(
  x = scale_grid,
  y = shape_grid,
  z = LL_rel,
  levels = c(-60,-50,-20,-10,-5,-2,-1,-0.5),
  xlab = "Scale",
  ylab = "Shape",
  main = expression(paste("Relative log-likelihood, ", theta," = 0.2"))
)

## Diagnostic 1: Profile likelihood for shape -----
# For the reviewer, the most important figure remains
# ℓp(α)=max⁡θ,σℓ(θ,α,σ).\ell_p(\alpha) = \max_{\theta,\sigma} \ell(\theta,\alpha,\sigma).ℓp​(α)=θ,σmax​ℓ(θ,α,σ).
# Because of the correlation, you should profile over BOTH nuisance parameters.
# For a fixed dataset:
#   α↦ℓ(α,β^(α))\alpha \mapsto \ell(\alpha,\hat\beta(\alpha))α↦ℓ(α,β^​(α))
# where all remaining parameters are optimized.
# A flat profile visually explains everything.
# If the curve is wide and almost horizontal near the maximum, you have exactly the figure the reviewer wants.

# version 1: not very well
# shape.grid <- seq(0.5, 3, length.out = 100)
#
# prof_shape <- numeric(length(shape.grid))
#
# for(i in seq_along(shape.grid)){
#   sh <- shape.grid[i]
#
#   opt <- optim(
#     par = c(theta = theta.hat, scale = scale.hat),
#     fn = function(par)
#     {
#       p <- list(theta = par[1],
#              shape = sh,
#              scale = par[2])
#
#       -logLik_fun_rec(data_rec, p)
#     },
#     method = "L-BFGS-B",
#     lower = c(0.01, 1e-6),
#     upper = c(1, 10)
#   )
#
#   prof_shape[i] <- -opt$value
# }
#
# plot(
#   shape.grid,
#   prof_shape - max(prof_shape),
#   type = "l",
#   lwd = 3,
#   xlab = expression(alpha),
#   ylab = "Profile log-likelihood"
#   #title = "profile likelihood for the Fréchet shape parameter showing weak curvature."
# )
#
# abline(h = -1.92, lty = 2, col = 2) ## for confidence interval

## better version for publication
library(ggplot2)

shape.grid <- c(seq(
  max(0.2, shape.hat*0.4),
  shape.hat*2,
  length.out = 120
), shape.hat)

profLL <- numeric(length(shape.grid))

for(i in seq_along(shape.grid))
{
  sh <- shape.grid[i]

  fit <- optim(
    par = c(theta.hat, scale.hat),
    fn = function(par)
    {
      p <- list(
        theta = par[1],
        shape = sh,
        scale = par[2]
      )

      -logLik_fun_rec(data_rec, p)
    },
    method = "L-BFGS-B",
    lower = c(0.01, 1e-6),
    upper = c(1, 10)
  )

  profLL[i] <- -fit$value
}

profile.df <- data.frame(
  shape = shape.grid,
  logLik = profLL
)

profile.df$relativeLL <-
  profile.df$logLik - max(profile.df$logLik)

profile = ggplot(profile.df,
       aes(shape, relativeLL)) +
  geom_line(linewidth = 1.3,
            colour = "#2C7BB6") +
  geom_hline(yintercept = -1.92,
             linetype = 2,
             colour = "red") +
  geom_hline(yintercept = 0,
             linetype = 1,
             colour = "black") +
  geom_vline(xintercept = shape.hat,
             colour = "black",
             linewidth = 0.8) +
  theme_minimal(base_size = 15) +
  labs(
    x = expression(alpha),
    y = expression(ell[alpha] - max[ell[p]]),
    title = "" #"Profile likelihood for the shape parameter"
  )

profile

cutoff <- max(profLL) - 1.92

CI_shape <- range(
  shape.grid[
    profLL >= cutoff
  ]
)

CI_shape

# Diagnostic 2: Likelihood contours ------
# This is actually more informative for your model.
# The reviewer asked:
#   likelihood profiles or curvature diagnostics
# A contour plot is a curvature diagnostic.
# Fix shape at MLE:
#   α^\hat{\alpha}α^
#   and evaluate
# ℓ(θ,α^,σ).\ell(\theta,\hat{\alpha},\sigma).ℓ(θ,α^,σ).
# (an elongated ellipse)
# then:
#   theta and scale are highly correlated;
# information matrix nearly singular;
# shape becomes unstable.
# This provides a direct explanation.

# theta.seq <- seq(MLE_C_record["theta"] - 0.5,
#                  MLE_C_record["theta"] + 0.5,
#                  length = 80)
#
# scale.seq <- seq(MLE_C_record["scale"]*0.5,
#                  MLE_C_record["scale"]*1.5,
#                  length = 80)
#
# LL <- matrix(NA,
#              nrow = length(theta.seq),
#              ncol = length(scale.seq))
#
# for(i in seq_along(theta.seq))
# {
#   for(j in seq_along(scale.seq))
#   {
#     pars <- list(
#       theta = theta.seq[i],
#       shape = MLE_C_record["shape"],
#       scale = scale.seq[j]
#     )
#
#     LL[i,j] <- logLik_fun_rec(data_rec, pars)
#   }
# }
#
# contour(
#   theta.seq,
#   scale.seq,
#   LL,
#   xlab = expression(theta),
#   ylab = "scale",
#   nlevels = 20
# )
#
# points(MLE_C_record["theta"],
#        MLE_C_record["scale"],
#        pch = 19,
#        col = "red")
#
#
# library(ggplot2)

## better visualization
theta.seq <- seq(
  theta.hat - 2*sd.theta,
  theta.hat + 2*sd.theta,
  length.out = 120
)

scale.seq <- seq(
  scale.hat*0.5,
  scale.hat*1.5,
  length.out = 120
)

grid <- expand.grid(
  theta = theta.seq,
  scale = scale.seq
)

grid$LL <- mapply(
  function(th, sc)
  {
    logLik_fun_rec(
      data_rec,
      list(theta = th,
        shape = shape.hat,
        scale = sc)
    )
  },
  grid$theta,
  grid$scale
)

grid$LLrel <- grid$LL - max(grid$LL)

contour = ggplot(grid, aes(theta, scale, z = LLrel)) +

  geom_contour_filled(
    breaks = c(-20,-10,-6,-4,-2,-1,-0.5,0)
  ) +

  geom_contour(
    colour = "white",
    linewidth = 0.4
  ) +

  geom_point(
    aes(theta.hat,
        scale.hat),
    inherit.aes = FALSE,
    size = 3,
    colour = "red"
  ) +

  theme_minimal(base_size = 15) +

  labs(
    title = "",
      # expression(
      #   paste("Likelihood contours: ",
      #         theta," vs scale")
      # ),
    x = expression(theta),
    y = expression(sigma),
  ) +
  # 4. Dynamic Text Annotation next to the red dot
  geom_text(
    aes(
      x = theta.hat,
      y = scale.hat,
      # Creates a clean label: "θ̂ = [val], σ̂ = [val]"
      label = paste0("hat(theta) == ", round(theta.hat, 2), "*\", \"*hat(sigma) = ", round(scale.hat, 2))
    ),
    inherit.aes = FALSE,
    parse = TRUE,          # Allows R to render mathematical symbols (hats and Greek letters)
    hjust = -0.15,         # Pushes the text slightly to the right of the red dot
    vjust = 0.5,           # Vertically centers the text with the dot
    colour = "white",      # Contrasts well against darker contour backgrounds
    fontface = "bold",
    size = 4.5
  )

contour

## combined plot -------
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


contour <- contour +
  labs(
    title = "Likelihood contours",
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

profile <- profile +
  labs(
    title = "Profile likelihood",
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

combined_plot = ( profile |contour) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  )

combined_plot

ggsave("vignettes/plots_profile_loglik/loglik_geometry_2.png",
       combined_plot,
       width = 6,
       height = 4,
       dpi = 300)

# Diagnostic 3: Hessian and condition number --------
# This is probably the cleanest quantitative diagnostic.
# You can report:
# The observed information matrix was ill-conditioned (condition number = XX), indicating weak identifiability of some parameter combinations.
# A large condition number (>100, >1000) is evidence of a likelihood ridge.

library(numDeriv)

curvature <- H$hessian

info <- -curvature
info

H$vcov

eig <- eigen(curvature)
eig

cond_number <- max(eig$values) /
  min(eig$values)

cond_number

# # 4. Curvature diagnostic -----
# A more quantitative diagnostic is
# I(α^)=−ℓ′′(α^).I(\hat\alpha) = -\ell''(\hat\alpha).I(α^)=−ℓ′′(α^).
# Small curvature means low information.
# "The observed information for α was small (I(α̂)=0.12), explaining the large variance observed in the Monte Carlo experiments."

## theta and scale vs shape --------

