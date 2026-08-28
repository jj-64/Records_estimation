frechet <- list(
  name = "frechet",
  dparams = 2,  # (A: scale, a: shape)
  logpdf = function(x, params) {
    if (any(x <= 0)) return(-Inf)
    sum(log(VGAM::dfrechet(x, scale = params$scale, shape = params$shape)))
  }
)

norm <- list(
  name = "normal",
  dparams = 2, # (mu, sigma)
  logpdf = function(x, params) {
    sum(log(dnorm(x, mean = params$mean, sd = params$sd)))
  }
)

gumbel <- list(
  name = "gumbel",
  dparams = 2, # (loc, scale)
  logpdf = function(x, params) {
    sum(log(VGAM::dgumbel(x, loc = params$loc, scale = params$scale)))
  }
)

weibull <- list(
  name = "weibull",
  dparams = 2, # (shape, scale)
  logpdf = function(x, params) {
    sum(log(dweibull(x, shape = params$shape, scale = params$scale)))
  }
)

cauchy <- list(
  name = "cauchy",
  dparams = 2, # (loc,  scale)
  logpdf = function(x, params) {
    sum(log(dcauchy(x, location = params$loc, scale = params$scale)))
  }
)

uniform <- list(
  name = "uniform",
  dparams = 2, # (min,  max)
  logpdf = function(x, params) {
    sum(log(dunif(x, min = params$min, max = params$max)))
  }
)
