#' @title Record Rate
#' @description Computes the record rate/Probability ob observing a record at times \eqn{P(R_t)} at time \eqn{t}
#'   (or its asymptotic limit as \eqn{t \to \infty})
#'   different record process.
#'
#' @param model character, model name one of "iid", "LDM", "YNM" or "DTRW"
#' @param t Time index (integer). If \code{t = Inf} or \code{NULL},
#'   the asymptotic record rate \eqn{P(R_\infty)} is returned.
#' @param approximate boolean (Default = FALSE).
#' @param ... additional arguments specific for each model
#' @details
#' For more details look at \code{\link{rec_rate_YNM}}, \code{\link{rec_rate_LDM}},
#' \code{\link{rec_rate_DTRW}}, \code{\link{rec_rate_iid}}
#' @return A numeric value representing the record rate probability.
#'
#' @examples
#' rec_rate(model = "iid", t=10)
#' rec_rate("DTRW",10, approximate = FALSE)
#' rec_rate("DTRW",10, approximate = TRUE)
#' rec_rate("LDM",t = 10,theta = 0.1, scale= 1)
#' rec_rate(model = "YNM", gamma = 1.1, t = 10)
#' rec_rate(model = "YNM", gamma = 1.1, t = Inf)
#' @export
rec_rate <- function(model = c("iid", "DTRW", "YNM", "LDM"), t, approximate = FALSE, ...){

  model <- match.arg(model)
  args <- list(...)

  if(model == "iid"){
    return(rec_rate_iid(t))
  } else if (model == "DTRW"){
    return(rec_rate_DTRW(t, approximate = approximate))
  } else if (model == "YNM"){
    return(rec_rate_YNM(gamma = args$gamma, t= t))
  } else if (model == "LDM"){
    return(rec_rate_LDM(theta = args$theta, t, scale = args$scale))
  }
}

##  IID ------------------------------
#' @title Record Rate in the Classical Model (iid)
#'
#' @description Computes the probability of observing a record at time \eqn{t} for an i.i.d. sequence.
#' @param t Time index (integer).
#' @return Probability of a record occurring at time \eqn{t}, given by:
#' \deqn{ P(R_t) = \frac{1}{t} }
#' #' @examples
#' rec_rate_iid(10)  # Finite-time record probability at t = 10
#' 0.1
rec_rate_iid = function(t){
  1/t
}

##  DTRW--------------------------
#' @title Record Rate for the Discrete-Time Random Walk (DTRW) Process
#'
#' @description Computes the record rate/probability of observing record \eqn{P(R_t)} at time \eqn{t}
#'   (or its asymptotic limit as \eqn{t \to \infty}) under symmetric discrete-time random walk (DTRW) process.
#'   Automatically uses the finite-time exact formula or the asymptotic approximation.
#'
#' @param t Integer or numeric time step (\code{t >= 1}). If \code{t = Inf}, the asymptotic form is used.
#' @param approximate Logical; if \code{TRUE}, the asymptotic form \eqn{1 / \sqrt{\pi t}} is used
#'   regardless of \code{t}. Default is \code{FALSE}.
#'
#' @return Numeric value giving the probability that time \code{t} is a record.
#'
#' @details
#' The **finite-time** record rate is based on the *survival probability* of a symmetric random walk:
#' \deqn{
#' P(R_t) = S(t)
#' }
#' where \eqn{S(t)} is the probability the walk has not crossed zero up to time \eqn{t}.
#'
#' The **asymptotic form** as \eqn{t \to \infty} is:
#' \deqn{
#' P(R_t) \approx \frac{1}{\sqrt{\pi t}}
#' }
#'
#' @seealso
#'   \code{\link{Survival}} for the finite-time survival probability.
#'
#' @examples
#' # Finite-time exact rate
#' rec_rate_DTRW(10)
#'
#' # Asymptotic approximation
#' rec_rate_DTRW(1000, approximate = TRUE)
#'
#' # Vectorized over t
#' t_seq <- 1:100
#' plot(t_seq, sapply(t_seq, rec_rate_DTRW), type = "l", col = "blue", lwd = 2,
#'      ylab = "Record Rate", xlab = "t")
#' lines(t_seq, sapply(t_seq, rec_rate_DTRW, approximate = TRUE),
#'       col = "red", lty = 2)
#' legend("topright", legend = c("Exact", "Asymptotic"), col = c("blue", "red"),
#'        lty = c(1,2), lwd = 2, bty = "n")
#'
#' @export
rec_rate_DTRW <- function(t, approximate = FALSE) {
  if (is.infinite(t) || approximate) {
    # Asymptotic form
    return(1 / sqrt(pi * t))
  } else {
    # Finite-time exact record probability
    return(Survival(t))
  }
}

## LDM ----------------------------
#' @title Record Rate in the Linear Drift Model (LDM)
#'
#' @description Computes the record rate/probability of observing record \eqn{P(R_t)} at time \eqn{t}
#'   (or its asymptotic limit as \eqn{t \to \infty}) under a Linear Drift Model (LDM)
#'   with an underlying Gumbel distribution.
#'
#' @param theta Drift parameter (\eqn{\theta > 0}).
#' @param t Time index (integer). If \code{t = Inf} or \code{NULL},
#'   the asymptotic record rate \eqn{P(R_\infty)} is returned.
#' @param scale Scale parameter of the Gumbel distribution (default = 1).
#' @param loc Optional location parameter (not used in the probability formula, included for consistency).
#' @details
#' For finite time \eqn{t}, the record rate is given by:
#' \deqn{
#' P(R_t) = \frac{1 - e^{-\theta / \text{scale}}}{1 - e^{-\theta t / \text{scale}}}
#' }
#'
#' As \eqn{t \to \infty}, the record rate converges to:
#' \deqn{
#' P(R_\infty) = 1 - e^{-\theta / \text{scale}}
#' }
#'
#' This model corresponds to linearly drifting observations where each new value
#' is drawn from a Gumbel distribution shifted by a constant drift \eqn{\theta}.
#'
#' @return A numeric value representing the record rate probability.
#'
#' @examples
#' # Finite-time record rate
#' rec_rate_LDM(theta = 0.5, t = 10, scale = 1, loc = 0)
#' # [1] 0.3961385
#'
#' # Asymptotic record rate (t -> infinity)
#' rec_rate_LDM(theta = 0.5, t = Inf, scale = 1, loc =0)
#' # [1] 0.3934693
#'
#' # Default behavior returns asymptotic rate
#' rec_rate_LDM(theta = 0.5, scale = 1, loc = 0)
#'  # [1] 0.3934693
#'
#' # Compare convergence
#' t_seq <- 1:100
#' rates <- sapply(t_seq, function(tt) rec_rate_LDM(0.5, tt, scale = 1))
#' plot(t_seq, rates, type = "l", col = "blue", lwd = 2,
#'      ylab = "Record Rate", xlab = "t")
#' abline(h = rec_rate_LDM(0.5, Inf, scale = 1), col = "red", lty = 2)
#'
#' @export
rec_rate_LDM <- function(theta, t = Inf, loc = 0, scale = 1) {
  if (is.infinite(t[1]) || is.null(t)) {
    # Asymptotic case
    return(1 - exp(-theta / scale))
  } else {
    # Finite-time record rate
    return((1 - exp(-theta / scale)) / (1 - exp(-theta * t / scale)))
  }
}



## YNM-Nevzorov -----------------------------
#' @title Record Rate in YNM-Nevzorov Model
#' @description Computes the record rate/Probability ob observing a record at times \eqn{P(R_t)} at time \eqn{t}
#'   (or its asymptotic limit as \eqn{t \to \infty})
#'   in the Yang–Nevzorov–Margolin (YNM) record process.
#'
#' @param gamma The memory parameter (\eqn{\gamma > 1}).
#' @param t Time index (integer). If \code{t = Inf} or \code{NULL},
#'   the asymptotic record rate \eqn{P(R_\infty)} is returned.
#'
#' @details
#' For finite \eqn{t}, the record rate is:
#' \deqn{ P(R_t) = \frac{\gamma^t}{\gamma} \cdot \frac{1 - \gamma}{1 - \gamma^t} }
#'
#' As \eqn{t \to \infty}, it converges to:
#' \deqn{ P(R_\infty) = 1 - \frac{1}{\gamma} }
#'
#' @return A numeric value representing the record rate probability.
#'
#' @examples
#' rec_rate_YNM(gamma = 1.1, t = 10)
#' rec_rate_YNM(gamma = 1.4, t = Inf)
#' rec_rate_YNM(gamma = 1.2)  # defaults to asymptotic
#'
#' @export
rec_rate_YNM <- function(gamma, t = Inf) {
  if (is.infinite(t[1]) || is.null(t)) {
    return(1 - 1/gamma)
  } else {
    return((gamma^t - gamma^(t-1)) / (gamma^t -1 ) )
  }
}

