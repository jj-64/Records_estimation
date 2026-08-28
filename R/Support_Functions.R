# helper for defaults
`%||%` <- function(a, b) if (!is.null(a)) a else b

########################################################################
##To compute the bounds/ CI intervals
#' bounds: Confidence-Interval
#'
#'Compute the confidence interval bounds
#' @param value the value or the mean of which we are computing the CI
#' @param z number for z-score from normal distribution
#' @param variance number for the variance (not standard deviation)
#'
#' @returns vector of two: lower and upper bounds of C.I.
#' @export
#'
#' @examples bounds(5, z=qnorm(1-0.05/2), variance = 2)
#' # 2.228192  7.771808
bounds = function(value,z,variance){
  upper = value+ z*sqrt(variance)
  lower = value- z*sqrt(variance)
  return(c(lower,upper))
}

############# Install packages ###################
## Install missing packages
#'
#' install missing packages if needed
#' @param pkgs character vector, names of the packages
#' @export
#'
#' @examples
#' # required_packages <- c(
#' #   "tidyverse",
#' #   "moments"
#' # )
#' # install_if_missing(required_packages)
install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!suppressWarnings(require(p, character.only = TRUE))) {
      install.packages(p, dependencies = TRUE)
      library(p, character.only = TRUE)
    }
  }
}

##################### Load packages #################
## Load needed packages
#'
#' load missing packages if needed
#' @param required_packages character vector, names of the packages
#' @export
#'
#' @examples
#' # required_packages <- c(
#' #  "tidyverse",
#' # "moments"
#' # )
#' # load_package(required_packages)
load_package = function(required_packages){
lapply(required_packages, require, character.only = TRUE)
}

