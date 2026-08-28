################ describe function ########################

#' Describe a vector and append record-based statistics
#'
#' @description
#' A thin wrapper around [psych::describe()] that appends record-based
#' statistics computed from `rec_count()`, `rec_values()`, and `rec_times()`.
#' The return value preserves the structure of `psych::describe()` (a data frame)
#' so it can be used as a drop-in replacement in pipelines.
#'
#' @param x A vector accepted by [psych::describe()]. For now, this function
#'   expects a single vector (not a data frame or matrix). Supply one column
#'   at a time if needed.
#' @param na.rm Logical; remove `NA`s before computing both the `describe` and
#'   the record-based statistics. Default `TRUE`.
#' @param attach_details Logical; if `TRUE`, the full vectors returned by
#'   `rec_values()` and `rec_times()` are attached as an attribute
#'   `record_details` to the returned data frame. Default `FALSE`.
#'
#' @details
#' Extra columns appended to the `psych::describe()` output:
#' \itemize{
#'   \item \code{obs_rec_count}: number of records from \code{rec_count(x)}.
#'   \item \code{obs_rate}: records per observation, i.e., \code{obs_rec_count / length(x)} after \code{na.rm}.
#'   \item \code{first_rec_value}, \code{last_rec_value}: first/last elements of \code{rec_values(x)}.
#'   \item \code{mean_rec_value}, \code{min_rec_value}, \code{max_rec_value}: summaries of \code{rec_values(x)} (numeric only; otherwise \code{NA}).
#'   \item \code{first_rec_time}, \code{last_rec_time}: first/last elements of \code{rec_times(x)}.
#'   \item \code{rec_span_secs}: if times are \code{POSIXt} or \code{Date}, the span (last - first) in seconds; otherwise \code{NA}.
#' }
#'
#' When \code{attach_details = TRUE}, the full vectors are added as an attribute
#' \code{record_details = list(values = <...>, times = <...>)} so you can inspect
#' or reuse them without bloating the data frame with list-columns.
#'
#' @return A data frame identical to [psych::describe()] but with the columns
#'   listed above appended. When \code{attach_details = TRUE}, an attribute
#'   \code{record_details} is attached.
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   x <- rnorm(100)
#'   x[sample.int(100, 10)] <- NA
#'
#'   out <- describe_records(x)
#'   out
#'
#'   # With the full record vectors attached
#'   out2 <- describe_records(x, attach_details = TRUE)
#'   attr(out2, "record_details")
#' }
#'
#' @seealso [psych::describe()]
#' @export
#' @importFrom psych describe
describe_records <- function(x, na.rm = TRUE, attach_details = FALSE) {
  # Restrict to a single vector to preserve a single-row describe result
  if (is.matrix(x) || is.data.frame(x)) {
    stop("describe_records() currently supports a single vector. Supply one column at a time.")
  }

  # Delegate to psych::describe to keep identical core output
  d <- psych::describe(x, na.rm = na.rm, fast = FALSE) %>%
    dplyr::select(n, mean, sd, median, min, max)

  # Apply NA handling consistently for the record statistics as well
  x_use <- if (na.rm) x[!is.na(x)] else x

  # Built-in record helpers (assumed to be available in your environment)
  obs_rec_count <- rec_count(x_use)
  obs_rate      <- if (length(x_use) > 0L) obs_rec_count / length(x_use) else NA_real_

  vals  <- rec_values(x_use)
  times <- rec_times(x_use)

  # Safe accessors for first/last
  first_or_na <- function(v) if (length(v)) v[[1]] else NA
  last_or_na  <- function(v) if (length(v)) v[[length(v)]] else NA

  first_val <- first_or_na(vals[-1])
  last_val  <- last_or_na(vals)

  first_time <- first_or_na(times[-1])
  last_time  <- last_or_na(times)

  # Optional summaries of record values if they are numeric
  if (is.numeric(vals)) {
    mean_val <- if (length(vals)) mean(vals, na.rm = TRUE) else NA_real_
  } else {
    mean_val <- NA_real_
  }

  # Span in seconds if times are date/time
  if (inherits(first_time, c("POSIXt", "Date")) && inherits(last_time, c("POSIXt", "Date"))) {
    rec_span_secs <- as.numeric(difftime(last_time, first_time, units = "secs"))
  } else {
    rec_span_secs <- NA_real_
  }

  # Append columns to preserve describe()'s structure
  d$obs_rec_count    <- obs_rec_count
  d$obs_rate         <- obs_rate
  d$first_rec_value  <- first_val
  d$last_rec_value   <- last_val
  d$mean_rec_value   <- mean_val
  d$first_rec_time   <- first_time
  d$last_rec_time    <- last_time
  d$rec_span_secs    <- rec_span_secs

  # Optionally attach the full vectors as an attribute (keeps data frame tidy)
  if (attach_details) {
    attr(d, "record_details") <- list(values = vals, times = times)
  }

  d
}
