## To add a new model
## register_record_model("new_model", rec_count_mean_new, rec_count_var_new)

#' Registry for record models
#'
#' Internal registry that stores model handlers. Use \code{register_record_model} to add models.
#' @keywords internal
#' @noRd
record_registry <- new.env(parent = emptyenv())
#' Register a record model
#'
#' Register functions for mean, variance and distribution for a named model.
#'
#' @param model Character name.
#' @param mean_fun Function to compute expected number of records.
#' @param var_fun Function to compute variance of record count.
#' @param dist_fun Function to compute probability mass of having m records.
#' @param required_args Character vector of required extra arguments (e.g. 'theta' or 'gamma').
#' @export
register_record_model <- function(model, mean_fun, var_fun, dist_fun, required_args = NULL) {
  record_registry[[model]] <- list(
    mean_fun = mean_fun,
    var_fun = var_fun,
    dist_fun = dist_fun,
    required_args = required_args
  )
}

register_record_model(
  "iid",
  mean_fun = rec_count_mean_iid,
  var_fun  = rec_count_var_iid,
  dist_fun = rec_count_dist_iid
)

register_record_model(
  "DTRW",
  mean_fun = rec_count_mean_DTRW,
  var_fun  = rec_count_var_DTRW,
  dist_fun = rec_count_dist_DTRW
)

register_record_model(
  "LDM",
  mean_fun = rec_count_mean_LDM,
  var_fun  = rec_count_var_LDM,
  dist_fun = rec_count_dist_LDM,
  required_args = c("theta")
)

register_record_model(
  "YNM",
  mean_fun = rec_count_mean_YNM,
  var_fun  = rec_count_var_YNM,
  dist_fun = rec_count_dist_YNM,
  required_args = c("gamma")
)

