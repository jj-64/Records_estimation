library(numDeriv)
#
# compute_var_hessian <- function(loglik_fun, data, par_vec) {
#
#   hess <- tryCatch({
#     numDeriv::hessian(
#       func = function(p) {
#         params <- list(theta = p["theta"], shape = p["shape"], scale = p["scale"])
#         loglik_fun(data = data, params = params)
#       },
#       x = par_vec
#     )
#   }, error = function(e) return(matrix(NA, 3, 3)))
#
#   # Check if valid
#   if (any(is.na(hess)) || any(!is.finite(hess))) {
#     return(c(NA, NA, NA))
#   }
#
#   # Variance = inverse negative Hessian
#   vcov <- tryCatch({
#     solve(-hess)
#   }, error = function(e) return(matrix(NA, 3, 3)))
#
#   if (any(is.na(vcov))) {
#     return(c(NA, NA, NA))
#   }
#
#   return(diag(vcov))  # return variances only
# }

# par_hat <- c(
#   MLE_C["theta"],
#   MLE_C["shape"],
#   MLE_C["scale"]
# )
#
# var_emp <- compute_var_hessian(logLik_fun_rec, data_rec, par_hat)
compute_var_hessian <- function(loglik_fun, data, par_vec, transform = NULL, eps = 1e-6, method = "Richardson") {

  # Ensure par_vec is named
  if (is.null(names(par_vec))) {
    stop("par_vec must be a named vector")
  }

  # Wrapper function for numDeriv
  func <- function(p) {
    names(p) <- names(par_vec) # Guarantee names stay intact during numerical steps

    # Apply transformations if provided
    if (!is.null(transform)) {
      # Ensure transform is named or ordered identically to p
      params <- Map(function(val, f) {
        if (is.null(f)) val else f(val)
      }, as.list(p), transform[names(p)]) # Match by name explicitly
    } else {
      params <- as.list(p)
    }

    val <- tryCatch({
      loglik_fun(data = data, params = params)
    }, error = function(e) NA_real_)

    if (!is.finite(val)) return(NA_real_)
    return(val)
  }

  # Compute Hessian
  hess <- tryCatch({
    numDeriv::hessian(
      func = func,
      x = par_vec,
      method = method, # "Simple" is less prone to boundary NA issues than Richardson
      method.args = list(eps = eps)
    )
  }, error = function(e) {
    warning("Hessian calculation threw an error: ", e$message)
    return(NULL)
  })

  # Diagnostic check: which element is causing the NA?
  if (is.null(hess) || any(!is.finite(hess))) {
    warning("Hessian contains non-finite values (NA/NaN/Inf). Surface might be flat or at a strict boundary.")
    return(setNames(rep(NA_real_, length(par_vec)), names(par_vec)))
  }

  # Try inversion
  # vcov <- tryCatch({
  #   solve(-hess)
  # }, error = function(e) {
  #   warning("Hessian inversion failed (matrix might be singular). Trying pseudo-inverse next time.")
  #   return(NULL)
  # })

  vcov <- tryCatch({
    solve(-hess)
  }, error = function(e) {
    warning("Matrix singular. Falling back to Moore-Penrose pseudo-inverse.")
    # ginv handles near-singular matrices smoothly
    MASS::ginv(-hess)
  })

  # Validate vcov
  if (is.null(vcov) || any(!is.finite(vcov))) {
    return(setNames(rep(NA_real_, length(par_vec)), names(par_vec)))
  }

  # Return named variances
  return(setNames(diag(vcov), names(par_vec)))
}
#__________________________#

compute_var_hessian_optim <- function(loglik_fun, data, par_vec) {

  hess <- tryCatch({
    optimHess(
      par = par_vec,
      func = function(p) {
        loglik_fun(data = data, params = list(p) )
      }
    )
  }, error = function(e) return(matrix(NA, 3, 3)))

  # Check if valid
  if (any(is.na(hess)) || any(!is.finite(hess))) {
    return(c(NA, NA, NA))
  }

  # Variance = inverse negative Hessian
  vcov <- tryCatch({
    solve(-hess)
  }, error = function(e) return(matrix(NA, 3, 3)))

  if (any(is.na(vcov))) {
    return(c(NA, NA, NA))
  }

  return(diag(vcov))  # return variances only
}

## compute hessian matrix and vcov and correlation ------------
# var_emp <- compute_var_hessian(logLik_fun_rec, data_rec, par_hat)
compute_hessian <- function(loglik_fun, data, par_vec, transform = NULL, eps = 1e-6, method = "Richardson") {

  # Ensure par_vec is named
  if (is.null(names(par_vec))) {
    stop("par_vec must be a named vector")
  }

  # Wrapper function for numDeriv
  func <- function(p) {
    names(p) <- names(par_vec) # Guarantee names stay intact during numerical steps

    # Apply transformations if provided
    if (!is.null(transform)) {
      # Ensure transform is named or ordered identically to p
      params <- Map(function(val, f) {
        if (is.null(f)) val else f(val)
      }, as.list(p), transform[names(p)]) # Match by name explicitly
    } else {
      params <- as.list(p)
    }

    val <- tryCatch({
      loglik_fun(data = data, params = params)
    }, error = function(e) NA_real_)

    if (!is.finite(val)) return(NA_real_)
    return(val)
  }

  # Compute Hessian
  hess <- tryCatch({
    numDeriv::hessian(
      func = func,
      x = par_vec,
      method = method, # "Simple" is less prone to boundary NA issues than Richardson
      method.args = list(eps = eps)
    )
  }, error = function(e) {
    warning("Hessian calculation threw an error: ", e$message)
    return(NULL)
  })

  # Diagnostic check: which element is causing the NA?
  if (is.null(hess) || any(!is.finite(hess))) {
    warning("Hessian contains non-finite values (NA/NaN/Inf). Surface might be flat or at a strict boundary.")
    return(setNames(rep(NA_real_, length(par_vec)), names(par_vec)))
  }

  # Try inversion
  # vcov <- tryCatch({
  #   solve(-hess)
  # }, error = function(e) {
  #   warning("Hessian inversion failed (matrix might be singular). Trying pseudo-inverse next time.")
  #   return(NULL)
  # })

  colnames(hess) = rownames(hess) = names(par_vec)

  vcov <- tryCatch({
    solve(-hess)
  }, error = function(e) {
    warning("Matrix singular. Falling back to Moore-Penrose pseudo-inverse.")
    # ginv handles near-singular matrices smoothly
    MASS::ginv(-hess)
  })

  # Validate vcov
  if (is.null(vcov) || any(!is.finite(vcov))) {
    return(setNames(rep(NA_real_, length(par_vec)), names(par_vec)))
  }

  ## correlation
  covcor = cov2cor(vcov)

  # Return named variances
  return(list(hessian = hess, vcov = vcov, corr = covcor ))
}
