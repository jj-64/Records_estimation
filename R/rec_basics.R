
## Function to get the indicator vector
#' Indicator Function
#'
#' Compute the indicator function of records in a series. "1" will indicate the entry is a record, "0" is not.
#'
#' @param X the series Xt  a vector of length greater or equal to one
#'
#' @returns the indicator vector of same length of the original X vector
#' @export
#'
#' @examples
#' Xt = c(-0.5, -0.4,  0.2, 0.2, -1.0,  0.1,  0.8, -0.2, -0.9,  0.8)
#' is_rec(X=Xt)
#' # [1] 1 1 1 0 0 0 1 0 0 0
is_rec = function(X){
  d=1
  if (length(X) >1 ){
  for(i in 2:length(X)){
    d[i] = ifelse(X[i]>max(X[1:(i-1)]),1,0)
  }}
  return(d)
}
################################################################################

#' Number of records
#'
#' Number of records in a give series, computed by summing the indicator vector
#'
#' @param X the series Xt  a vector of length greater or equal to one
#'
#' @returns integer: the number of records
#' @export
#'
#' @examples
#' Xt = c(-0.5, -0.4,  0.2, 0.2, -1.0,  0.1,  0.8, -0.2, -0.9,  0.8)
#' rec_count(X=Xt)
#' # 4
rec_count = function(X){
  sum(is_rec(X))
}

################################################################################
#' rec_values: Records values
#'
#' Records themselves present in the series.
#' @param X the series Xt  a vector of length greater or equal to one
#'
#' @returns a vector of the record values available in the series
#' @export
#'
#' @examples
#' Xt = c(-0.5, -0.4,  0.2, 0.2, -1.0,  0.1,  0.8, -0.2, -0.9,  0.8)
#' rec_values(X=Xt)
#' # [1] -0.5 -0.4  0.2  0.8
rec_values = function(X){  ## previously Rn
  X[which(is_rec(X)==1)]
}
########################################################################
#' rec_times: Record occurrences
#'
#' The position of records in a series. Always start by 1 referring to initial trivial record
#' @param X the series Xt  a vector of length greater or equal to one
#'
#' @returns a vector of the position of records in a series
#' @export
#'
#' @examples
#' Xt = c(-0.5, -0.4,  0.2, 0.2, -1.0,  0.1,  0.8, -0.2, -0.9,  0.8)
#' rec_times(X=Xt)
#' # [1] 1 2 3 7
rec_times = function(X){  ## previously Ln
  which(is_rec(X)==1)
}
########################################################################
## Inter records times
#' rec_gaps: Inter-record times
#'
#' Inter-record times, i.e the time between two consecutive records. "1" means that there exist one non-record value between two records.
#' @param X the series Xt  a vector of length greater or equal to one
#'
#' @returns a vector of the inter-record times
#' @export
#'
#' @examples
#' Xt = c(-0.5, -0.4,  0.2, 0.2, -1.0,  0.1,  0.8, -0.2, -0.9,  0.8)
#' rec_gaps(X=Xt)
#' # [1] 1 1 4
rec_gaps =  function(X){   ## previously rec_gaps
  #Ln_pos = Rn_time(X)
  # Ln_delta_pos=0
  # for(i in 2:length(Ln_pos)){
  #   Ln_delta_pos[i-1] = Ln_pos[i]-Ln_pos[i-1]}
  return(diff(rec_times(X)))
}
########################################################################
