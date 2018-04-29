
#' Return posterior summary of `rstan::sampling()` output
#' @param ... arguments to `rstan::sampling(...)`
#' @export
sampling_summary <- function(...){
  fit <- rstan::sampling(...)
  return(rstan::summary(fit)$summary)
}
