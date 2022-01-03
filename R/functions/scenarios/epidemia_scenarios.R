#' Fit a single generation time scenario based upon a mean and standard deviation
#' of a discretized gamma distribution. 
#' 
#' The functions used for discretization is [EpiEstim::discr_si()].
#' 
#' @param data Preprocessed data
#' @param gen_time_mean Mean generation time (assumed discretized gamma distributed)
#' @param gen_time_sd Standard deviation of generation time (assumed discretized gamma distributed)
#' @param scenario_name Name of the scenario 
#' 
#' @return Enriched fit object as in [epidemia_fit_variant()]. Also includes a scenario list element with parameters set in the function call.
epidemia_generation_time_scenario <- function(data, 
                                              gen_time_mean, 
                                              gen_time_sd, 
                                              scenario_name = "",
                                              ...) {
  gen_time_dist <- EpiEstim::discr_si(0:99, mu=gen_time_mean, sigma = gen_time_sd)
  
  # Workaround for tiny remaining probability after 99 days and rounding errors with
  # epidemia sum(si) == 1 check
  gen_time_dist <- gen_time_dist/sum(gen_time_dist)
  
  fit <- epidemia_fit_variant(variants_all = data,
                              generation_time = gen_time_dist,
                              ...)
  
  fit$scenario <- list(
    name = scenario_name,
    gen_time_dist = gen_time_dist,
    gen_time_mean = gen_time_mean,
    gen_time_sd = gen_time_sd
  )
  
  fit
}