#' Fit a case only epidemia model to a specific variant 
#' 
#' @param variants_investigated Name of the variant as in the full variant report / data frame
#' @param variants_all Full data frame with all variant data as per the AGES 
#' 
#' @return Enriched epidemia fit object 
epidemia_fit_variant <- function(variants_investigated = "B.1.1.529 (Omikron)",
                        variants_all,
                        ostart = NA,
                        min_cases = 10,
                        generation_time = EuropeCovid$si) {
  #' Converts a daily infection to observation distribution to a weekly one 
  #' 
  #' @param i2o Daily infection to observation distribution 
  i2o2week <- function(i2o)
    rowSums(sapply(0:6, function (k) c(rep(0,k),i2o,rep(0,6-k))))
  
  if(is.na(ostart)) {
    # Calculate start date 
    variant_first_cases <- variants_all %>% 
      filter(variant == variants_investigated) %>% 
      filter(cases >= min_cases) %>% 
      summarise(date=min(Kalenderwoche))
    
    # Start of cases in model
    ostart <- variant_first_cases$date - 3
  }
  
  # First day in model
  epistart <- ostart - 20
  
  # Extrapolate the number of variant cases from the number of total cases
  # and the proportion of the variant cases of all cases assigned to any variant
  variants_extrapolated <- variants_all %>%
    filter(Kalenderwoche >= epistart) %>%
    mutate(cases_extrapolated = cases_all*proportion) %>%
    filter(variant == variants_investigated)
  
  # List of all days to fill in missing days later
  all_days <- data.frame(date=seq(min(variants_extrapolated$Kalenderwoche), max(variants_extrapolated$Kalenderwoche), 1)) %>%
    merge(data.frame(variant=unique(variants_extrapolated$variant)))
  
  variants_model <- variants_extrapolated %>%
    rename(date=Kalenderwoche) %>%
    right_join(all_days) %>%
    mutate(week = as.integer(ceiling((date - max(date))/7))) %>%
    arrange(date, variant) %>%
    mutate(cases_extrapolated = round(cases_extrapolated)) %>%
    mutate(sampled=cases_assigned/cases_all) %>%
    mutate(cases_extrapolated = replace(cases_extrapolated, which(date < ostart), NA))  %>%
    mutate(cases = replace(cases, which(date < ostart), NA)) 
  
  # Start random walk 2 weeks before first variant cases
  variants_model$week <- pmax(variants_model$week,variants_model$week[which(variants_model$date==ostart - 14)])
  
  
  rt <- epirt(
    formula = R(variant, date) ~ 1 + rw(time = week, prior_scale = 0.3) +
      I(date>="2021-04-01") +
      I(date>="2021-11-08") +
      I(date>="2021-11-15") +
      I(date>="2021-11-22"),
    prior_intercept = rstanarm::normal(log(1), log(3)),
    link = 'log' # scaled_logit(7)
  )
  
  # Observations model input
  i2o_cases <- read_csv(here("data/i2o/i2o_cases.csv"))
  
  obs <-  epiobs(
    formula = cases ~ 0 + sampled,
    prior = rstanarm::normal(location=0.6, scale=0.2), # IAR
    prior_intercept = rstanarm::normal(location=0.6, scale=0.01),
    link = "identity",
    i2o = i2o2week(i2o_cases$i2o)
  )
  
  inf <- epiinf(gen = generation_time, seed_days = 6, pop_adjust = FALSE)
  
  args <- list(rt=rt, inf=inf, obs=obs, data=variants_model, seed=12345)
  
  options(mc.cores = parallel::detectCores())
  
  args$algorithm <- "sampling"
  args$iter <- 3000
  args$control = list(max_treedepth = 12)
  
  fm <- do.call(epim, args)
  
  list(
    fit=fm,
    variant=variants_investigated,
    original_data=variants_model,
    epi_params=list(
      ostart=ostart,
      epistart=epistart
    )
  )
}

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
                                              scenario_name = "") {
  gen_time_dist <- EpiEstim::discr_si(0:99, mu=gen_time_mean, sigma = gen_time_sd)
  
  # Workaround for tiny remaining probability after 99 days and rounding errors with
  # epidemia sum(si) == 1 check
  gen_time_dist <- gen_time_dist/sum(gen_time_dist)
  
  fit <- epidemia_fit_variant(variants_all = data,
                       generation_time = gen_time_dist,
                       min_cases = 1)
  
  fit$scenario <- list(
      name = scenario_name,
      gen_time_dist = gen_time_dist,
      gen_time_mean = gen_time_mean,
      gen_time_sd = gen_time_sd
    )
  
  fit
}

