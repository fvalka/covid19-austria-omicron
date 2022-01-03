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
                        generation_time = EuropeCovid$si,
                        pop_adjust = TRUE,
                        pop = 8932664,
                        i2o_cases = epidemia_param_data$i2o_cases_original,
                        pop_immunity = NULL) {
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
  epistart <- ostart - 14
  
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
    mutate(cases = replace(cases, which(date < ostart), NA)) %>%
    mutate(pop = pop)
  
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
  obs <-  epiobs(
    formula = cases ~ 0 + sampled,
    prior = rstanarm::normal(location=0.6, scale=0.2), # IAR
    prior_intercept = rstanarm::normal(location=0.6, scale=0.01),
    link = "identity",
    i2o = i2o2week(i2o_cases)
  )
  
  # Previous immunity
  if(!is.null(pop_immunity)) {
    start_immunity <- pop_immunity %>%
      filter(date == epistart)
    
    change_immunity <- pop_immunity %>%
      mutate(rm_col = pmax((1 - susceptible) - lag(1 - susceptible, 1), 0)) %>%
      filter(!is.na(rm_col))
    
    variants_model <- variants_model %>%
      left_join(change_immunity)
    
    rm_col <- "rm_col"
    
    inf <- epiinf(gen = generation_time, 
                  seed_days = 6, 
                  pop_adjust = pop_adjust, 
                  pops = "pop", 
                  prior_susc = rstanarm::normal(location = start_immunity$susceptible, scale = 0.1), 
                  rm = rm_col)
  } else {
    inf <- epiinf(gen = generation_time, 
                  seed_days = 6, 
                  pop_adjust = pop_adjust, 
                  pops = "pop")
  }
  
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


epidemia_param_data <- list(
  i2o_cases_original = read_csv(here("data/i2o/i2o_cases.csv"))$i2o,
  # Fit to original and recuded incubation time distr. mean and sd by factor 
  # from South Korea serial interval factor between Omicron and Delta
  i2o_omicron_si_shortened = read_csv(here("data/i2o/i2o_cases_omicron.csv"))$i2o 
)
