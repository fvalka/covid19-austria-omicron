source(here("R/functions/model/multinomial_models.R"))

#' Extrapolate daily 7d cases from multinomial fit
#' 
#' @param data Preprocessed data
#' @param fit Multinomial fit (brms)
#' @param cases_ems_austria Daily case data
#' 
#' @return Dataframe in format of data with extrapolated cases
multinomial_extrapolated_cases <- function(data, fit, cases_ems_austria, tail_offset = 1) {
  multinomial_samples <- draw_samples_from_multinomial_fit(data, fit, extend_days = 90)
  
  cases_all <- cases_ems_austria %>%
    rename(Kalenderwoche=Time,
           cases_all=AnzahlFaelle7Tage) %>%
    select(Kalenderwoche, cases_all) %>%
    filter(Kalenderwoche <= max(Kalenderwoche) - tail_offset)
  
  extrapolate_variant <- function(variant_name, category_name) {
    cases_proportion <- multinomial_samples %>% 
      filter(.category == category_name) %>%
      group_by(Kalenderwoche) %>%
      summarise(proportion = median(prop)) %>%
      mutate(variant = variant_name)
    
    variant_data <- cases_proportion %>%
      inner_join(cases_all) %>%
      mutate(cases_assigned=cases_all,
             cases=round(proportion*cases_all),
             proportion_low = proportion,
             proportion_high = proportion) %>%
      mutate(week = as.integer(ceiling((Kalenderwoche - max(Kalenderwoche))/7))) %>%  
      group_by(week) %>%
      filter(Kalenderwoche == max(Kalenderwoche)) %>%
      ungroup() %>%
      select(Kalenderwoche, variant, cases_all, cases, cases_assigned, proportion, proportion_low, proportion_high)
  }
  
  omicron <- extrapolate_variant("B.1.1.529 (Omikron)", "omicron") 
  
  delta <-extrapolate_variant("B.1.617.2 (Delta)", "delta")
  
  data %>%
    filter(Kalenderwoche < min(delta$Kalenderwoche)) %>%
    rbind(delta) %>%
    rbind(omicron)
  
}

#' Draw samples from the proportions of the brms multinomial fit and 
#' obtain median and 90% eti credible intervals for each variant and day
#' 
#' @param data Preprocessed data object
#' @param fit Multinomial fit object
#' @param extend_days Days to extend past the last day in the fitted data 
#' 
#' @return Data frame with date, variant, proportions and 90% ci
multinomial_proportions <- function(data, fit, extend_days=90) {
  samples <- draw_samples_from_multinomial_fit(data, tar_read("multinomial_fit"), extend_days=extend_days)
  
  samples %>%
    group_by(Kalenderwoche, variant_full_name) %>%
    summarise(prop_median=median(prop),
              prop_ci_low=quantile(prop, 0.05),
              prop_ci_high=quantile(prop, 0.95)) %>%
    rename(variant=variant_full_name,
           date=Kalenderwoche)
}
