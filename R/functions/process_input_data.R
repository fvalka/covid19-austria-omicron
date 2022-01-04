library(dplyr)
library(here)
library(readr)

process_variant_data <- function(raw_input_data) {
  variant_cases <- raw_input_data %>% 
    select(-filename) %>%
    mutate(Kalenderwoche = stringr::str_replace(Kalenderwoche, "2021-W", "")) %>%
    filter(!is.na(as.numeric(Kalenderwoche))) %>%
    reshape2::melt(id.vars=c("Kalenderwoche", "Fälle gesamt"),
                   variable.name="variant",
                   value.name="cases") 
  
  variant_cases_assigned <- variant_cases %>%
    group_by(Kalenderwoche) %>%
    summarise(cases_assigned=sum(cases))
  
  variant_cases %>%
    left_join(variant_cases_assigned) %>%
    mutate(proportion=cases/cases_assigned) %>%
    mutate(cases_extrapolated = proportion * `Fälle gesamt`) %>%
    mutate(Kalenderwoche = as.numeric(Kalenderwoche)) %>%
    mutate(Kalenderwoche = ISOweek::ISOweek2date(sprintf("2021-W%02d-7", Kalenderwoche))) %>% 
    purrr::pmap_dfr(~ estimate_proportion(..1, ..2, ..3, ..4, ..5))
}

process_ems_cases <- function(file) {
  read_csv2(file) %>%
    filter(BundeslandID==10) %>%
    mutate(Time=as.Date(lubridate::dmy_hms(Time)))
}


estimate_proportion <- function(Kalenderwoche, cases_all, variant,  cases, cases_assigned) {
  binom_estimate <- binom.test(cases, cases_assigned)
  
  data.frame(
    Kalenderwoche=Kalenderwoche,
    variant=variant,
    cases_all=cases_all,
    cases=cases,
    cases_assigned=cases_assigned,
    proportion = binom_estimate$estimate,
    proportion_low = binom_estimate$conf.int[1],
    proportion_high = binom_estimate$conf.int[2]
  )
}

#' Stores the raw data locally for archiving purposes
#' 
#' @param raw_data Raw data as read by read_csv
#' 
#' @return Path where the csv file was written to locally
export_raw_data <- function(raw_data,
                            target_dir="data/ages_varianten/") {
  output_path <- paste0(target_dir, head(raw_data$filename, 1))
  
  write_csv(raw_data %>% select(-filename), file = output_path)
  
  output_path
}


#' Adjust data for assumed sampling bias by a time constant sampling adjustment factor
#' 
#'  @param data Preprocssed data 
#'  @param sampling_adjust Sampling adjustment factor between 0 and 1. At 1 all non-assigned variant are assigned to the selected variant. At 0 none are. 
#'  @param variant_adjusted Variant which is to be imputed based upon the sampling_adjust factor
data_adjust_sampling <- function(data,
                                 sampling_adjust=0,
                                 variant_adjusted="B.1.617.2 (Delta)"
) {
  
  stopifnot(sampling_adjust >= 0, sampling_adjust <= 1)
  
  adjusted_cases <- data %>% 
    mutate(cases=case_when(variant == variant_adjusted ~ cases + round((cases_all - cases_assigned) * sampling_adjust),
                           T ~ cases)) %>%
    select(-cases_assigned)
  
  adjusted_cases_assigned <- adjusted_cases %>%
    group_by(Kalenderwoche) %>%
    summarise(cases_assigned = sum(cases))
  
  adjusted_cases %>%
    left_join(adjusted_cases_assigned)
  
}

