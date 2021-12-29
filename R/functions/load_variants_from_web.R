library(rvest)
library(dplyr)
library(here)

load_variants_from_web <- function(variants_url = "https://www.ages.at/themen/krankheitserreger/coronavirus/sars-cov-2-varianten-in-oesterreich/") {
  ages_varianten <- read_html(variants_url)
  
  tables <- ages_varianten %>%
    html_table()
  
  variant_cases <- tables[[1]] %>% 
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
    mutate(Kalenderwoche = ISOweek::ISOweek2date(sprintf("2021-W%02d-7", Kalenderwoche))) 
}

load_all_variants <- function() {
  variants_current <- load_variants_from_web()
  variants_archive <- load_variants_from_web("https://web.archive.org/web/20211208215839/https://www.ages.at/themen/krankheitserreger/coronavirus/sars-cov-2-varianten-in-oesterreich/")
  
  variants_all <- rbind(variants_current,
                        variants_archive %>% 
                          filter(Kalenderwoche < min(variants_current$Kalenderwoche)))
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
