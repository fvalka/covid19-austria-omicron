estimate_ve_from_study <- function(ve_file, days, start_day = 7) {
  ve_study_data <- read_csv(ve_file) %>%
    mutate(t=t_start + (t_end - t_start)/2) %>%
    mutate(t=t*7) %>%
    mutate(ve = ve/100)
  
  ve_fit <- glm(ve ~ t, data = ve_study_data, family = gaussian("logit")) 
  
  estimate_ve <- data.frame(t=0:(days - 1),
                            ve = c(rep(0, start_day), predict.glm(ve_fit, newdata = data.frame(t=start_day:(days - 1)), type = "response"))) 
  
  estimate_ve %>%
    left_join(ve_study_data %>% 
                select(t, ve) %>% 
                rename(ve_original_data=ve))
}


estimate_population_immunity <- function(vaccination_data) {
  pop_austria <- 8932664
  
  bnt_bnt_ve <- estimate_ve_from_study(here("data/ve_data/bnt_dose_2.csv"), length(vaccination_data$date), start_day = 14)
  bnt_bnt_bnt_ve <- estimate_ve_from_study(here("data/ve_data/bnt_bnt_boost.csv"), length(vaccination_data$date))
  
  unvaccinated_data <- vaccination_data %>% 
    filter(state_id != 10) %>%
    group_by(date, dose_number) %>%
    summarise(doses_administered_cumulative=sum(doses_administered_cumulative)) %>%
    mutate(dose_number = sprintf("cumulative_dose_%d", dose_number)) %>%
    reshape2::dcast(date ~ dose_number, value.var = "doses_administered_cumulative") %>%
    mutate(unvaccinated = pop_austria - cumulative_dose_1)  %>%
    arrange(date)
  
  estimate_population_immunity_day <- function(vaccination_data, start_date) {
    # Estimate waned immunity
    dose_status <- vaccination_data %>% 
      filter(as.Date(date) <= start_date) %>%
      filter(state_id != 10) %>%
      group_by(date, dose_number) %>%
      summarise(doses_administered_cumulative=sum(doses_administered_cumulative)) %>%
      group_by(dose_number) %>%
      mutate(doses_daily_change = doses_administered_cumulative - lag(doses_administered_cumulative, 1)) %>%
      mutate(doses_daily_change_7d = RcppRoll::roll_meanr(doses_daily_change, 7)) %>%
      mutate(dose_number = sprintf("dose_%d", dose_number)) %>%
      reshape2::dcast(date ~ dose_number, value.var = "doses_daily_change") %>%
      filter(!is.na(dose_1)) %>%
      mutate(dose_1_unadj = dose_1,
             dose_2_unadj = dose_2)
    
    for (t in 2:length(dose_status$date)) {
      new_2nd_doses <- dose_status$dose_2[t]
      
      for (j in 1:(t-14)) {
        min_doses <- dose_status$dose_1_unadj[j] * 0.05
        adjustment <- max(min(dose_status$dose_1[j] - min_doses, new_2nd_doses), 0)
        
        dose_status$dose_1[j] <- dose_status$dose_1[j] - adjustment
        new_2nd_doses <- new_2nd_doses - adjustment
      }
      
      if(new_2nd_doses > 0) {
        warning("2nd Doses not fully accounted for")
      }
      
      new_3rd_doses <- dose_status$dose_3[t]
      for (j in 1:t) {
        min_doses <- dose_status$dose_2_unadj[j] * 0.2
        adjustment <- max(min(dose_status$dose_2[j] - min_doses, new_3rd_doses), 0)
        
        dose_status$dose_2[j] <- dose_status$dose_2[j] - adjustment
        new_3rd_doses <- new_3rd_doses - adjustment
      }
      
      if(new_3rd_doses > 0) {
        warning("3rd Doses not fully accounted for")
      }
      
    }
    
    infection_estimate <- 0.24
    infection_protection <- 0.19
    
    immunity_dose_3 <- sum(dose_status$dose_3 * rev(bnt_bnt_bnt_ve$ve[1:length(dose_status$date)])) # No increase for prior infection
    immunity_dose_2_uninfected <- sum((1-infection_estimate) * dose_status$dose_2 * rev(bnt_bnt_ve$ve[1:length(dose_status$date)]))
    immunity_dose_2_infected <- sum(infection_estimate * dose_status$dose_2 * rev(bnt_bnt_bnt_ve$ve[1:length(dose_status$date)])) # Assume booster VE
    immunity_dose_1_infected <- sum(infection_estimate * dose_status$dose_1 * rev(bnt_bnt_ve$ve[1:length(dose_status$date)])) # Assume two dose VE
    immunity_infection <- tail(unvaccinated_data$unvaccinated, 1) * infection_estimate * infection_protection
    
    data.frame(
      date = as.Date(max(dose_status$date)),
      susceptible = 1 - (immunity_dose_3 + immunity_dose_2_uninfected + 
                           immunity_dose_2_infected + immunity_dose_1_infected + 
                           immunity_infection)/pop_austria
    )
  }
  
  immunity_over_time <- seq(as.Date("2021-11-01"), as.Date(max(vaccination_data$date)), 1) %>%
    purrr::map_dfr(~ estimate_population_immunity_day(vaccination_data, start_date = .))
  
  immunity_over_time
}
