#' Fit a Bayesian multinomial GLMM with a weekly 
#' random effect accounting for overdispersion.
#' 
#' Currently only supports the Alpha, Delta and Omicron variants. 
#' 
#' @param data Preprocessed variant data 
#' @param model_start Start date of variant data fitted in the model
#' 
#' @return Fit object, enriched with model specific information
fit_bayesian_multinomial_model <- function(data, 
                                       model_start = as.Date("2021-09-01")) {
  # Variants to include in the model, also determines the order in the multinomial model, first is reference 
  variants_included <- c(delta="B.1.617.2 (Delta)", omicron="B.1.1.529 (Omikron)", alpha="B.1.1.7 (Alpha)")
  
  model_data <- data %>%
    filter(variant %in% variants_included) %>%
    filter(Kalenderwoche >= model_start) %>%
    select(Kalenderwoche, variant, cases) %>%
    reshape2::dcast(Kalenderwoche ~ variant) %>%
    arrange(Kalenderwoche) %>%
    mutate(days=as.numeric(Kalenderwoche)) %>% 
    rename(alpha=`B.1.1.7 (Alpha)`, delta=`B.1.617.2 (Delta)`, omicron=`B.1.1.529 (Omikron)`) %>%
    mutate(cases_assigned = alpha + delta + omicron)
  
  
  fit <- brm(cbind(delta, omicron, alpha) | trials(cases_assigned) ~ days + (1 | days), 
      data=model_data,
      family = multinomial(link = "logit"))
  
  list(
    fit=fit,
    model_data=model_data,
    model_start=model_start
  )
}


#' Draws posterior samples from a multinomial fit object using 
#' [tidybayes::add_epred_draws()].
#' 
#' @param data Preprocessed data
#' @param fit [brms::brm()] multinomial fit
#' @param ndraws Number of draws to draw
#' 
#' @return Dataframe with all draws in the .epred column
draw_samples_from_multinomial_fit <- function(data, fit, ndraws=1000, extend_days=0) {
  full_variant_names <- 
    data.frame(
      .category = c("delta", "omicron", "alpha"),
      variant_full_name = c("B.1.617.2 (Delta)", "B.1.1.529 (Omikron)", "B.1.1.7 (Alpha)"))
  
  extended_dates <- 
    data.frame(
      Kalenderwoche =seq(min(fit$model_data$Kalenderwoche), max(fit$model_data$Kalenderwoche) + extend_days, 1)
    ) %>%
    mutate(days=as.numeric(Kalenderwoche))
  
  data %>%
    filter(Kalenderwoche >= fit$model_start) %>%
    full_join(extended_dates) %>%
    arrange(Kalenderwoche) %>%
    tidyr::fill(cases_assigned) %>%
    mutate(days=as.numeric(Kalenderwoche)) %>%
    tidybayes::add_epred_draws(fit$fit, ndraws = ndraws, allow_new_levels=T, re_formula = NA) %>%
    mutate(prop = .epred/cases_assigned) %>%
    left_join(full_variant_names)
}


#' Extracts the daily growth for Omicron over Delta from the Bayesian multinomial model 
#' 
#' @param fit Fit object of the bayesian multinomial model.
#' @return Data frame with central estimates and CIs
extract_bayesian_multinomial_growth_rate <- function(fit) {
  bayestestR::describe_posterior(fit$fit, ci=c(.3, .6, .9), 
                                 parameters="muomicron_days") %>%
    as.data.frame()
}


