#' Fit a bayesian multinomial model to the weekly AGES variant data 
#' 
#' Currently only supports specific variants 
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

#' Extracts the daily growth for Omicron over Delta from the bayesian multinomial model 
#' 
#' @param fit Fit object of the bayesian multinomial model.
#' @return Data frame with central estimates and CIs
extract_bayesian_multinomial_growth_rate <- function(fit) {
  bayestestR::describe_posterior(fit$fit, ci=c(.3, .6, .9), 
                                 parameters="muomicron_days") %>%
    as.data.frame()
}


