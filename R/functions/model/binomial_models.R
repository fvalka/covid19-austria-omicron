#' Fit a Bayesian binomial GLM
#' 
#' @param data Preprocessed variant data
#' 
#' @return BRMS Binomial fit object
fit_binomial_glm <- function(data) {
  model_data <- binomial_model_data_convert(data)
  
  binomial_fit <- brm(cases | trials(cases_assigned) ~ date,
                      family = binomial(logit),
                      data=model_data)
  
  binomial_fit
}

#' Classical quasi binomal (overdispersed binomial) GLM
#' 
#' @param data Preprocessed data
#' 
#' @return GLM fit object
fit_quasibinomial_glm <- function(data) {
  model_data <- binomial_model_data_convert(data)
  
  glm_fit <- glm(cbind(cases, non_omicron_cases) ~ date,
                 family = quasibinomial(logit),
                 data=model_data)
  
  glm_fit
}

#' Helper function for extracting omicron case data for binomial modelling
#' 
#' @param data Preprocessed data
#' 
#' @return Model input data
binomial_model_data_convert <- function(data) {
  data %>% 
    filter(variant == "B.1.1.529 (Omikron)") %>%
    mutate(non_omicron_cases = cases_assigned - cases) %>%
    rename(date=Kalenderwoche) %>%
    filter(cases > 0)
}