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

fit_betabinomial_glm <- function(data) {
  beta_binomial2 <- 
    custom_family(
      "beta_binomial2", dpars = c("mu", "phi"),
      links = c("logit", "log"), lb = c(NA, 0),
      type = "int", vars = "trials[n]"
    )
  
  stan_funs <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }
"
  
  stanvars <- 
    stanvar(scode = stan_funs, block = "functions")
  
  binomial_fit <- brm(cases | trials(cases_assigned) ~ date,
                      family = beta_binomial2,
                      stanvars = stanvars,
                      data=model_data)
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