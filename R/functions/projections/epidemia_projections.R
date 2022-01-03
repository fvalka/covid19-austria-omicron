#' Produce a new data frame for usage in newdata for epidemia
#' e.g. in [epidemia::posterior_predict()].
#'
#' @param fit Enriched epidemia fit
#' @param extend_days Days to extend the data frame after the last date included in the modeled data
#' @param sampled_cases Proportion of cases sampled for this variant, use 1 for extrapolating all cases
#'
#' @return Dataframe for usage in newdata argument of epidemia plot_\* and predict_\* functions functions
epidemia_extend_cases <-
  function(fit,
           extend_days = 60,
           sampled_cases = 1,
           rm_per_day = 0) {
    data <- fit$original_data
    
    pop_modeled <- fit$fit$data$pop[1]
    
    newdata <- data.frame(
      date = max(data$date) + seq_len(extend_days),
      variant = fit$variant,
      cases_all = NA,
      cases = -1,
      cases_assigned = NA,
      proportion = NA,
      proportion_low = NA,
      proportion_high = NA,
      cases_extrapolated = -1,
      sampled = sampled_cases,
      pop = pop_modeled
    ) 
    
    if("rm_col" %in% colnames(data)) {
      newdata$susceptible <- NA
      newdata$rm_col <- rm_per_day
    }
    
    newdata %>%
      mutate(week = as.integer(ceiling((
        date - max(data$date)
      ) / 7))) %>%
      rbind(data) %>%
      mutate(sampled = sampled_cases)
  }

epidemia_projection_cases_generator <- function(fit, newdata) {
  posterior_predict(fit$fit, type = "cases", newdata = newdata, draws=2000)
}

epidemia_projection_infections_generator <- function(fit, newdata) {
  posterior_infections(fit$fit, newdata = newdata, draws = 2000)
}

#' Project for a specific posterior draw function 
#' 
#' @param fit_reference Enriched epidemia fit object for reference variant
#' @param fit_new_variant Enriched epidemia fit object for new variant 
#' @param extend_days Days to project after the last model input data date
#' @param ci_levels Credible interval levels for quantiles, e.g. 90 for 90% CrI
#' 
#' @return Nested list object with quantiles and draws for reference variant, new variant and a combination of both
epidemia_project <- function(fit_reference,
                             fit_new_variant,
                             extend_days = 60,
                             ci_levels = c(30, 60, 90),
                             predictor = epidemia_projection_cases_generator) {
  reference_newdata <- epidemia_extend_cases(fit = fit_reference, extend_days = extend_days)
  new_variant_newdata <- epidemia_extend_cases(fit = fit_new_variant, extend_days = extend_days)
  
  reference_pred <- epidemia_projection_cases_generator(fit_reference)
  new_variant_pred <- epidemia_projection_cases_generator(fit_new_variant)
  
  t_fill <- length(reference_pred$time) - length(new_variant_pred$time)
  draws <- dim(reference_pred$draws)[1]
  
  new_variant_draws_filled <- 
    cbind(matrix(rep(0, draws * t_fill), nrow = draws, ncol = t_fill),
          new_variant_pred$draws)
  
  combined_pred <- reference_pred
  combined_pred$draws <- combined_pred$draws + new_variant_draws_filled
  combined_pred$group <- "Combined"
  
  reference_quantiles_pred <- get_quantiles(reference_pred, levels = ci_levels)
  new_variant_quantiles_pred <- get_quantiles(new_variant_pred, levels = ci_levels)
  combined_quantiles_pred <- get_quantiles(combined_pred, levels = ci_levels)
  
  list(
    quantiles = rbind(
      reference = reference_quantiles_pred,
      new_variant = new_variant_quantiles_pred
    ),
    quantiles_combined = combined_quantiles_pred,
    variants = list(
      reference = reference_quantiles_pred$group[1],
      new_variant = new_variant_quantiles_pred$group[1]
    ),
    draws = list(
      reference = reference_pred,
      new_variant = new_variant_pred,
      combined = combined_pred
    )
  )
}


# Compute quantiles for all levels
#
# Taken from epidemia open source package plots_epi.R
#
# @param object Result of a posterior_ function
# @param levels A numeric vector defining levels
get_quantiles <- function(object, levels, dates=NULL, date_format=NULL) {
  levels <- levels[order(levels)]
  f <- function(level) {
    res <- apply(
      object$draws,
      2,
      function(x) quantile(x, 0.5 + level * c(-1, 1) / 200)
    )
    return(
      data.frame(
        date = object$time,
        lower = res[1, ],
        upper = res[2, ],
        group = object$group,
        tag = paste0(level, "% CI"),
        level = level
      )
    )
  }
  out <- lapply(levels, f)
  out <- do.call(rbind, out)
  out$tag <- factor(out$tag, ordered = T, levels = rev(levels(factor(out$tag))))
  if (!is.null(dates)){
    out <- subset_for_dates(
      out,
      dates,
      date_format
    )
  }
  return(out)
}
