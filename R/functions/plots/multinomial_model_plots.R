source(here("R/functions/plots/utils.R"))

#' Plot a [brms::brm()] multinomial fit.
#' 
#' Including a Bayesian credible interval based upon sampling from the posterior 
#' distribution.
#' 
#' @param data Preprocessed data 
#' @param fit Multinomial model fit result object
#' @param extend_days Project by this amount of days into the future
#' @param show_zoom Show a zoom inset on the plot
#' 
#' @return A ggplot object with the multinomial fit 
plot_multinomial_fit <- function(data, 
                                 fit, 
                                 extend_days = 0,
                                 show_zoom = T,
                                 ci_levels = c(0.3, 0.6, 0.9)) {
  
  samples <- draw_samples_from_multinomial_fit(data, fit, extend_days=extend_days)

  complete_plot <- samples %>%
    ggplot(aes(x=Kalenderwoche)) +
    tidybayes::stat_lineribbon(aes(y=prop, color=variant_full_name), .width = ci_levels, alpha=0.25) +
    geom_point(aes(y=proportion, color=variant), alpha=0.7) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("#adb5bd", "#6c757d", "#495057")) +
    ggsci::scale_color_nejm(limits=unique(data$variant)) +
    ggpubr::theme_pubr() +
    labs(color=NULL, fill="CI",
         x="Date, week ending",
         y="Weekly proportion of cases with reported variant") 
  
  if(show_zoom) {
    complete_plot <- complete_plot +
      ggforce::facet_zoom(xlim = c(as.Date("2021-11-01"), max(data$Kalenderwoche)), ylim = c(0, determine_zoom_y_limit(data)))
  }
  
  complete_plot
}

#' Plots the posterior distribution of the growth advantage of Omicron over Delta 
#' from the Bayesian multinomial model
#' 
#' @param fit Enriched multinomial fit object
#' 
#' @return Ggplot object of the growth advantage posterior 
plot_multinomial_growth_advantage <- function(fit, ci_levels=c(0.3, 0.6, 0.9)) {
  plot(bayestestR::hdi(as.data.frame(fit$fit)$b_muomicron_days, ci_levels)) +
    scale_fill_brewer(palette = "Blues", direction = -1) +
    ggpubr::theme_pubr() +
    labs(fill="CI",
         title = "Growth advantage of Omicron over Delta",
         x="Growth rate advantage",
         y=NULL)
}

plot_multinomial_weekly_fold_increase <- function(fit, ci_levels=c(0.3, 0.6, 0.9)) {
  plot(bayestestR::hdi(exp(as.data.frame(fit$fit)$b_muomicron_days*7), ci_levels)) +
    scale_fill_brewer(palette = "Blues", direction = -1) +
    ggpubr::theme_pubr() +
    labs(fill="CI",
         title = "Weekly fold increase of Omicron over Delta",
         x="Weekly fold increase",
         y=NULL)
}

#' Plot ECDF of more than 50% Omicron cases in the posterior draws
#' 
#' @param data Preprocssed data
#' @param fit [brms::brm()] multinomial model fit 
#' 
#' @return Ggplot object of ECDF of more than 50% Omicron cases#' 
plot_probablity_more_than_half <- function(data, fit)  {
  draw_samples_from_multinomial_fit(data, fit, extend_days=90) %>%
    filter(variant_full_name == "B.1.1.529 (Omikron)") %>%
    filter(prop >= 0.5) %>%
    group_by(.draw) %>%
    summarise(Kalenderwoche = min(Kalenderwoche)) %>%
    ggplot(aes(x=Kalenderwoche)) +
    stat_ecdf(geom = "line") +
    labs(x="Date, week ending",
         y="Probability") +
    scale_x_date(date_labels =  "%d %b %Y") +
    scale_y_continuous(labels = scales::percent) +
    coord_cartesian(xlim = c(as.Date("2021-12-7"), as.Date("2022-01-21"))) +
    ggpubr::theme_pubr() +
    ggtitle("Probability of more than half Omicron cases")
  
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

