source(here("R/functions/plots/ggplot2_step_ribbon.R"))
source(here("R/functions/projections/epidemia_projections.R"))

#' Plot a combined plot of Rt, cases, extrapolated cases and daily infections
#' 
#' @param Enriched epidemia fit
#' @param ci_levels Credible interval levels to include inplots
#' 
#' @return Ggplot object 
plot_epidemia_full_combination <- function(fit, ci_levels=c(30, 60, 90)) {
  plot_rt <- epidemia::plot_rt(fit$fit, step = T, levels = ci_levels) + 
    ggpubr::theme_pubr() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("Reproduction Number") +
    labs(y="Effective reproduction number")
  
  plot_obs <- epidemia::plot_obs(fit$fit, type = "cases", levels = ci_levels)  + 
    scale_y_continuous(labels = sitools::f2si) +
    ggpubr::theme_pubr() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("Weekly assigned cases") +
    labs(y="Weekly cases")
  
  plot_obs_extrapolated <- epidemia::plot_obs(fit$fit, 
                                              type = "cases", 
                                              levels = ci_levels, 
                                              newdata = fit$original_data %>%
                                                mutate(sampled=1)) +
    scale_y_continuous(labels = sitools::f2si) +
    ggtitle("Weekly extrapolated cases") +
    ggpubr::theme_pubr() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(y="Weekly cases")
  
  plot_infections <- epidemia::plot_infections(fit$fit, 
                                               levels = ci_levels) + 
    scale_y_continuous(labels = sitools::f2si) +
    ggpubr::theme_pubr() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("Daily infections") +
    labs(y="Daily infections")
  
  (plot_rt | plot_obs) /
    (plot_infections | plot_obs_extrapolated) + 
    plot_annotation(tag_levels = 'A') &
    coord_cartesian(xlim=c(fit$epi_params$ostart, NA))
}


#' Plot an epidemia projection result for two variants
plot_epidemia_projection <- function(data, 
                                     fit_reference, 
                                     fit_new_variant, 
                                     projection, 
                                     cases_ems = NULL, 
                                     start_date = as.Date("2021-10-14"),
                                     by_100k = T,
                                     annotate_max_cases = T,
                                     ylim = NULL) {
  pop <- fit_reference$fit$data$pop[1]
  
  if(by_100k) {
    data_scaling_factor <- (1e5/pop)
  } else {
    data_scaling_factor <- 1
  }
  
  plot_result <- projection$quantiles %>%
    mutate(lower=lower * data_scaling_factor,
           upper =  upper * data_scaling_factor) %>%
    ggplot(aes(x=date)) +
    geom_ribbon(aes(ymin= lower, ymax=upper, alpha=factor(level), fill=group)) +
    geom_vline(xintercept = max(data$Kalenderwoche), linetype="dashed") +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x="Date", y="Weekly cases per 100k", fill="Variant", alpha="CrI", color=NULL) +
    ggsci::scale_fill_nejm(limits=unique(data$variant)) +
    scale_alpha_manual(values = rev(c(0.3, 0.5, 0.7)), limits=factor(c(30, 60, 90))) +
    scale_color_manual(values = c("black", "purple")) +
    scale_x_date(date_breaks = "week") +
    coord_cartesian(xlim = c(start_date, NA), ylim = ylim)
  
  if(!is.null(cases_ems)) {
    plot_result <- plot_result + 
      geom_point(data=cases_ems, aes(x=Time, y=AnzahlFaelle7Tage * data_scaling_factor, color="All observed cases"), stat="identity")
  }
  
  if(annotate_max_cases & !is.null(cases_ems)) {
    plot_result <- plot_result  +
      annotate("rect", xmin = start_date, 
               xmax = max(projection$quantiles$date), 
               ymin = max(cases_ems$AnzahlFaelle7Tage * data_scaling_factor),
               ymax = max(projection$quantiles$upper  * data_scaling_factor),
               fill = "black",
               alpha = 0.1) +
      geom_hline(yintercept = max(cases_ems$AnzahlFaelle7Tage * data_scaling_factor), linetype="dotted")
  }
  
  plot_result
}

#' Plots all generation time scenarios in a single facet_wraped plot
#' 
#' @param data Preprocessed data
#' @param fits All generation time scenario fits
#' @param extend_days Days to project after the last date of data in the model
#' @param ci_levels Credible interval levels shown on plot
plot_epidemia_generation_time_scenarios <- function(data, fits, 
                                                    extend_days = 30, 
                                                    ci_levels = c(30, 60, 90)) {
  first_fit <- fits[[1]]
  pop <- first_fit$fit$data$pop[1]
  
  newdata = epidemia_extend_cases(fit = first_fit, extend_days = extend_days)
  
  project_infections <- function(fit) {
    projected_infections <- posterior_infections(fit$fit, newdata = newdata)
    
    get_quantiles(projected_infections, levels=ci_levels) %>%
      mutate(gen_time_mean = fit$scenario$gen_time_mean,
             gen_time_sd = fit$scenario$gen_time_sd)
  }
  
  infections_all <- fits %>%
    purrr::map_dfr(~ project_infections(.))
  
  infections_all %>% 
    mutate(lower = lower * (1e5/pop)) %>%
    mutate(upper = upper * (1e5/pop)) %>%
    mutate(plot_label = sprintf("μ=%.2f, σ=%.2f", gen_time_mean, gen_time_sd)) %>%
    ggplot(aes(x=date)) +
    geom_ribbon(aes(ymin= lower, ymax=upper, alpha=factor(level), fill=group)) +
    scale_x_date(limits = c(as.Date("2021-12-15"), max(newdata$date))) +
    ggsci::scale_fill_nejm(limits=unique(data$variant)) +
    scale_alpha_manual(values = rev(c(0.3, 0.5, 0.7)), limits=factor(c(30, 60, 90))) +
    facet_wrap(~ plot_label) +
    ggpubr::theme_pubr() +
    labs(x="Date", y="Daily infections per 100k", fill="Variant", alpha="CrI") 
}

plot_epidemia_rt_projection <- function(fit, extend_days = 30, ci_levels=c(30, 60, 90)) {
  newdata <-  epidemia_extend_cases(fit = fit, extend_days = extend_days)
  
  rt_projections <- epidemia::posterior_rt(fit$fit, newdata=newdata)
  
  rt_quantiles <- get_quantiles(rt_projections, levels=ci_levels)  
  
  rt_quantiles %>% 
    ggplot(aes(x=date)) +
    geom_stepribbon(aes(ymin= lower, ymax=upper, alpha=factor(level), fill=group)) +
    geom_hline(yintercept = 1, linetype="dashed", alpha=0.6) +
    scale_x_date(limits = c(as.Date("2021-12-15"), max(newdata$date))) +
    ggsci::scale_fill_nejm(limits=unique(data$variant)) +
    scale_alpha_manual(values = rev(c(0.3, 0.5, 0.7)), limits=factor(c(30, 60, 90))) +
    ggpubr::theme_pubr() +
    labs(x="Date", y="Reproduction number", fill="Variant", alpha="CrI") 
  
}

plot_epidemia_scenarios <- function() {
  tar_load("epidemia_fit_generation_time_scenarios")
  tar_load("epidemia_fit_delta")
  
  fit_delta <- tar_read("epidemia_fit_delta")
  fit_omicron <- tar_read("epidemia_fit_generation_time_scenarios")
  
  rt_offset <- -10
  
  rt_delta_draws <- epidemia::posterior_rt(fit_delta$fit)
  
  process_single_fit <- function(fit) {
    rt_omicron_draws <- epidemia::posterior_rt(fit$fit)
    
    multiplicative_advantage <- rt_omicron_draws$draws[, length(rt_omicron_draws$time) + rt_offset] / rt_delta_draws$draws[, length(rt_delta_draws$time) + rt_offset]
    
    data.frame(
      name = sprintf("mu: %.2f, sigma: %.2f", fit$scenario$gen_time_mean, fit$scenario$gen_time_sd),
      advantage = multiplicative_advantage
    )
  }
  
  advantages <- fit_omicron %>%
    purrr::map_dfr(~ process_single_fit(.))
  
  advantages %>%
    ggplot(aes(x=advantage, y=name)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    coord_cartesian(xlim = c(0, NA))
}