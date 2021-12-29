source(here("R/functions/plots/utils.R"))

#' Plots the proportion of variants among all cases assigned to variants per week
#' 
#' @param data Full preprocessed variant data set 
#' 
#' @return Ggplot object of extraploated variant proportion
plot_extraploted_variant_proportion <- function(data) {
  labeled_variants <- c("B.1.1.529 (Omikron)")
  
  data %>% 
    filter(Kalenderwoche >= as.Date("2021-04-01")) %>%
    mutate(label = if_else(Kalenderwoche == max(Kalenderwoche) & variant %in% labeled_variants, sprintf("%.1f%%", proportion*100), NA_character_)) %>%
    ggplot(aes(x=Kalenderwoche, y=proportion)) +
    geom_line(aes( color=variant)) +
    geom_ribbon(aes(ymin=proportion_low, ymax=proportion_high, fill=variant), alpha=0.3) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggsci::scale_color_nejm() +
    ggsci::scale_fill_nejm() +
    labs(x="Date, week ending",
         y="Proportion of cases with known/specified variant",
         color=NULL, 
         fill=NULL) +
    ggforce::facet_zoom(xlim = c(as.Date("2021-11-01"), max(data$Kalenderwoche)), ylim = c(0, determine_zoom_y_limit(data))) +
    ggpubr::theme_pubr()
}

#' Plot an area plot of the variants proportion of the total proportion of
#' all cases assigned to any variant. 
#' 
#' @param data Preprocssed data
#' 
#' @return Gggplot object
plot_variants_area <- function(data) {
  data %>% 
    filter(Kalenderwoche >= as.Date("2021-04-01")) %>%
    ggplot(aes(x=Kalenderwoche, y=proportion, color=variant, fill=variant)) +
    geom_area(alpha=0.8) +
    geom_hline(yintercept = 1) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x="Date, week ending",
         y="Proportion of cases with known/specified variant",
         color=NULL, 
         fill=NULL, 
         caption = glue::glue("Data up to {max(data$Kalenderwoche)}
           Source: AGES, https://www.ages.at/themen/krankheitserreger/coronavirus/sars-cov-2-varianten-in-oesterreich/
         95% ci's estimated from two-sided binomial test")) +
    ggsci::scale_color_nejm() +
    ggsci::scale_fill_nejm() +
    ggforce::facet_zoom(xlim = c(as.Date("2021-11-01"), max(data$Kalenderwoche)), ylim = c(0, determine_zoom_y_limit(data))) +
    ggpubr::theme_pubr()
}

#' Plot an area plot of the proportion of cases assigned to each variant of all
#' cases in that week. Unassigned cases will be left white/blank in the plot. 
#' 
#' @param data Preprocssed data 
#' 
#' @return Ggplot object
plot_cases_assigned_voc <- function(data) {
  data %>% 
    mutate(prop_assigned = cases_assigned/cases_all) %>%
    mutate(prop = cases/cases_all) %>%
    ggplot(aes(x=Kalenderwoche, y=prop, fill=variant)) +
    geom_area() +
    geom_hline(yintercept = 1) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(x="Date, week ending",
         y="Proportion of cases assigned to any variant",
         color=NULL, 
         fill=NULL, 
         caption = glue::glue("Data up to {max(data$Kalenderwoche)}
         Source: AGES, https://www.ages.at/themen/krankheitserreger/coronavirus/sars-cov-2-varianten-in-oesterreich/
       95% ci's estimated from two-sided binomial test")) +
    ggsci::scale_fill_nejm() +
    ggpubr::theme_pubr()
}