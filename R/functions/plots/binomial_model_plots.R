#' Plot the classical binomial [stats::glm()] fit with confidence interval.
#' 
#' Uses [emmeans::emmeans()] for the response projection CI. 
#' 
#' @param data Preprocessed variant data
#' @param fit GLM fit object
#' 
#' @return Ggplot object for the plot configuration
plot_glm_with_ci <- function(data, fit) {
  dates <- seq(min(fit$data$date)-30, max(fit$data$date), 1)
  
  projection <- emmeans::emmeans(fit, ~date, at=list(date=dates), type="response") %>%
    as.data.frame()
  
  data %>% 
    filter(Kalenderwoche >= as.Date("2021-09-01")) %>%
    ggplot(aes(x=Kalenderwoche, y=proportion)) +
    geom_point(aes(color=variant)) +
    geom_line(data=projection, aes(x=date, y=prob, color="B.1.1.529 (Omikron)")) +
    geom_ribbon(data=projection, aes(x=date, y=prob, ymin=asymp.LCL, ymax=asymp.UCL, fill="B.1.1.529 (Omikron)"), alpha=0.4) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggsci::scale_color_nejm() +
    ggsci::scale_fill_nejm(limits=unique(data$variant)) +
    labs(x="Date, week ending",
         y="Proportion of cases with known/specified variant",
         color=NULL, 
         fill=NULL) +
    ggforce::facet_zoom(xlim = c(as.Date("2021-11-01"), max(data$Kalenderwoche)), ylim = c(0, determine_zoom_y_limit(data))) +
    ggpubr::theme_pubr()

}

#' Plot a [stats::glm()] fit on the log-odds (logit) scale. 
#' 
#' Uses [emmeans::emmeans()] for projection. 
#' 
#' @param data Preprocssed variant data
#' @param fit GLM fit object
#' 
#' @return Ggplot object of the log-odds plot with CI's
plot_glm_log_odds <- function(data, fit) {
  dates <- seq(min(fit$data$date), max(fit$data$date), 1)
  
  projection <- emmeans::emmeans(fit, ~date, at=list(date=dates)) %>%
    as.data.frame()
  
  data %>% 
    filter(variant=="B.1.1.529 (Omikron)") %>%
    filter(cases > 0) %>%
    mutate(log_odds = log(proportion/(1-proportion))) %>%
    ggplot(aes(x=Kalenderwoche, y=log_odds)) +
    geom_point(aes(color=variant, size=cases)) +
    geom_line(data=projection, aes(x=date, y=emmean, color="B.1.1.529 (Omikron)")) +
    geom_ribbon(data=projection, aes(x=date, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, fill="B.1.1.529 (Omikron)"), alpha=0.4) +
    ggsci::scale_color_nejm(limits=unique(data$variant)) +
    ggsci::scale_fill_nejm(limits=unique(data$variant)) +
    scale_size_area(limits = c(0, 1e3),
                    oob=scales::squish,
                    guide="none") +
    labs(x="Date, week ending",
         y="Log odds of Omicron variant among all VOC cases",
         color=NULL, 
         fill=NULL) +
    ggpubr::theme_pubr()
}