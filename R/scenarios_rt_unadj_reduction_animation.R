epidemia_fit_omicron <- tar_read("epidemia_fit_omicron")
extended_data <- epidemia_extend_cases(epidemia_fit_omicron)

data_for_scenario <- function(t_switch=as.Date("2022-01-08") , rt_adj_new = 0.33) {
  newdata <- extended_data %>%
    mutate(rt_adj = case_when(date >= t_switch ~ rt_adj_new,
                              T ~ rt_adj))
}

plot_scenario <- function(t_switch, hide_y_axis=F, estimate_rt_unadj_reduction=F) {
  newdata_scenario = data_for_scenario(t_switch)
  
  if(estimate_rt_unadj_reduction) {
    rt_unadj <- posterior_rt(epidemia_fit_omicron$fit, newdata = newdata_scenario, adjusted = FALSE)
    
    switch_idx <- which(rt_unadj$time == t_switch)
    message(median(rt_unadj$draws[, switch_idx + 1]/rt_unadj$draws[, switch_idx - 1]))
  }
  
  plot_rt_scenario <- function() {
    plot_rt(epidemia_fit_omicron$fit, newdata=newdata_scenario, draws=3000) + 
      geom_vline(xintercept = t_switch, linetype="dashed", alpha=0.6) +
      coord_cartesian(xlim=c(as.Date("2022-01-01"), NA), ylim=c(0, 4)) +
      ggtitle("Effective reproduction number") +
      labs(y="Effective reproduction number") 
  }
  
  plot_cases_scenario <- function() {
    plot_obs(epidemia_fit_omicron$fit, newdata=newdata_scenario, type="cases", by_100k = T, draws=3000) +
      annotate("rect", xmin = as.Date("2021-12-14"), 
               xmax = as.Date("2022-03-01"), 
               ymin = 1100,
               ymax = 11000,
               fill = "black",
               alpha = 0.1) +
      ggtext::geom_richtext(x=as.Date("2022-01-05"), y=5550, label="Might exceed testing capacity", angle = 90, size=3, alpha=0.7) +
      geom_vline(xintercept = t_switch, linetype="dashed", alpha=0.6) +
      coord_cartesian(xlim=c(as.Date("2022-01-01"), NA), ylim=c(0, 1e4)) +
      ggtitle("Weekly cases per 100k") +
      labs(y="Weekly cases per 100k")
  }
  
  plot_infections_scenario <- function() {
    plot_infections(epidemia_fit_omicron$fit, newdata=newdata_scenario, by_100k = T, draws=3000) +
      geom_vline(xintercept = t_switch, linetype="dashed", alpha=0.6) +
      coord_cartesian(xlim=c(as.Date("2022-01-01"), NA), ylim=c(0,3e3)) +
      ggtitle("Daily infections per 100k") +
      labs(y="Daily infections per 100k")
  }
  
  p_result <- (plot_rt_scenario() |
      plot_infections_scenario() | 
      plot_cases_scenario()) +
    plot_layout(heights = c(2, 4, 4)) &
    theme(legend.position = "none",
          plot.margin=unit(c(6,2,4,2), "mm"),
          plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(face = "bold", size = 18))
  
  if(hide_y_axis) {
    p_result <- p_result &
      theme(axis.text.y=element_blank(), 
            axis.ticks.y=element_blank())  
  }
  
  p_result
}

plot_single_date <- function(plot_date) {

  plot_scenario(plot_date) +
    plot_annotation(title = "Assumed 50% reduction in reproduction number before immunity on:",
                    subtitle = plot_date,
                    caption = "Black line: Median. Colored shaded areas: Credible intervals: 30%, 60% and 90%
  Grey shaded area: Actual case numbers likely to be limited by testing capacity, higher estimates likely never to be ascertainable") 
  
  ggsave(glue::glue("output/reduction_anim/{plot_date}.png"), width = 10, height = 5)
}

seq(as.Date("2022-01-03"), as.Date("2022-02-15"), by=1) %>%
  purrr::map(~ plot_single_date(.))

