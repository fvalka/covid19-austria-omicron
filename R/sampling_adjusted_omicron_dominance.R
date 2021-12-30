library(targets)
library(dplyr)
library(ggplot2)
library(here)
library(brms)

source(here("R/functions/model/multinomial_models.R"))
source(here("R/functions/process_input_data.R"))

process_cumulative_prop <- function(data, fit, sampling_adjust) { 
  draw_samples_from_multinomial_fit(data, fit, ndraws = 2000, extend_days = 60) %>%
  filter(variant_full_name == "B.1.1.529 (Omikron)") %>%
  filter(prop >= 0.5) %>%
  group_by(.draw) %>%
  summarise(Kalenderwoche = min(Kalenderwoche)) %>%
  mutate(sampling_adjust=sampling_adjust)
}

tar_load("data")

fit_quarter_adj <- fit_bayesian_multinomial_model(data_adjust_sampling(data, 0.25))
fit_half_adj <- fit_bayesian_multinomial_model(data_adjust_sampling(data, 0.5))
fit_three_quarter_adj <- fit_bayesian_multinomial_model(data_adjust_sampling(data, 0.75))
fit_fully_adj <- fit_bayesian_multinomial_model(data_adjust_sampling(data, 1))

process_cumulative_prop(data, tar_read("multinomial_fit"), 0) %>%
  rbind(process_cumulative_prop(data, fit_quarter_adj, 0.25)) %>%
  rbind(process_cumulative_prop(data, fit_half_adj, 0.5)) %>%
  rbind(process_cumulative_prop(data, fit_three_quarter_adj, 0.75)) %>%
  rbind(process_cumulative_prop(data, fit_fully_adj, 1)) %>%
  ggplot(aes(x=Kalenderwoche, group=sampling_adjust, color=sampling_adjust)) +
  stat_ecdf(geom = "line") +
  labs(x="Date, week ending",
       y="Probability Omicron Dominance",
       color="Sampling Adjustment") +
  scale_x_date(date_labels =  "%d %b %Y") +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(xlim = c(as.Date("2021-12-7"), as.Date("2022-01-21"))) +
  ggpubr::theme_pubr() +
  ggtitle("Probability of Omicron dominance (more than half Omicron cases)")


ggsave("output/austria-more-than-half-omicron-sample_adj.png", width = 10, height = 5)
