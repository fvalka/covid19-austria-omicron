library(targets)
library(future)
library(future.callr)
plan(callr)

# Data processing, ETL
source("R/functions/process_input_data.R")

# Models
source("R/functions/model/epidemia_models.R")
source("R/functions/model/binomial_models.R")
source("R/functions/model/multinomial_models.R")

# Vaccinations
source("R/functions/vaccination/population_immunity.R")

# Projections
source("R/functions/projections/epidemia_projections.R")

# Scenarios
source("R/functions/scenarios/epidemia_scenarios.R")

# Visualization/Plots
source("R/functions/plots/proportion_plots.R")
source("R/functions/plots/twitter_plots.R")
source("R/functions/plots/binomial_model_plots.R")
source("R/functions/plots/multinomial_model_plots.R")
source("R/functions/plots/epidemia_plots.R")

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("readr", "brms", "dplyr", "ggplot2", "epidemia", "here", "emmeans"))
list(
  tar_target(
    raw_data_file,
    "https://www.ages.at/fileadmin/AGES2015/Themen/Krankheitserreger_Dateien/Coronavirus/Varianten_ab_Mai/Varianten_nach_KWs_2022-01-04.csv",
    format = "url"
  ),
  tar_target(
    raw_ems_case_file,
    "data/CovidFaelle_Timeline.csv",
    format = "file"
  ),
  tar_target(
    raw_vaccination_data_file,
    "data/COVID19_vaccination_doses_timeline.csv",
    format = "file"
  ),
  tar_target(
    raw_data,
    read_csv2(raw_data_file,  locale = readr::locale(encoding = "latin1")) %>%
      mutate(filename=basename(raw_data_file))
  ),
  tar_target(
    data,
    process_variant_data(raw_data)
  ),
  tar_target(
    cases_ems_austria,
    process_ems_cases(raw_ems_case_file)
  ),
  tar_target(
    vaccination_data,
    read_csv2(raw_vaccination_data_file)
  ),
  tar_target(
    bnt_bnt_ve_file,
    "data/ve_data/omicron/bnt162b2/dose_2.csv",
    format = "file"
  ),
  tar_target(
    bnt_bnt_bnt_ve_file,
    "data/ve_data/omicron/bnt162b2/dose_3.csv",
    format = "file"
  ),
  tar_target(
    bnt_bnt_ve,
    estimate_ve_from_study(bnt_bnt_ve_file, length(vaccination_data$date), start_day = 14)
    
  ),
  tar_target(
    bnt_bnt_bnt_ve,
    estimate_ve_from_study(bnt_bnt_bnt_ve_file, length(vaccination_data$date))
  ),
  tar_target(
    pop_immunity_estimate,
    estimate_population_immunity(vaccination_data, 
                                 ve_dose_1 = estimate_population_immunity_assumptions$zero_ve,
                                 ve_dose_2 = bnt_bnt_ve, 
                                 ve_dose_3 = bnt_bnt_bnt_ve)
  ),
  tar_target(
    epidemia_fit_omicron, 
    epidemia_generation_time_scenario(data = data,
                                      variants_investigated = "B.1.1.529 (Omikron)",
                                      min_cases = 1,
                                      gen_time_mean = 3.13,
                                      gen_time_sd = 2.22,
                                      i2o_cases = epidemia_param_data$i2o_omicron_si_shortened,
                                      pop_immunity = pop_immunity_estimate)
  ),
  tar_target(
    epidemia_fit_delta, 
    epidemia_fit_variant(variants_investigated = "B.1.617.2 (Delta)",
                         variants_all = data,
                         ostart = as.Date("2021-05-01"),
                         generation_time = EpiEstim::discr_si(0:99, mu = 4.6, sigma = 3.1)/sum(EpiEstim::discr_si(0:99, mu = 4.6, sigma = 3.1)))
  ),
  tar_target(
    epidemia_projection_central, 
    epidemia_project(epidemia_fit_delta,
                     epidemia_fit_omicron,
                     extend_days = 30)
  ),
  tar_target(
    epidemia_projection_central_infections,
    epidemia_project(epidemia_fit_delta,
                     epidemia_fit_omicron,
                     extend_days = 30,
                     predictor = epidemia_projection_infections_generator)
  ),
  tar_target(
    scenarios_generation_times_mean_omicron,
    c(3.13, 4.6, 3.13, 2.84, 2.84, 2.22)
  ),
  tar_target(
    scenarios_generation_times_sd_omicron,
    c(2.22, 3.1, 1.6, 1.6, 2.2, 1.57)
  ),
  tar_target(
    epidemia_fit_generation_time_scenarios,
    epidemia_generation_time_scenario(data=data,
      gen_time_mean = scenarios_generation_times_mean_omicron, 
      gen_time_sd = scenarios_generation_times_sd_omicron,
      min_cases = 1,
      i2o_cases = epidemia_param_data$i2o_omicron_si_shortened,
      pop_immunity = pop_immunity_estimate),
    pattern = map(scenarios_generation_times_mean_omicron, scenarios_generation_times_sd_omicron),
    iteration = "list"
  ),
  tar_target(
    binomial_fit,
    fit_binomial_glm(data)
  ),
  tar_target(
    quasibinomial_fit,
    fit_quasibinomial_glm(data)
  ),
  tar_target(
    multinomial_fit,
    fit_bayesian_multinomial_model(data)
  ),
  tar_target(
    multinomial_growth_rate,
    extract_bayesian_multinomial_growth_rate(multinomial_fit)
  ),
  tar_target(
    epidemia_generation_time_scenarios_plot,
    plot_epidemia_generation_time_scenarios(data, epidemia_fit_generation_time_scenarios)
  ),
  tar_target(
    plot_quasibinomial_fit_ci_zoom,
    plot_glm_with_ci(data, quasibinomial_fit)
  ),
  tar_target(
    plot_quasibinomial_log_odds,
    plot_glm_log_odds(data, quasibinomial_fit)
  ),
  tar_target(
    plot_omicron_proportion_assigned,
    plot_extraploted_variant_proportion(data)
  ),
  tar_target(
    plot_variants_area_zoomed,
    plot_variants_area(data)
  ),
  tar_target(
    plot_assigned,
    plot_cases_assigned_voc(data)
  ),
  tar_target(
    plot_multinomial_fit_zoomed,
    plot_multinomial_fit(data, multinomial_fit)
  ),
  tar_target(
    plot_multinomial_projection, 
    plot_multinomial_fit(data, multinomial_fit, extend_days = 60, show_zoom = F) + 
      coord_cartesian(xlim = c(as.Date("2021-12-01"), as.Date("2022-02-01")))
  ),
  tar_target(
    plot_multinomial_growth_advantage_posterior,
    plot_multinomial_growth_advantage(multinomial_fit)
  ),
  tar_target(
    plot_multinomial_half,
    plot_probablity_more_than_half(data, multinomial_fit)
  ),
  tar_target(
    file_raw_data_export,
    export_raw_data(raw_data),
    format = "file"
  ),
  tar_target(
    file_extraploted_variant_plot,
    export_plot_with_title(data, plot_omicron_proportion_assigned),
    format = "file"
  ),
  tar_target(
    file_all_assigned_plot,
    export_plot_with_title(data, plot_assigned, "austria-variants-assigned.png"),
    format = "file"
  ),
  tar_target(
    file_variants_area_plot,
    export_plot_with_title(data, plot_variants_area_zoomed, "austria-variants-area.png"),
    format = "file"
  ),
  tar_target(
    file_multinomial_zoomed,
    export_plot_with_title(data, plot_multinomial_fit_zoomed, "austria-variants-multinomial-zoom.png", subtitle="Multinomial GLMM"),
    format = "file"
  ),
  tar_target(
    file_multinomial_projection,
    export_plot_with_title(data, plot_multinomial_projection, "austria-variants-multinomial-projection.png", subtitle="Multinomial GLMM"),
    format = "file"
  ),
  tar_target(
    file_plot_quasibinomial_fit_ci_zoom,
    export_plot_with_title(data, plot_quasibinomial_fit_ci_zoom, "austria-variants-binomial-glm.png",
                           subtitle = "Overdispersed binomial GLM fit\nShaded area showing 95% CI of quasibinomial model")
  ),
  tar_target(
    file_plot_quasibinomial_fit_log_odds,
    export_plot_with_title(data, plot_quasibinomial_log_odds, "austria-variants-binomial-glm-log-odds.png",
                           subtitle = "Overdispersed binomial GLM fit\nShaded area showing 95% CI of quasibinomial model")
  )
)