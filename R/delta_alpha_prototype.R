library(rvest)
library(dplyr)
library(EpiEstim)
library(ggplot2)

ages_varianten <- read_html("https://www.ages.at/themen/krankheitserreger/coronavirus/sars-cov-2-varianten-in-oesterreich/")

tables <- ages_varianten %>%
  html_table()


variant_cases <- tables[[1]] %>% 
  filter(!is.na(Kalenderwoche)) %>%
  reshape2::melt(id.vars=c("Kalenderwoche", "Fälle gesamt"),
                 variable.name="variant",
                 value.name="cases") 

variant_cases_assigned <- variant_cases %>%
  group_by(Kalenderwoche) %>%
  summarise(cases_assigned=sum(cases))

variant_cases %>%
  left_join(variant_cases_assigned) %>%
  mutate(proportion=cases/cases_assigned) %>%
  mutate(cases_extrapolated = proportion * `Fälle gesamt`) %>%
  ggplot(aes(x=Kalenderwoche, y=cases_extrapolated, color=variant)) +
  geom_line()

variant_cases %>%
  left_join(variant_cases_assigned) %>%
  mutate(proportion=cases/cases_assigned) %>%
  mutate(cases_extrapolated = proportion * `Fälle gesamt`) %>%
  mutate(proportion_assigned=cases_assigned/`Fälle gesamt`) %>%
  ggplot(aes(x=Kalenderwoche, y=proportion, color=variant)) +
  geom_line() +
  scale_y_continuous(labels = scales::percent) +
  ggtitle("Anteil der Fälle ")


variant_cases %>%
  left_join(variant_cases_assigned) %>%
  mutate(proportion_assigned=cases_assigned/`Fälle gesamt`) %>%
  ggplot(aes(x=Kalenderwoche, y=proportion_assigned)) +
  geom_line() +
  scale_y_continuous(labels = scales::percent) +
  ggtitle("Anteil der VOC Fälle an den Gesamtfällen", subtitle = "Österreich") +
  labs(y="Anteil VOCs an Gesamtfällen", caption = "Quelle: SARS-CoV-2-Varianten in Österreich, AGES")
