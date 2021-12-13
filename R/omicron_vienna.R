library(here)
library(dplyr)
library(readr)
library(ggplot2)
library(lubridate)

omicron_cases <- read_csv(here("data/omicron_cases_vienna.csv")) %>%
  mutate(date=lubridate::dmy(date),
         cases_new=cases-lag(cases,1))

all_cases <- read_csv2(here("data/CovidFaelle_Timeline.csv")) %>%
  mutate(Time=as.Date(dmy_hms(Time))) %>%
  rename(date=Time) %>%
  filter(Bundesland=="Wien")

omicron_cases %>%
  ggplot(aes(x=date, y=cases_new)) +
  geom_bar(stat="identity")

omicron_cases %>%
  left_join(all_cases) %>%
  mutate(prop=cases_new/AnzahlFaelle) %>%
  ggplot(aes(x=date, y=prop)) +
  geom_point()
