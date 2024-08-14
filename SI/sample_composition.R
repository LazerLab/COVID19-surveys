setwd("COVID19/CODE/PROJECTS/Validation/paper") ### set to repository location

source("functions_read_data.R")
source("functions_analysis.R")

library(tidyverse)
library(survey)
library(srvyr)

path <- "PROJECTS/Validation/paper/"


waves <- c(5,7,9:14, 16:26)

### reads CSP data as a list of datasets, one per wave. The variables to keep from the full datasets are 
### provided as arguments
datasets_csp <- read_CSP_data(waves, wave, StartDate, EndDate, weight, weight_state, state, state_code, id, covid, 
                              starts_with("cov_test"), age, education, urbanicity, party, region, service,
                              starts_with("cov_month"), gender, age_cat_4, age_cat_6, income, race, 
                              education_cat, urban_type, party3, employ) 

datasets_csp <- lapply(datasets_csp, function(d) weight_vars(d))
datasets_csp <- generate_var_method2(datasets_csp)


pop <- weight_pop()   
pop$edu_dist$Freq = c(sum(pop$edu_age_dist$Freq[pop$edu_age_dist$edu_age %in% 
                                          c("age_18_24_no_hs",
                                            "age_25_44_no_hs",
                                            "age_45_64_no_hs",
                                            "age_65plus_no_hs")]),
                  sum(pop$edu_age_dist$Freq[pop$edu_age_dist$edu_age %in% 
                                              c("age_18_24_hs",
                                                "age_25_44_hs",
                                                "age_45_64_hs",
                                                "age_65plus_hs")]),
                  sum(pop$edu_age_dist$Freq[pop$edu_age_dist$edu_age %in% 
                                              c("age_18_24_some_col",
                                                "age_25_44_some_col",
                                                "age_45_64_some_col",
                                                "age_65plus_some_col")]),
                  sum(pop$edu_age_dist$Freq[pop$edu_age_dist$edu_age %in% 
                                              c("age_18_24_ba",
                                                "age_25_44_ba",
                                                "age_45_64_ba",
                                                "age_65plus_ba")]),
                  sum(pop$edu_age_dist$Freq[pop$edu_age_dist$edu_age %in% 
                                              c("age_18_24_grad",
                                                "age_25_44_grad",
                                                "age_45_64_grad",
                                                "age_65plus_grad")]))


dat_fin <- readRDS("../../../../DATA/RDS/CSP_W26.RDS") %>%
  weight_vars()

dat_fin <- generate_var_method2(list(dat_fin))[[1]]


race <- lapply(datasets_csp, function(d) {
  d %>% as_survey_design(weight = weight) %>% 
    group_by(race_cat) %>%
    summarise(Sample = n()/nrow(.),
              Sample_wtd = survey_total()/nrow(.),
              Wave = first(wave)) %>%
    select(-Sample_wtd_se)
}) %>% 
  bind_rows() %>%
  group_by(race_cat) %>%
  summarise(Sample_mean = round(100*mean(Sample)),
            Sample_sd = round(100*sd(Sample)), 
            Sample_wtd_mean = round(100*mean(Sample_wtd))) %>%
  left_join(pop$race_dist) %>% 
  rename(Population = Freq) %>%
  mutate(Population = round(100*Population)) %>%
  left_join(
    dat_fin %>% 
      as_survey_design(weight = weight) %>%
      group_by(race_cat) %>%
      summarise(vax_rate = 100*survey_mean(vaccine_1, na.rm = T),
                inf_rate = 100*survey_mean(count_cases, na.rm = T)) %>%
      mutate(vax_rate = str_c(round(vax_rate), " (", round(vax_rate_se*1.96), ")"),
             inf_rate = str_c(round(inf_rate), " (", round(inf_rate_se*1.96), ")")) %>%
      select(race_cat, vax_rate, inf_rate)
  )


education <- lapply(datasets_csp, function(d) {
  d %>% as_survey_design(weight = weight) %>% 
    group_by(edu_cat) %>%
    summarise(Sample = n()/nrow(.),
              Sample_wtd = survey_total()/nrow(.),
              Wave = first(wave)) %>%
    select(-Sample_wtd_se)
}) %>% 
  bind_rows() %>%
  group_by(edu_cat) %>%
  summarise(Sample_mean = round(100*mean(Sample)),
            Sample_sd = round(100*sd(Sample)), 
            Sample_wtd_mean = round(100*mean(Sample_wtd))) %>%
  left_join(pop$edu_dist) %>% 
  rename(Population = Freq) %>%
  mutate(Population = round(100*Population)) %>%
  slice(4,3,5,1,2) %>%
  left_join(
    dat_fin %>% 
      as_survey_design(weight = weight) %>%
      group_by(edu_cat) %>%
      summarise(vax_rate = 100*survey_mean(vaccine_1, na.rm = T),
                inf_rate = 100*survey_mean(count_cases, na.rm = T)) %>%
      mutate(vax_rate = str_c(round(vax_rate), " (", round(vax_rate_se*1.96), ")"),
             inf_rate = str_c(round(inf_rate), " (", round(inf_rate_se*1.96), ")")) %>%
      select(edu_cat, vax_rate, inf_rate)
  )


urbanicity <- lapply(datasets_csp, function(d) {
  d %>% as_survey_design(weight = weight) %>% 
    group_by(urban_cat) %>%
    summarise(Sample = n()/nrow(.),
              Sample_wtd = survey_total()/nrow(.),
              Wave = first(wave)) %>%
    select(-Sample_wtd_se)
}) %>% 
  bind_rows() %>%
  group_by(urban_cat) %>%
  summarise(Sample_mean = round(100*mean(Sample)),
            Sample_sd = round(100*sd(Sample)), 
            Sample_wtd_mean = round(100*mean(Sample_wtd))) %>%
  left_join(pop$urban_dist) %>% 
  rename(Population = Freq) %>%
  mutate(Population = round(100*Population)) %>%
  left_join(
    dat_fin %>% 
      as_survey_design(weight = weight) %>%
      group_by(urban_cat) %>%
      summarise(vax_rate = 100*survey_mean(vaccine_1, na.rm = T),
                inf_rate = 100*survey_mean(count_cases, na.rm = T)) %>%
      mutate(vax_rate = str_c(round(vax_rate), " (", round(vax_rate_se*1.96), ")"),
             inf_rate = str_c(round(inf_rate), " (", round(inf_rate_se*1.96), ")")) %>%
      select(urban_cat, vax_rate, inf_rate)
  )
  

gender <- lapply(datasets_csp, function(d) {
  d %>% as_survey_design(weight = weight) %>% 
    group_by(gender) %>%
    summarise(Sample = n()/nrow(.),
              Sample_wtd = survey_total()/nrow(.),
              Wave = first(wave)) %>%
    select(-Sample_wtd_se)
}) %>% 
  bind_rows() %>%
  group_by(gender) %>%
  summarise(Sample_mean = round(100*mean(Sample)),
            Sample_sd = round(100*sd(Sample)), 
            Sample_wtd_mean = round(100*mean(Sample_wtd))) %>%
  left_join(pop$gender_dist) %>% 
  rename(Population = Freq) %>%
  mutate(Population = round(100*Population)) %>%
  left_join(
    dat_fin %>% 
      as_survey_design(weight = weight) %>%
      group_by(gender) %>%
      summarise(vax_rate = 100*survey_mean(vaccine_1, na.rm = T),
                inf_rate = 100*survey_mean(count_cases, na.rm = T)) %>%
      mutate(vax_rate = str_c(round(vax_rate), " (", round(vax_rate_se*1.96), ")"),
             inf_rate = str_c(round(inf_rate), " (", round(inf_rate_se*1.96), ")")) %>%
      select(gender, vax_rate, inf_rate)
  )
  

age <- lapply(datasets_csp, function(d) {
  d %>% as_survey_design(weight = weight) %>% 
    group_by(age_cat) %>%
    summarise(Sample = n()/nrow(.),
              Sample_wtd = survey_total()/nrow(.),
              Wave = first(wave)) %>%
    select(-Sample_wtd_se)
}) %>% 
  bind_rows() %>%
  group_by(age_cat) %>%
  summarise(Sample_mean = round(100*mean(Sample)),
            Sample_sd = round(100*sd(Sample)), 
            Sample_wtd_mean = round(100*mean(Sample_wtd))) %>%
  left_join(pop$age_dist) %>% 
  rename(Population = Freq) %>%
  mutate(Population = round(100*Population)) %>%
  left_join(
    dat_fin %>% 
      as_survey_design(weight = weight) %>%
      group_by(age_cat) %>%
      summarise(vax_rate = 100*survey_mean(vaccine_1, na.rm = T),
                inf_rate = 100*survey_mean(count_cases, na.rm = T)) %>%
      mutate(vax_rate = str_c(round(vax_rate), " (", round(vax_rate_se*1.96), ")"),
             inf_rate = str_c(round(inf_rate), " (", round(inf_rate_se*1.96), ")")) %>%
      select(age_cat, vax_rate, inf_rate)
  )


region <- lapply(datasets_csp, function(d) {
  d %>% as_survey_design(weight = weight) %>% 
    group_by(region) %>%
    summarise(Sample = n()/nrow(.),
              Sample_wtd = survey_total()/nrow(.),
              Wave = first(wave)) %>%
    select(-Sample_wtd_se)
}) %>% 
  bind_rows() %>%
  group_by(region) %>%
  summarise(Sample_mean = round(100*mean(Sample)),
            Sample_sd = round(100*sd(Sample)), 
            Sample_wtd_mean = round(100*mean(Sample_wtd))) %>%
  left_join(pop$region_dist) %>% 
  rename(Population = Freq) %>%
  mutate(Population = round(100*Population)) %>%
  left_join(
    dat_fin %>% 
      as_survey_design(weight = weight) %>%
      group_by(region) %>%
      summarise(vax_rate = 100*survey_mean(vaccine_1, na.rm = T),
                inf_rate = 100*survey_mean(count_cases, na.rm = T)) %>%
      mutate(vax_rate = str_c(round(vax_rate), " (", round(vax_rate_se*1.96), ")"),
             inf_rate = str_c(round(inf_rate), " (", round(inf_rate_se*1.96), ")")) %>%
      select(region, vax_rate, inf_rate)
  )


party <- lapply(datasets_csp, function(d) {
  d %>% as_survey_design(weight = weight) %>% 
    group_by(party_pid) %>%
    summarise(Sample = n()/nrow(.),
              Sample_wtd = survey_total()/nrow(.),
              Wave = first(wave)) %>%
    select(-Sample_wtd_se)
}) %>% 
  bind_rows() %>%
  group_by(party_pid) %>%
  summarise(Sample_mean = round(100*mean(Sample)),
            Sample_sd = round(100*sd(Sample)), 
            Sample_wtd_mean = round(100*mean(Sample_wtd))) %>%
  left_join(pop$party_pid_dist, by = c("party_pid" = "party_w")) %>% 
  rename(Population = Freq) %>%
  mutate(Population = round(100*Population)) %>%
  left_join(
    dat_fin %>% 
      as_survey_design(weight = weight) %>%
      group_by(party_pid) %>%
      summarise(vax_rate = 100*survey_mean(vaccine_1, na.rm = T),
                inf_rate = 100*survey_mean(count_cases, na.rm = T)) %>%
      mutate(vax_rate = str_c(round(vax_rate), " (", round(vax_rate_se*1.96), ")"),
             inf_rate = str_c(round(inf_rate), " (", round(inf_rate_se*1.96), ")")) %>%
      select(party_pid, vax_rate, inf_rate)
  )


