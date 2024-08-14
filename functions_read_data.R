#source("DATA PROCESS/csp_functions.R")
library(tidyverse)
library(lubridate)
library(haven)



read_CSP_data <- function(waves, ...){
  subset = missing(...)
  
  path = "../../../../DATA/RDS/CSP_W"
  
  datasets <- lapply(waves, function(i) {
    d <- readRDS(paste0(path, i, ".RDS"))
    
    if (d$wave[1] == 14) d <- d %>% filter(StartDate < as.Date("2021-01-11"))
    if (d$wave[1] == 19) d <- d %>% filter(StartDate < as.Date("2021-09-28"))
    if (d$wave[1] == 21) d <- d %>% filter(StartDate < as.Date("2022-01-26"))
    if (d$wave[1] == 23) d <- d %>% filter(StartDate < as.Date("2022-07-06"))
    
    if(!subset) {
      d <- select(d, ...)
    }
    
    return(d)
  }) 
  
  return(datasets)
}

read_preprocess_pop_data <- function(census_pop_path = "Data/18_pop.csv",
                                     state_codes_path = "https://raw.githubusercontent.com/jasonong/List-of-US-States/master/states.csv") {
  
  state_codes <- read_csv(state_codes_path,
                          col_names = TRUE) %>%
    rename(state_code = Abbreviation)
  
  pop18 <- read_csv(census_pop_path) %>%
    rename(State = NAME,
           Population = value) %>%
    left_join(state_codes)
  pop18[pop18$State == "United States", "state_code"] <- "US"
  
  return(pop18)
  
}

read_preprocess_CDC_vax_data <- function(cdc_vax_data_path = "../COVID-19_Vaccination_Age_and_Sex_Trends_in_the_United_States__National_and_Jurisdictional.csv"){
  

  
  age_cats <- c("Ages_<5yrs", "Ages_5-11_yrs", "Ages_12-17_yrs", "Ages_18-24_yrs",
                "Ages_25-49_yrs", "Ages_50-64_yrs", "Ages_65+_yrs")
  
  pop18 <- read_preprocess_pop_data()
  
  cdc <- read_csv(cdc_vax_data_path) %>%
    select(Date, Location, Demographic_Category, Administered_Dose1) %>%
    mutate(Date = as.Date(mdy_hms(Date))) %>%
    filter(Demographic_Category %in% age_cats) %>%
    mutate(under_18_flag = as.numeric(Demographic_Category %in%
                                        c("Ages_<5yrs", "Ages_5-11_yrs", "Ages_12-17_yrs")),
           over_18_flag = as.numeric(Demographic_Category %in%
                                       c("Ages_18-24_yrs", "Ages_25-49_yrs", "Ages_50-64_yrs", 
                                         "Ages_65+_yrs"))) %>%
    group_by(Date, Location) %>%
    summarize(n_vax = sum(Administered_Dose1 * (over_18_flag), na.rm = T)) %>%
    left_join(pop18, by = c("Location" = "state_code")) %>%
    filter(!is.na(State)) %>%
    mutate(pct_vax = 100 * n_vax/Population) %>%
    ungroup()
  
  return(cdc)
} 
    

read_preprocess_IPSOS_data <- function(ipsos_path = "../../../../DATA/AUX DATA/Ipsos/",
                                       ipsos_waves_dates = "Data/Ipsos_dates.csv",
                                       date_return = "mid") {
  
  
  dfs <- lapply(list.files(ipsos_path), 
                function(f) {
                  if(tools::file_ext(f) == "por"){
                    return(read_por(str_c(ipsos_path, f)))
                  } 
                  if(tools::file_ext(f) == "sav"){
                    return(read_sav(str_c(ipsos_path, f)))
                  }
                })
  
  
  dfs[[length(dfs)]]$WAVE <- 71
  dfs[[length(dfs)]]$WT_FINAL <- dfs[[length(dfs)]]$WEIGHTS 
  
  if(date_return == "mid"){
    dates <- read_csv(ipsos_waves_dates) %>%
      mutate(Date1 = mdy(Date1),
             Date2 = mdy(Date2)) %>%
      mutate(Date = Date1 + (Date2 - Date1)/2) %>%
      rename(WAVE = Wave) %>%
      select(WAVE, Date)
  }
  
  if(date_return == "last"){
    dates <- read_csv(ipsos_waves_dates) %>%
      mutate(Date = mdy(Date2)) %>%
      rename(WAVE = Wave) %>%
      select(WAVE, Date)
  }

  
  return(list(dfs, dates))
}

read_infections_national_admin_data <- function(nyt_path = "https://raw.githubusercontent.com/nytimes/covid-19-data/master/us.csv",
                                       write_path = "../nyt_us_data.csv",
                                       census_path = "Data/total_pop.csv") {
  
  us_cumulative_nyt <- fread(nyt_path)
  ### save copy of data for security
  write_csv(us_cumulative_nyt, write_path)
  
  us_cumulative_nyt <- data.frame(us_cumulative_nyt)
  us_cumulative_nyt$date <- as.Date(us_cumulative_nyt$date)
  
  census_pop <- read_csv(census_path)
  us_pop_2021 <- census_pop %>% filter(NAME == "United States") %>% pull(value)
  
  us_cumulative_nyt <- us_cumulative_nyt %>% 
    mutate(pct_cumulative_cases = cases*100/us_pop_2021)
  
  return(us_cumulative_nyt)
  
}


read_infections_state_admin_data <- function(nyt_path = "https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv",
                                                write_path = "../nyt_states_data.csv",
                                                census_path = "Data/total_pop.csv") {

  states_cumulative_nyt <- fread(nyt_path)
  ### save copy of data for security
  write_csv(states_cumulative_nyt, write_path)
  
  states_cumulative_nyt <- data.frame(states_cumulative_nyt)
  states_cumulative_nyt$date <- as.Date(states_cumulative_nyt$date)
  
  census_pop <- read_csv(census_path)
  
  states_cumulative_nyt <- states_cumulative_nyt %>%
    left_join(census_pop, by = c("state" = "NAME")) %>%
    mutate(pct_cumulative_cases = cases*100/value) %>%
    filter(!is.na(value), !state == "Puerto Rico")
  
  return(states_cumulative_nyt)
}
