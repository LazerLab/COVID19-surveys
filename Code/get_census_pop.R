library(tidycensus)
library(tidyverse)
setwd("COVID19/CODE")


state <- get_decennial(geography = "state", 
              variables = "DP1_0021C", 
              year = 2020,
              sumfile = "dp")

nat <- get_decennial(geography = "us", 
              variables = "DP1_0021C", 
              year = 2020,
              sumfile = "dp")

df <- bind_rows(nat, state) %>%
  select(NAME, value) 

write_csv(df, file = "Data/18_pop.csv")



state_total <- get_decennial(geography = "state", 
                       variables = "DP1_0001C", 
                       year = 2020,
                       sumfile = "dp")

nat_total <- get_decennial(geography = "us", 
              variables = "DP1_0001C", 
              year = 2020,
              sumfile = "dp")


df <- bind_rows(nat_total, state_total) %>%
  select(NAME, value) 

write_csv(df, file = "Data/total_pop.csv")
