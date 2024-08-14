setwd("COVID19/CODE/PROJECTS/Validation/paper") ### set to repository location

source("functions_read_data.R")
source("functions_analysis.R")

library(tidyverse)
library(haven)
library(survey)
library(lubridate)
library(srvyr)
library(hrbrthemes)


#### READ CSP DATA ####

waves <- c(5,7,9:14, 16:22, 24)

demo_vars <- c("gender", "race", "age", "education", 
               "region", "state", "state_code", "urbanicity")

### reads CSP data as a list of datasets, one per wave. The variables to keep from the full datasets are 
### provided as arguments
datasets_csp <- read_CSP_data(waves, wave, StartDate, EndDate, weight, weight_state, state, state_code, id, covid, 
                              starts_with("cov_test"), cov_trust_cdc, starts_with("cov_month"),
                              gender, race, age, education, region, urbanicity, party) %>%
  lapply(function(d) mutate(d, wave_date = max(EndDate),
                            EndDate_numeric = as.numeric(as.Date(EndDate)))) 


#### CALCULATE METHOD 1 ####


# Using a hand-coded table to save the info on which waves should be used for
# which months for method 1:
out <- read_wave_month_table("Data/")
wave_month_match_table <- out[[1]]
months_per_wave <- out[[2]]

### generate the necessary variables to calculate method 1 estimates
datasets_csp_method1 <- get_datasets_method1(datasets_csp, wave_month_match_table)

### apply 'calculate_infected_method1' to each dataset, calculate the cdf and
### the cumulative error and join the dates of each block of months
method1_results <- lapply(datasets_csp_method1, calculate_infected_method1) %>%
  bind_rows() %>%
  mutate(cdf_estimate = cumsum(estimate),
         cdf_std = sqrt(cumsum(std_error^2))) %>%
  left_join(months_per_wave) %>%
  mutate(date_plot  = (month_max %m+% months(1)) - 1)



trust_downsample_iter <- function(i, fraction_keep) {
  out <- downsample_csp_trust(datasets_csp_method1, fraction_keep)
  
  print("Fraction kept:")
  print(out[[2]])
  
  datasets_trust <- out[[1]]
  
  csp_inf_trust <- lapply(datasets_trust, calculate_infected_method1) %>%
    bind_rows() %>%
    mutate(cdf_estimate = cumsum(estimate),
           cdf_std = sqrt(cumsum(std_error^2))) %>%
    left_join(months_per_wave) %>%
    mutate(date_plot  = (month_max %m+% months(1)) - 1) %>%
    mutate(iteration = i)
  
  return(list(csp_inf_trust, out[[2]]))
}


calc_inf_sims <- function(sims) {
  dfs <- lapply(sims, function(d) d[[1]])
  
  dfs %>% 
    bind_rows %>%
    group_by(wave) %>%
    summarise(inf_mean = mean(cdf_estimate),
              std_error_mean = mean(cdf_std),
              inf_sd = sd(cdf_estimate),
              inf_range = max(abs(range(cdf_estimate)-cdf_estimate)),
              date_plot = first(date_plot)) %>%
    select(-wave)  
  
}

calc_diff <- function(sims_sum) {
  sims_sum %>%
    left_join(method1_results, by = c("date_plot")) %>%
    mutate(diff = inf_mean - cdf_estimate) %>%
    select(diff, std_error_mean, std_error, date_plot)
}

n_iter <- 100

fraction_keep1 <- 
  tibble(cov_trust_cdc = c(1,2,3,4, -99), 
         fraction_keep=c(0.7,0.8,0.9,1, 1))
#sims_trust1 <- lapply(c(1:n_iter), trust_downsample_iter, fraction_keep1) 
sims_trust1 <- readRDS("Simulations/sims_trust1_inf.RDS")


fraction_keep2 <- 
  tibble(cov_trust_cdc = c(1,2,3,4, -99), 
         fraction_keep=c(0.5,0.6,0.9,1, 1))
#sims_trust2 <- lapply(c(1:n_iter), trust_downsample_iter, fraction_keep2) 
sims_trust2 <- readRDS("Simulations/sims_trust2_inf.RDS")


fraction_keep3 <- 
  tibble(cov_trust_cdc = c(1,2,3,4, -99), 
         fraction_keep=c(0.3,0.4,0.9,1, 1))
#sims_trust3 <- lapply(c(1:n_iter), trust_downsample_iter, fraction_keep3) 
sims_trust3 <- readRDS("Simulations/sims_trust3_inf.RDS")


fraction_keep4 <-
  tibble(cov_trust_cdc = c(1,2,3,4, -99),
         fraction_keep=c(0.2,0.2,1,1, 1))
#sims_trust4 <- lapply(c(1:n_iter), trust_downsample_iter, fraction_keep4) 
sims_trust4 <- readRDS("Simulations/sims_trust4_inf.RDS")


# 
# saveRDS(sims_trust1, "PROJECTS/Validation/paper/Simulations/sims_trust1_inf.RDS")
# saveRDS(sims_trust2, "PROJECTS/Validation/paper/Simulations/sims_trust2_inf.RDS")
# saveRDS(sims_trust3, "PROJECTS/Validation/paper/Simulations/sims_trust3_inf.RDS")
# saveRDS(sims_trust4, "PROJECTS/Validation/paper/Simulations/sims_trust4_inf.RDS")



diff_sims1 <- calc_inf_sims(sims_trust1) %>% 
  calc_diff() %>% 
  mutate(`Trust based downsample` = "Low")

diff_sims2 <- calc_inf_sims(sims_trust2) %>% 
  calc_diff() %>%
  mutate(`Trust based downsample` = "Mild")

diff_sims3 <- calc_inf_sims(sims_trust3) %>% 
  calc_diff() %>%
  mutate(`Trust based downsample` = "High")

diff_sims4 <- calc_inf_sims(sims_trust4) %>% 
  calc_diff() %>%
  mutate(`Trust based downsample` = "Very High")


diff_df <- bind_rows(diff_sims1, diff_sims2, diff_sims3, diff_sims4) %>%
  mutate(`Trust based downsample` =
           factor(`Trust based downsample`,
                  levels = c("Low", "Mild", "High", "Very High")))

plot_diff <- ggplot(diff_df, 
                    aes(x = date_plot, y = diff, colour = `Trust based downsample`)) + 
  geom_pointrange(aes(ymin = (diff - std_error_mean), 
                      ymax = (diff + std_error_mean)),
                  position = position_dodge2(width = 10)) +
  scale_colour_viridis_d()+
  geom_hline(yintercept = 0) +
  scale_x_date(breaks = "2 months", 
               date_labels = "%b%y") +
  xlab("Median date of survey wave") + 
  ylab("Difference between simulation and whole sample \nestimates of cumulative infection rates")+
  theme_ipsum()+
  theme(legend.position=c(0.12, 0.2), 
        legend.justification='left',
        legend.direction='vertical',
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(color='black',face='bold', size=14,
                                    margin = margin(0,8,0,0)),
        axis.title.x = element_text(color='black',face='bold', size=14,
                                    margin = margin(5,0,0,0)),
        legend.text = element_text(color='black', size=15),
        legend.title = element_text(color = "black", face = "bold",
                                    size = 15),
        panel.border = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt"))

ggsave("Plots/trust_downsample_infections.png",
       plot_diff, 
       height = 6.5, width = 9,
       bg = "white", dpi = 300)



