setwd("COVID19/CODE/PROJECTS/Validation/paper") ### set to repository location

source("functions_read_data.R")
source("functions_analysis.R")

library(tidyverse)
library(haven)
library(survey)
library(lubridate)
library(srvyr)
library(hrbrthemes)


### read CSP data ####


waves <- c(14, 16:26)

datasets_csp <- read_CSP_data(waves, state, state_code, vaccine_1,
                              weight, EndDate, StartDate, wave, gender,
                              race, age, education, urbanicity, region) 


### read and process CDC and pop data ####

cdc <- read_preprocess_CDC_vax_data()

state_pop <- read_preprocess_pop_data() %>%
  filter(!State == "Puerto Rico") %>%
  mutate(Population_perc = Population/first(Population))


#### Calculate vax rate ####

csp_vax <- calculate_csp_vax(datasets_csp, "weight") 

datasets_csp <- lapply(datasets_csp, function(dat) mutate(dat, weight_unwt = 1))
csp_vax_unwt <- calculate_csp_vax(datasets_csp, "weight_unwt") 


#### Downsample based on state pop data ####

population_perc <- state_pop %>% select(state_code, Population_perc)


state_downsample_iter <- function(i) {
  out <- downsample_csp_state(datasets_csp, population_perc)
  
  print("Fraction kept:")
  print(out[[2]])
  
  datasets_state <- out[[1]]
  
  csp_vax_state <- calculate_csp_vax(datasets_state) %>%
    mutate(iteration = i)
  
  csp_vax_state_unwt <- lapply(
    datasets_state, function(d) 
      d %>%
      mutate(EndDate = as.Date(EndDate),
             StartDate = as.Date(StartDate),
             days = EndDate - min(StartDate),
             Date = min(StartDate) + median(days)) %>%
      summarise(vax = 100*mean(vaccine_1, na.rm = T), 
                Date = first(Date), 
                Wave = first(wave))
  ) %>%
    bind_rows()
  
  return(list(csp_vax_state, csp_vax_state_unwt))
}


calc_vax_sims <- function(sims) {
  
  sims %>% 
    bind_rows %>%
    group_by(Wave) %>%
    summarise(vax_mean = mean(vax),
              moe_mean = mean(moe),
              vax_sd = sd(vax),
              vax_range = max(abs(range(vax)-vax))) 
}

calc_diff <- function(sims_sum) {
  sims_sum %>%
    left_join(csp_vax, by = c("Wave")) %>%
    mutate(diff = vax_mean - vax,
           moe = moe + moe_mean) %>%
    select(diff, moe, Date)
}



n_iter <- 100
#sims_state <- lapply(c(1:n_iter), state_downsample_iter) 
#saveRDS(sims_state, "PROJECTS/Validation/paper/Simulations/sims_state.RDS")
sims_state<- readRDS("Simulations/sims_state.RDS")



#### Diff plot ####


vax_state_wt_diff <- lapply(sims_state, function(x) x[[1]]) %>%
  calc_vax_sims() %>% calc_diff()


plot_diff <- ggplot(vax_state_wt_diff, 
                    aes(x = Date, y = diff)) + 
  geom_pointrange(aes(ymin = (diff - moe), 
                      ymax = (diff + moe)),
                  position = position_dodge2(width = 10)) +
  geom_hline(yintercept = 0) +
  scale_x_date(breaks = "2 months", 
               date_labels = "%b%y") +
  xlab("Median date of survey wave") + 
  ylab("Difference between simulation whole sample \nestimates of vaccination rate")+
  theme_ipsum()+
  theme(legend.position.inside=c(0.12, 0.2), 
        legend.justification='left',
        legend.direction='vertical',
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(color='black',face='bold', size=14,
                                    margin = margin(0,8,0,0)),
        axis.title.x = element_text(color='black',face='bold', size=14,
                                    margin = margin(5,0,0,0)),
        legend.text = element_text(color='black', size=15),
        legend.title = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt"))


ggsave("Plots/state_downsample.png",
       plot_diff, 
       height = 5, width = 7,
       bg = "white", dpi = 300)



