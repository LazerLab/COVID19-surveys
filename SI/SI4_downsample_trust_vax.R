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


waves <- c(14, 16:22, 24)

datasets_csp <- read_CSP_data(waves)


#### Calculate vax rate ####

csp_vax <- calculate_csp_vax(datasets_csp, "weight") 


#### downsample based on trust ####


trust_downsample_iter <- function(i, fraction_keep) {
  out <- downsample_csp_trust(datasets_csp, fraction_keep)
  
  print("Fraction kept:")
  print(out[[2]])
  
  datasets_trust <- out[[1]]
  
  csp_vax_trust <- calculate_csp_vax(datasets_trust) %>%
    mutate(iteration = i)
  
  return(list(csp_vax_trust, out[[2]]))
}


calc_vax_sims <- function(sims) {
  
  sims %>% 
    bind_rows() %>%
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

fraction_keep1 <-
  tibble(cov_trust_cdc = c(1,2,3,4, -99),
       fraction_keep=c(0.7,0.8,0.9,1, 1))
#sims_trust1 <- lapply(c(1:n_iter), trust_downsample_iter, fraction_keep1)
sims_trust1 <- readRDS("Simulations/sims_trust1.RDS")

fraction_keep2 <-
  tibble(cov_trust_cdc = c(1,2,3,4, -99),
         fraction_keep=c(0.5,0.6,0.9,1, 1))
#sims_trust2 <- lapply(c(1:n_iter), trust_downsample_iter, fraction_keep2)
sims_trust2 <- readRDS("Simulations/sims_trust2.RDS")


fraction_keep3 <-
  tibble(cov_trust_cdc = c(1,2,3,4, -99),
         fraction_keep=c(0.3,0.4,0.9,1, 1))
#sims_trust3 <- lapply(c(1:n_iter), trust_downsample_iter, fraction_keep3)
sims_trust3 <- readRDS("Simulations/sims_trust3.RDS")

fraction_keep4 <-
  tibble(cov_trust_cdc = c(1,2,3,4, -99),
         fraction_keep=c(0.2,0.2,1,1, 1))
#sims_trust4 <- lapply(c(1:n_iter), trust_downsample_iter, fraction_keep4)
sims_trust4 <- readRDS("Simulations/sims_trust4.RDS")



# saveRDS(sims_trust1, "PROJECTS/Validation/paper/Simulations/sims_trust1.RDS")
# saveRDS(sims_trust2, "PROJECTS/Validation/paper/Simulations/sims_trust2.RDS")
# saveRDS(sims_trust3, "PROJECTS/Validation/paper/Simulations/sims_trust3.RDS")
# saveRDS(sims_trust4, "PROJECTS/Validation/paper/Simulations/sims_trust4.RDS")



vax_sims1_diff <- lapply(sims_trust1, function(x) x[[1]]) %>%
  calc_vax_sims() %>%
  calc_diff() %>%
  mutate(`Trust based downsample` = "Low")

vax_sims2_diff <- lapply(sims_trust2, function(x) x[[1]]) %>%
  calc_vax_sims() %>% 
  calc_diff() %>%
  mutate(`Trust based downsample` = "Mild")

vax_sims3_diff <- lapply(sims_trust3, function(x) x[[1]]) %>%
  calc_vax_sims() %>%
  calc_diff()  %>%
  mutate(`Trust based downsample` = "High")

vax_sims4_diff <- lapply(sims_trust4, function(x) x[[1]]) %>%
  calc_vax_sims() %>%
  calc_diff()  %>%
  mutate(`Trust based downsample` = "Very High")

diff_df <- bind_rows(vax_sims1_diff, vax_sims2_diff,
                     vax_sims3_diff, vax_sims4_diff) %>%
  mutate(`Trust based downsample` =
           factor(`Trust based downsample`,
                   levels = c("Low", "Mild", "High", "Very High")))

plot_diff <- ggplot(diff_df, 
                    aes(x = Date, y = diff, 
                        colour = `Trust based downsample`)) + 
  geom_pointrange(aes(ymin = (diff - moe), 
                      ymax = (diff + moe)),
                  position = position_dodge2(width = 10)) +
  scale_colour_viridis_d()+
  geom_hline(yintercept = 0) +
  scale_x_date(breaks = "2 months", 
               date_labels = "%b%y") +
  xlab("Median date of survey wave") + 
  ylab("Difference between simulation whole sample \nestimates of vaccination rate")+
  theme_ipsum()+
  theme(legend.position=c(0, 0.8), 
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

ggsave("Plots/trust_downsample_vax.png",
       plot_diff, 
       height = 6.5, width = 9,
       bg = "white", dpi = 300)

