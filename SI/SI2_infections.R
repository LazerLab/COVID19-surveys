setwd("COVID19/CODE/PROJECTS/Validation/paper") ### set to repository location

source("functions_read_data.R")
source("functions_analysis.R")

library(tidyverse)
library(survey)
library(lubridate)
library(srvyr)
library(hrbrthemes)
library(ggpubr)
library(patchwork)


#### READ CSP DATA ####

waves <- c(5,7,9:14, 16:26)

### reads CSP data as a list of datasets, one per wave. The variables to keep from the full datasets are 
### provided as arguments
datasets_csp <- read_CSP_data(waves, wave, StartDate, EndDate, weight, weight_state, state, state_code, id, covid, 
                              starts_with("cov_test"), starts_with("cov_trust"), age, urbanicity, region,
                              starts_with("cov_month"), gender, age_cat_4, age_cat_6, income, race, 
                              education_cat, urban_type, party3, employ, education) %>%
  lapply(function(d) mutate(d, wave_date = max(EndDate),
                            EndDate_numeric = as.numeric(as.Date(EndDate)))) 


#### READ ADMIN DATA ####

us_cumulative_nyt <- read_infections_national_admin_data()

states_cumulative_nyt <- read_infections_state_admin_data()


#### IPSOS Data ####


### read and process

out <- read_preprocess_IPSOS_data(date_return = "last")

datasets_ipsos <- out[[1]]
dates_ipsos <- out[[2]]


waves1 <- sapply(seq_along(datasets_ipsos), function(i) {
  if ("Q174" %in% colnames(datasets_ipsos[[i]])) {
    return(i)
  }
}) %>% unlist()

datasets_ipsos_n_inf <- datasets_ipsos[waves1]

datasets_ipsos_n_inf <- lapply(datasets_ipsos_n_inf, function(d) d %>% 
                                 mutate(n_infections = case_when(!is.na(Q174) & Q174 != -1 & Q21 == 1 ~ as.numeric(Q174),
                                                                 Q13 == -1 | Q21 == -1 | Q174 == -1 ~ NA_real_, 
                                                                 T ~  0))
)


datasets_ipsos_n_inf <- bind_rows(
  lapply(datasets_ipsos_n_inf, function(df)
    df %>% select(WT_FINAL, WAVE, n_infections))
) %>%
  left_join(dates_ipsos, by = "WAVE")

### calculate average infections per respondent

n_infections_ipsos <- calculate_ipsos_infections(datasets_ipsos_n_inf)



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



####  Robustness check ####

#Sampling step (picking only 1 response from each returning respondent)
nrow(bind_rows(datasets_csp_method1))
#408,515 responses over 17 waves (method 1)
sampled_dataset_method1 <- bind_rows(datasets_csp_method1) %>% group_by(id) %>% sample_n(1)
nrow(sampled_dataset_method1)
#306,799 after de-duplication

#Assert: check if there exists duplicates. the following should yield 0:
sum(duplicated(sampled_dataset_method1 %>% select(id))) 

#Calculate new weights after discarding returning respondents:
split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col])

sampled_datasets_method1 <- split_tibble(sampled_dataset_method1, 'wave') %>%
  lapply(function(d) d %>% add_column( robustness_national_weights =  weight_csp(dat=d,
                                                                                 #wave=2,
                                                                                 state=FALSE))
  )


#Run the same steps as we did above, but with weights = robustness_national_weights this time:
method1_robustness_results <- lapply(sampled_datasets_method1, calculate_infected_method1, 
                                     weight_var = "robustness_national_weights") %>%
  bind_rows() %>%
  mutate(cdf_estimate = cumsum(estimate),
         cdf_std = sqrt(cumsum(std_error^2))) %>%
  left_join(months_per_wave) %>%
  mutate(date_plot  = (month_max %m+% months(1)) - 1)
  
#Plots:
plot_cdf_robustness_method1 <- ggplot()+
  geom_line(data = method1_results,
            aes(x = date_plot,
                y = cdf_estimate,
                color = 'CSP Method 1'), 
            lwd=0.95) +
  geom_line(data = method1_robustness_results,
            aes(x = date_plot,
                y = cdf_estimate,
                color='CSP Method 1 - robustness'),
            lwd=0.95) +
  geom_line(data = n_infections_ipsos,
            aes(x = Date,
                y = avg_infections,
                color='Ipsos'),
            lwd=0.95) +
  geom_line(data = us_cumulative_nyt,
            aes(x = date,
                y = pct_cumulative_cases,
                color='Official cumulative infection rate'),
            alpha = .85,
            lwd=0.95)+
  geom_point(data = method1_results,
             aes(x = date_plot,
                 y = cdf_estimate,
                 color= 'CSP Method 1'),
             size=2) +
  geom_point(data = method1_robustness_results,
             aes(x = date_plot,
                 y = cdf_estimate,
                 color= 'CSP Method 1 - robustness'),
             size=2) +
  geom_point(data = n_infections_ipsos,
             aes(x = Date,
                 y = avg_infections,
                 color = 'Ipsos'),
             size = 2) +
  geom_errorbar(data=method1_results,
                aes(x= date_plot, 
                    ymin= cdf_estimate - cdf_std*1.96,
                    ymax= cdf_estimate + cdf_std*1.96,
                    color='CSP Method 1'),
                width=4) +
  geom_errorbar(data = method1_robustness_results,
                aes(x=date_plot,
                    ymin=(cdf_estimate)-(cdf_std*1.96),
                    ymax=(cdf_estimate)+(cdf_std*1.96),
                    color='CSP Method 1 - robustness'),
                width=4)+
  geom_errorbar(data = n_infections_ipsos,
                aes(x=Date,
                    ymin=(avg_infections)-(std_error*1.96),
                    ymax=(avg_infections)+(std_error*1.96),
                    color='Ipsos'),
                width=4)+
  geom_vline(xintercept = as.numeric(as.Date("2022-02-01")),
             linetype='dashed', color = "chartreuse4", lwd = 0.6) +
  annotate("text", x = as.Date("2022-07-15"), y = 10, 
           label = "Mass distribution of rapid tests",color = "black") + 
  theme_ipsum()+
  theme(legend.position=c(0.05, 0.85), 
        legend.justification='left',
        legend.direction='vertical',
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(color='black',face='bold', size=14),
        axis.title.x = element_text(color='black',face='bold', size=14),
        legend.text = element_text(color='black', size=15),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt"))+ #Note that 5.5pt is the default plot margin of a ggplot2 plot, i.e. in the code below we are only changing the area margin at the bottom of our graph.
  labs(x = "Date",
       y = "Cumulative rate of COVID-19 Infections (%)")+
  scale_color_manual(values=c('steelblue1', 'steelblue4', 'darkred', 'black'))+
  scale_x_date(date_breaks = "2 months",
               date_labels = "%b%y",
               limits = as.Date(c("2020-01-01","2023-01-31"))) 
plot_cdf_robustness_method1




#Calculate and plot the differences to see the effect of removing returners
#(CSP method 1, CSP method 1 but returners removed, and Ipsos; minus the official numbers:

robustness_diff_method1 <- us_cumulative_nyt %>%
  rename(nyt_pct_cumulative_cases = pct_cumulative_cases) %>%
  full_join(select(method1_results, date_plot, cdf_estimate, cdf_std) %>%
              rename(method1_pct_cumulative_cases = cdf_estimate,
                     method1_std = cdf_std,
                     date = date_plot), 
            by = "date")  %>%
  full_join(method1_robustness_results %>%
              rename(method1_robustness_pct_cumulative_cases = cdf_estimate,
                     method1_robustness_std = cdf_std,
                     date = date_plot), 
            by = "date") %>% 
  full_join(n_infections_ipsos %>% 
              rename(ipsos_pct_cumulative_cases = avg_infections,
                     ipsos_std = std_error, 
                     date = Date), 
            by = "date") %>%
  mutate(method1_diff =  method1_pct_cumulative_cases - nyt_pct_cumulative_cases,
         method1_robustness_diff = method1_robustness_pct_cumulative_cases - nyt_pct_cumulative_cases,
         ipsos_diff = ipsos_pct_cumulative_cases - nyt_pct_cumulative_cases) %>%
  filter(!is.na(method1_diff) | !is.na(method1_robustness_diff) | !is.na(ipsos_diff)) %>%
  select(date, method1_diff, method1_robustness_diff, ipsos_diff, method1_std, method1_robustness_std, ipsos_std)

plot_diff_robustness_method1 <- ggplot()+
  geom_point(data = robustness_diff_method1,
             aes(x = date,
                 y = method1_diff,
                 color= 'CSP Method 1'),
             size=2) +
  geom_point(data = robustness_diff_method1,
             aes(x = date,
                 y = method1_robustness_diff,
                 color= 'CSP Method 1 - robustness'),
             size=2) +
  geom_point(data = robustness_diff_method1,
             aes(x = date,
                 y = ipsos_diff,
                 color = 'Ipsos'),
             size = 2) +
  geom_errorbar(data=robustness_diff_method1,
                aes(x= date, 
                    ymin= method1_diff - method1_std*1.96,
                    ymax= method1_diff + method1_std*1.96,
                    color= 'CSP Method 1'),
                width=4) +
  geom_errorbar(data = robustness_diff_method1,
                aes(x=date,
                    ymin=(method1_robustness_diff)-(method1_robustness_std*1.96),
                    ymax=(method1_robustness_diff)+(method1_robustness_std*1.96),
                    color= 'CSP Method 1 - robustness'),
                width=4)+
  geom_errorbar(data = robustness_diff_method1,
                aes(x=date,
                    ymin=(ipsos_diff)-(ipsos_std*1.96),
                    ymax=(ipsos_diff)+(ipsos_std*1.96),
                    color = 'Ipsos'),
                width=4)+
  geom_vline(xintercept = as.numeric(as.Date("2022-02-01")),
             linetype='dashed', color = "chartreuse4", lwd = 0.6) +
  annotate("text", x = as.Date("2021-08-10"), y = 10, 
           label = "Mass distribution of rapid tests",color = "black") + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray', lwd = 0.6) +
  labs(x = "Date",
       y = "Difference between survey and\nclosest official estimate")+
  scale_color_manual(values=c('steelblue1', 'steelblue4', 'darkred'))+
  scale_x_date(date_breaks = "2 months",
               date_labels = "%b%y",
               limits = as.Date(c("2020-01-01","2023-01-31")))+
  scale_y_continuous(limits = c(-4,20)) +
  theme_ipsum()+
  theme(legend.position=c(0.05, 0.85), 
        legend.justification='left',
        legend.direction='vertical',
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(color='black',face='bold', size=14),
        axis.title.x = element_text(color='black',face='bold', size=14),
        legend.text = element_text(color='black', size=15),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")) #Note that 5.5pt is the
#default plot margin of a ggplot2 plot, i.e. here we are only
#changing the area margin at the bottom of our graph.

plot_diff_robustness_method1

#See both the estimate trendlines and the divergence from official numbers
#all together via patchwork:

plot_cdf_robustness_method1 / plot_diff_robustness_method1

ggsave("Plots/infections_robustness_method1_reweighted_panel.png", 
       plot_cdf_robustness_method1 / plot_diff_robustness_method1, 
       height = 12, width = 9,
       bg = "white", dpi = 300)




