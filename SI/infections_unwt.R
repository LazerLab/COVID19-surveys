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
library(labelled)


#### READ CSP DATA ####

waves <- c(5,7,9:14, 16:26)

### reads CSP data as a list of datasets, one per wave. The variables to keep from the full datasets are 
### provided as arguments

### Set all individual weights to 1 for unweighted analysis

datasets_csp <- read_CSP_data(waves, wave, StartDate, EndDate, weight, weight_state, state, state_code, id, covid, 
                              starts_with("cov_test"),
                              starts_with("cov_month"), gender, age_cat_4, age_cat_6, income, race, 
                              education_cat, urban_type, party3, employ) %>%
  lapply(function(d) mutate(d, wave_date = max(EndDate),
                            EndDate_numeric = as.numeric(as.Date(EndDate)),
                            weight = 1,
                            weight_state = 1)) 


#### RUN code as in infections.R ####

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


### apply 'calculate_infected_method1_states' to each dataset and calculate the cdf and
### the cumulative error
method1_results_states <- lapply(datasets_csp_method1, calculate_infected_method1_states) %>%
  bind_rows() %>%
  group_by(state_code) %>%
  mutate(cdf_estimate = cumsum(estimate),
         cdf_std = sqrt(cumsum(std_error^2))) %>%
  ungroup() %>%
  left_join(months_per_wave) %>%
  mutate(date_plot  = (month_max %m+% months(1)) - 1)



#### CALCULATE METHOD 2 ####

datasets_csp <- generate_var_method2(datasets_csp)

method2_results <- lapply(datasets_csp, calculate_infected_method2) %>%
  bind_rows() %>%
  mutate(weighted_median_date = as.Date(weighted_median_date))




#### NATIONAL TRENDS PLOT ####


plot_cdf_all <- ggplot()+
  geom_line(data = method1_results,
            aes(x = date_plot,
                y = cdf_estimate,
                color = 'CSP Method 1'), 
            lwd=0.95) +
  geom_line(data = method2_results,
            aes(x = weighted_median_date,
                y = estimate,
                color='CSP Method 2'),
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
  geom_point(data = method2_results,
             aes(x = weighted_median_date,
                 y = estimate,
                 color= 'CSP Method 2'),
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
  geom_errorbar(data = method2_results,
                aes(x=weighted_median_date,
                    ymin=(estimate)-(std_error*1.96),
                    ymax=(estimate)+(std_error*1.96),
                    color='CSP Method 2'),
                width=4)+
  geom_errorbar(data = n_infections_ipsos,
                aes(x=Date,
                    ymin=(avg_infections)-(std_error*1.96),
                    ymax=(avg_infections)+(std_error*1.96),
                    color='Ipsos'),
                width=4)+
  geom_vline(xintercept = as.numeric(as.Date("2022-02-01")),
             linetype='dashed', color = "chartreuse4", lwd = 0.6) +
  annotate("text", x = as.Date("2022-07-05"), y = 10, 
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
  scale_color_manual(values=c('steelblue1', 'darkorchid1', 'darkred', 'black'))+
  scale_x_date(date_breaks = "2 months",
               date_labels = "%b%y",
               limits = as.Date(c("2020-01-01","2023-01-31"))) 


#### NATIONAL DIFFERENCES PLOT ####


national_diff <- calculate_infections_diff(us_cumulative_nyt, method1_results,
                                           method2_results, n_infections_ipsos)


plot_diff <- ggplot()+
  
  geom_point(data = national_diff,
             aes(x = date,
                 y = method1_diff,
                 color= 'CSP Method 1'),
             size=2) +
  geom_point(data = national_diff,
             aes(x = date,
                 y = method2_diff,
                 color= 'CSP Method 2'),
             size=2) +
  geom_point(data = national_diff,
             aes(x = date,
                 y = ipsos_diff,
                 color = 'Ipsos'),
             size = 2) +
  
  geom_errorbar(data=national_diff,
                aes(x= date, 
                    ymin= method1_diff - method1_std*1.96,
                    ymax= method1_diff + method1_std*1.96,
                    color= 'CSP Method 1'),
                width=4) +
  geom_errorbar(data = national_diff,
                aes(x=date,
                    ymin=(method2_diff)-(method2_std*1.96),
                    ymax=(method2_diff)+(method2_std*1.96),
                    color= 'CSP Method 2'),
                width=4)+
  geom_errorbar(data = national_diff,
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
       y = "Difference between survey and\nclosest nofficial estimate")+
  scale_color_manual(values=c('steelblue1', 'darkorchid1', 'darkred'))+
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
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")) #Note that 5.5pt is the default plot margin of a ggplot2 plot, i.e. in the code below we are only changing the area margin at the bottom of our graph.



## states diffs
method1_results_states_join <- method1_results_states %>% 
  left_join(states_cumulative_nyt, 
            by = c("date_plot" = "date",
                   "state" = "state")) %>%
  mutate(cdf_estimate_diff = cdf_estimate - pct_cumulative_cases,
         `Difference within 95% CI` = (pct_cumulative_cases < (cdf_estimate + (cdf_std*1.96)) ) & 
           (pct_cumulative_cases > (cdf_estimate - (cdf_std*1.96)) ) ) %>%
  mutate(Period = case_when(date_plot > as.Date("2021-11-30") ~ "Post-rapid",
                            date_plot <= as.Date("2021-11-30") ~ "Pre-rapid"))


state_histogram_infections <- ggplot(method1_results_states_join,
                                     aes(x = cdf_estimate_diff,
                                         y = 100*after_stat(density),
                                         fill = Period,
                                         color = Period)) +
  geom_histogram(position = "identity",
                 alpha = 0.4, 
                 key_glyph = "rect") +
  geom_vline(data = method1_results_states_join %>% group_by(Period) %>% 
               summarise(average = median(cdf_estimate_diff)),
             aes(xintercept = average, color = Period),
             linetype = "dashed", 
             linewidth = 1) + 
  scale_color_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2") +
  xlab("Difference between unweighted CSP and official state estimate of cumulative infection rate") +
  ylab("Percentage of state-wave pairs") +
  theme_ipsum() +
  scale_y_continuous(#minor_breaks = c(-5,5),
    limits = c(0, 30)) +
  theme(
    legend.position=c(0.01, 0.701), 
    legend.justification='left',
    legend.direction='vertical',
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(color='black', size=12),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    axis.title.y = element_text(color='black',face='bold', size=14,
                                margin = margin(0,10,0,0)),
    axis.title.x = element_text(color='black',face='bold', size=14),
    panel.border = element_blank(),
    plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))

state_histogram_infections


# MERGE PLOTS
full_plot <- plot_cdf_all / plot_diff / state_histogram_infections + 
  plot_annotation(tag_levels = list(c('a)', 'b)', 'c)')),
  ) &
  theme(plot.title = element_text(size = 16),
        plot.tag = element_text(face = 'bold',
                                #size = 8,
                                hjust = -5,
                                vjust = 5))

ggsave("Plots/infections_full_plot_unwt.png",
       full_plot, 
       height = 14, width = 10,
       bg = "white", dpi = 300)

