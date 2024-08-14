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
datasets_csp <- read_CSP_data(waves, wave, StartDate, EndDate, weight, weight_state, state, state_code, id, covid, 
                              starts_with("cov_test"), starts_with("cov_trust"), source,
                              starts_with("cov_month"), gender, age_cat_4, age_cat_6, income, race, 
                              education_cat, urban_type, party3, employ) %>%
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

method2_results_states <- lapply(datasets_csp, calculate_infected_method2_states) %>%
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
plot_cdf_all

## Text numbers: Diff between method1 and method 2

library(slider)

slide(method1_results, function(row){
  est1 <- method2_results %>% mutate(d_diff = abs(weighted_median_date - row$date_plot)) %>%
    filter(d_diff == min(d_diff)) %>% pull(estimate)
  
  return(abs(est1 - row$cdf_estimate))
})


#### NATIONAL DIFFERENCES PLOT ####


national_diff <- calculate_infections_diff(us_cumulative_nyt, method1_results,
                                           method2_results, n_infections_ipsos)


plot_diff<- ggplot()+
  
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

plot_diff



#### STATE ANALYSYS ####


## states diffs
method1_results_states_join <- method1_results_states %>% 
  left_join(states_cumulative_nyt, 
            by = c("date_plot" = "date",
                   "state" = "state")) %>%
  mutate(cdf_estimate_diff = cdf_estimate - pct_cumulative_cases,
         `Difference within 95% CI` = (pct_cumulative_cases < (cdf_estimate + (cdf_std*1.96)) ) & 
                                      (pct_cumulative_cases > (cdf_estimate - (cdf_std*1.96)) ) ) 

method2_results_states_join <- method2_results_states %>% 
  left_join(states_cumulative_nyt, 
            by = c("weighted_median_date" = "date",
                   "state" = "state")) %>%
  mutate(cdf_estimate_diff = estimate - pct_cumulative_cases,
         `Difference within 95% CI` = (pct_cumulative_cases < (estimate + (std_error*1.96)) ) & 
           (pct_cumulative_cases > (estimate - (std_error*1.96)) ) ) 


#pre and post rapid separation
method1_results_states_join_prerapid <- method1_results_states_join %>%
  filter(date_plot <= as.Date("2021-11-30")) ### wave 22 has monthmin at 2021-12-01 and monthmax 2022-02-01. we're leaving it to postrapid

method1_results_states_join_postrapid <- method1_results_states_join %>%
  filter(date_plot >= as.Date("2021-12-01"))


method2_results_states_join_prerapid <- method2_results_states_join %>%
  filter(weighted_median_date <= as.Date("2021-11-30")) ### wave 22 has monthmin at 2021-12-01 and monthmax 2022-02-01. we're leaving it to postrapid

method2_results_states_join_postrapid <- method2_results_states_join %>%
  filter(weighted_median_date >= as.Date("2021-12-01"))


### Text numbers: SUMMARY Tables

sum_tab_prerapid_M1 <- method1_results_states_join_prerapid %>%
  group_by(state) %>%
  summarise(`Number of waves` = n(),
            `Differences within margin of error` = sum(1.96*cdf_std >= abs(cdf_estimate_diff)),
            `Estimates within 10% of official estimate` = sum(abs(cdf_estimate_diff) <= 10),
            `Estimates within 5% of official estimate` = sum(abs(cdf_estimate_diff) <= 5),
            `Estimates within 2.5% of official estimate` = sum(abs(cdf_estimate_diff) <= 2.5),
            `Estimates within 1% of official estimate` = sum(abs(cdf_estimate_diff) <= 1))

sum_tab_prerapid_M2 <- method2_results_states_join_prerapid %>%
  group_by(state) %>%
  summarise(`Number of waves` = n(),
            `Differences within margin of error` = sum(1.96*std_error >= abs(cdf_estimate_diff)),
            `Estimates within 10% of official estimate` = sum(abs(cdf_estimate_diff) <= 10),
            `Estimates within 5% of official estimate` = sum(abs(cdf_estimate_diff) <= 5),
            `Estimates within 2.5% of official estimate` = sum(abs(cdf_estimate_diff) <= 2.5),
            `Estimates within 1% of official estimate` = sum(abs(cdf_estimate_diff) <= 1))


#reporting mean and median vals for both pre and postrapid histograms
median(method1_results_states_join_prerapid$cdf_estimate_diff)
median(method1_results_states_join_postrapid$cdf_estimate_diff)

mean(method1_results_states_join_prerapid$cdf_estimate_diff)
mean(method1_results_states_join_postrapid$cdf_estimate_diff)

### Text numbers: method 1 vs method 2

sum_tab_prerapid_M1 %>% summarise(across(3:7, ~sum(.)/
                                           (12*nrow(sum_tab_prerapid_M1))))

sum_tab_prerapid_M2 %>% summarise(across(3:7, ~sum(.)/
                                           (13*nrow(sum_tab_prerapid_M2))))

method1_results_states_join_prerapid %>% summarise(mean(abs(cdf_estimate_diff)))
method2_results_states_join_prerapid %>% summarise(mean(abs(cdf_estimate_diff)))

method1_results_states_join_prerapid %>% summarise(median(abs(cdf_estimate_diff)))
method2_results_states_join_prerapid %>% summarise(median(abs(cdf_estimate_diff)))

method1_results_states_join_prerapid %>% summarise(mean(cdf_std))
method2_results_states_join_prerapid %>% summarise(mean(std_error))



method1_results_states_join <- method1_results_states_join %>%
  mutate(Period = case_when(date_plot > as.Date("2021-11-30") ~ "Post-rapid",
                            date_plot <= as.Date("2021-11-30") ~ "Pre-rapid"))


# STATE HISTOGRAM

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
  xlab("Difference between CSP and official state estimate of cumulative infection rate") +
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

ggsave("Plots/infections_full_plot.png",
       full_plot, 
       height = 14, width = 10,
       bg = "white", dpi = 300)



### STATE TRENDS PLOT

plot_cdf_states <- ggplot()+
  
  geom_line(data = method1_results_states,
            aes(x = date_plot,
                y = cdf_estimate,
                color = 'CSP Method 1'), 
            lwd=0.75) +
  geom_line(data = method2_results_states,
            aes(x = weighted_median_date,
                y = estimate,
                color='CSP Method 2'),
            lwd=0.75) +
  
  geom_point(data = method1_results_states,
             aes(x = date_plot,
                 y = cdf_estimate,
                 color = 'CSP Method 1'),
             size=1) +
  geom_point(data = method2_results_states,
             aes(x = weighted_median_date,
                 y = estimate,
                 color = 'CSP Method 2'),
             size=0.9) +
  
  geom_errorbar(data=method1_results_states,
                aes(x= date_plot, 
                    ymin= cdf_estimate - cdf_std*1.96,
                    ymax= cdf_estimate + cdf_std*1.96,
                    color = 'CSP Method 1'),
                width=2) +
  geom_errorbar(data = method2_results_states,
                aes(x=weighted_median_date,
                    ymin=(estimate)-(std_error*1.96),
                    ymax=(estimate)+(std_error*1.96),
                    color = 'CSP Method 2'),
                width=2)+
  
  geom_line(data = states_cumulative_nyt,
            aes(x = date,
                y = pct_cumulative_cases,
                color='Official cumulative case rate'),
            lwd=0.75)+
  facet_wrap(~state, nrow = 9) +
  
  geom_vline(xintercept = as.numeric(as.Date("2022-02-01")),
             linetype='dashed', color = "chartreuse4", lwd = 0.6) +
  labs(x = "Last day of pandemic period (red) and survey wave (blue)",
       y = "COVID-19 Infections (%)")+
  scale_color_manual(values=c('steelblue1', 'darkorchid1', 'black'), name = NULL)+
  scale_x_date(date_breaks = "4 months",
               date_labels = "%b%y",
               limits = as.Date(c("2020-03-01","2023-01-31"))) +
  theme_ipsum()+
  theme(legend.position='top', 
        legend.justification='center',
        legend.direction='horizontal',
        legend.box.background = element_rect(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(color='black',face='bold', size=13),
        axis.title.x = element_text(color='black',face='bold', size=13),
        axis.title.y.right = element_text(color='blue',face='bold', size=15),
        #axis.title.x = element_blank(),
        legend.text = element_text(color='black', size=12),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.1, "lines"))

ggsave("Plots/infections_state.png",
       plot_cdf_states, 
       height = 14, width = 12,
       bg = "white", dpi = 300)


## STATE DIFFS PLOT

plot_cdf_states_diff <- ggplot() + 
  
  geom_line(data = method1_results_states_join_prerapid,
            aes(x = date_plot,
                y = cdf_estimate_diff,
                color = "CSP Method 1"),
            lwd=0.75) +
  
  geom_point(data = method1_results_states_join_prerapid,
             aes(x = date_plot,
                 y = cdf_estimate_diff,
                 color = "CSP Method 1"),
             size=1) +
  
  geom_errorbar(data=method1_results_states_join_prerapid,
                aes(x= date_plot, 
                    ymin= cdf_estimate_diff - cdf_std*1.96,
                    ymax= cdf_estimate_diff + cdf_std*1.96,
                    linetype = `Difference within 95% CI`,
                    color= "CSP Method 1"),
                width=2) +
  
  geom_line(data = method2_results_states_join_prerapid,
            aes(x = weighted_median_date,
                y = cdf_estimate_diff,
                color = "CSP Method 2"),
            lwd=0.75) +
  
  geom_point(data = method2_results_states_join_prerapid,
             aes(x = weighted_median_date,
                 y = cdf_estimate_diff,
                 color = "CSP Method 2"),
             size=1) +
  
  geom_errorbar(data=method2_results_states_join_prerapid,
                aes(x= weighted_median_date, 
                    ymin= cdf_estimate_diff - std_error*1.96,
                    ymax= cdf_estimate_diff + std_error*1.96,
                    linetype = `Difference within 95% CI`,
                    color = "CSP Method 2"),
                width=2) +
  
  facet_wrap(~state, ncol = 6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray")+
  
  labs(x = "Date",
       y = "Difference between CSP and official estimate of cumulative COVID-19 Infections (%)")+
  scale_color_manual(values=c('steelblue1', 'darkorchid1'), name = NULL) + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  scale_x_date(date_breaks = "3 months",
               date_labels = "%b%y") +
  scale_y_continuous(minor_breaks = c(-5,5),
                     limits = c(-8, 10),
                     oob = scales::squish) +
  
  theme_ipsum()+
  theme(legend.position='top', 
        legend.justification='center',
        legend.direction='horizontal',
        legend.box.background = element_rect(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(color='black',face='bold', size=13),
        axis.title.x = element_text(color='black',face='bold', size=13),
        axis.title.y.right = element_text(color='blue',face='bold', size=15),
        #axis.title.x = element_blank(),
        legend.text = element_text(color='black', size=12),
        panel.border = element_blank(),
        panel.spacing = unit(0.1, "lines")) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))

plot_cdf_states_diff

ggsave("Plots/infections_state_diff.png", plot_cdf_states_diff, 
       height = 14, width = 12,
       bg = "white", dpi = 300)




#### Extended data Figure 3: AT LEAST 1 INFECTION PLOT ####

### from ipsos: calculate percentage with at least 1 infection

datasets_ipsos_1_inf <-  Filter(function(d) "Q13" %in% colnames(d) 
                                & ("Q21" %in% colnames(d)) 
                                & ("WT_FINAL" %in% colnames(d)), datasets_ipsos)

datasets_ipsos_1_inf <- lapply(datasets_ipsos_1_inf, function(d) d %>%
                                 mutate(inf = case_when(Q13 == 1 & Q21 == 1 ~ 1,
                                                        Q13 == 2 | (Q13==1 & Q21==2) ~ 0,
                                                        T ~ NA_real_)) %>%
                                 select(WT_FINAL, WAVE, inf)) %>%
  bind_rows() %>%
  left_join(dates_ipsos, by = "WAVE")

infections_1_ipsos = datasets_ipsos_1_inf %>% group_by(WAVE) %>%
  group_map(~ as_survey_design(., weight = WT_FINAL) %>%
              summarise(pct_infected = survey_mean(inf, na.rm = T),
                        Date = as.Date(first(Date)),
                        Wave = first(WAVE),
                        wt = list(WT_FINAL)),
            .keep = T) %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(std_error = calc_err(pct_infected, wt)*100,
         pct_infected = 100*pct_infected) %>%
  select(-wt, -pct_infected_se) %>%
  ungroup()


### CSP: calculate percentage with at least 1 infection


infections_1_csp <- lapply(datasets_csp, function(d) d %>% as_survey_design(weights = weight) %>%
                    summarise(pct_infected = 100*survey_mean(cov_test == 1, na.rm = T),
                              weighted_median_date = survey_median(EndDate_numeric, na.rm = T),
                              wave = first(wave))) %>%
  bind_rows() %>%
  mutate(weighted_median_date = as.Date(weighted_median_date))


pct_tested_poz <- ggplot()+
  
  geom_line(data = infections_1_csp,
            aes(x = weighted_median_date,
                y = pct_infected,
                color = 'CSP'), 
            lwd=0.95) +
  geom_line(data = infections_1_ipsos,
            aes(x = Date,
                y = pct_infected,
                color='IPSOS'),
            lwd=0.95) +

  geom_point(data = infections_1_csp,
             aes(x = weighted_median_date,
                 y = pct_infected,
                 color= 'CSP'),
             size=1) +
  geom_point(data = infections_1_ipsos,
             aes(x = Date,
                 y = pct_infected,
                 color= 'IPSOS'),
             size=1) +

  geom_errorbar(data=infections_1_csp,
                aes(x= weighted_median_date, 
                    ymin= pct_infected - pct_infected_se*1.96,
                    ymax= pct_infected + pct_infected_se*1.96,
                    color='CSP'),
                width=3) +
  geom_errorbar(data = infections_1_ipsos,
                aes(x=Date,
                    ymin=(pct_infected)-(std_error*1.96),
                    ymax=(pct_infected)+(std_error*1.96),
                    color='IPSOS'),
                width=3)+
  
  theme_ipsum()+
  theme(#legend.position = 'right',
    #legend.position='top', 
    legend.position=c(0.05, 0.85), 
    legend.justification='left',
    legend.direction='vertical',
    #axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
    axis.title.y = element_text(color='black',face='bold', size=13),
    axis.title.x = element_text(color='black',face='bold', size=13),
    axis.title.y.right = element_text(color='blue',face='bold', size=15),
    #axis.title.x = element_blank(),
    legend.text = element_text(color='black', size=12),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt"))+ #Note that 5.5pt is the default plot margin of a ggplot2 plot, i.e. in the code below we are only changing the area margin at the bottom of our graph.)+
  labs(x = "Date",
       y = "Percentage who tested positive at least once")+
  scale_color_manual(values=c( 'blue', 'darkred'))+
  scale_x_date(date_breaks = "2 months",
               date_labels = "%b%y",
               limits = as.Date(c("2020-01-01","2023-01-31"))) 

pct_tested_poz

ggsave("Plots/at_least_1_infection.png", 
       height = 6, width = 9,
       bg = "white", dpi = 300)

## Text numbers ##

infections_1_comp <- infections_1_csp %>% 
  filter(wave <=23) %>%
  rename(csp_infected = pct_infected,
         csp_se = pct_infected_se) %>% 
  left_join(infections_1_ipsos %>%
              rename(ipsos_infected = pct_infected,
                     ipsos_se = std_error), join_by(closest(weighted_median_date >= Date))) %>%
  left_join(infections_1_ipsos %>%
              rename(ipsos_infected = pct_infected,
                     ipsos_se = std_error), join_by(closest(weighted_median_date <= Date))) %>%
  mutate(days_diff1 = abs(weighted_median_date - Date.x),
         days_diff2 = abs(weighted_median_date - Date.y),
         days_diff = pmin(days_diff1, days_diff2, na.rm = T),
         infected_diff = case_when(days_diff2 >= days_diff1 | is.na(days_diff2) ~
                                     csp_infected - ipsos_infected.x,
                              (days_diff1 > days_diff2)  ~ 
                                csp_infected - ipsos_infected.y), 
         moe_total = case_when(days_diff2 >= days_diff1 | is.na(days_diff2) ~ 
                                 csp_se*1.96 + ipsos_se.x*1.96,
                               (days_diff1 > days_diff2) ~
                                 csp_se*1.96 + ipsos_se.y*1.96)) %>%
  select(weighted_median_date, days_diff, infected_diff, moe_total)

infections_1_comp %>% 
  summarise(sum(abs(infected_diff)>moe_total),
            mean(infected_diff))
nrow(infections_1_comp)



##### Waves Info SI Table 1 ####

tab <- bind_rows(datasets_csp) %>% group_by(wave) %>% 
  summarise(wave = first(wave), 
            date1 = format(min(StartDate), "%m/%d/%y"),
            date2 = format(max(StartDate), "%m/%d/%y"), 
            n = n()) %>%
  mutate(Field = str_c(date1," - ", date2))



##### Infections by trust: Extended data Figure 9 ####

calculate_infected_method2_trust <- function(dat, var){
  dat %>% as_survey_design(ids = 1, weights = weight) %>%
    group_by(!!sym(var)) %>%
    summarize(wave = first(wave),
              weighted_median_date = survey_median(EndDate_numeric, na.rm = T),
              wave_end_date = as.Date(first(wave_date)),
              estimate = 100*survey_mean(count_cases, na.rm = T)) %>%
    rename(std_error = estimate_se) %>%
    filter(!is.na(!!sym(var)))
} 

method2_results_trust_cdc <- lapply(datasets_csp[1:17], calculate_infected_method2_trust,
                                    "cov_trust_cdc") %>%
  bind_rows() %>%
  mutate(weighted_median_date = as.Date(weighted_median_date)) %>% 
  mutate(cov_trust_cdc = to_factor(cov_trust_cdc)) 

trust_cdc <- ggplot(method2_results_trust_cdc, 
       aes(x = weighted_median_date,
           y = estimate,
           group = cov_trust_cdc)) +
  geom_pointrange(aes(ymin = (estimate - std_error*1.96), 
                      ymax = (estimate + std_error*1.96),
                      color = cov_trust_cdc),
                  position = position_dodge2(width = 10)) +
  ylim(c(-2, 43)) +
  xlab("Median date of survey wave") + 
  ylab("Cumulative infection rate") +
  scale_color_viridis_d(name = "Trust in CDC to handle \n the COVID-19 pandemic") +
  theme_bw() +
  theme(legend.position = "bottom")


method2_results_trust_science <- lapply(datasets_csp[1:17], calculate_infected_method2_trust, "cov_trust_science") %>%
  bind_rows() %>%
  mutate(weighted_median_date = as.Date(weighted_median_date)) %>% 
  mutate(cov_trust_science = to_factor(cov_trust_science)) 


trust_science <- ggplot(method2_results_trust_science, 
       aes(x = weighted_median_date,
           y = estimate,
           group = cov_trust_science)) +
  geom_pointrange(aes(ymin = (estimate - std_error*1.96), 
                      ymax = (estimate + std_error*1.96),
                      color = cov_trust_science),
                  position = position_dodge2(width = 10)) +
  ylim(c(-3, 45)) +
  xlab("Median date of survey wave") + 
  ylab("Cumulative infection rate") +
  scale_color_viridis_d(name = "Trust in scientists to handle \n the COVID-19 pandemic") +
  theme_bw() +
  theme(legend.position = "bottom")

trust_bw_plot <- trust_cdc + trust_science

ggsave("Plots/trust_bw_infections.png",
       trust_bw_plot,
       height = 7, width = 12,
       bg = "white", dpi = 300)







##### SI Table 5 (full state samples) ####

state_samples <- lapply(datasets_csp, function(d) d %>% group_by(state_code) %>% 
         summarise(n = n())) %>%
  reduce(left_join, by = "state_code")


