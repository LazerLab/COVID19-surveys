setwd("COVID19/CODE/PROJECTS/Validation/paper") ### set to repository location

source("functions_read_data.R")
source("functions_analysis.R")

library(haven)
library(tidyverse)
library(survey)
library(lubridate)
library(matrixStats)
library(srvyr)
library(hrbrthemes)
library(ggpubr)
library(patchwork)

### read CSP data ####


waves <- c(14, 16:26)

datasets_csp <- read_CSP_data(waves)


### read and process CDC data ####

cdc <- read_preprocess_CDC_vax_data()


#### IPSOS Data ####


### read and process

out <- read_preprocess_IPSOS_data()

datasets_ipsos <- out[[1]]
dates_ipsos <- out[[2]]


waves_inc <- sapply(seq_along(datasets_ipsos), function(i) {
  if ("Q107_1" %in% colnames(datasets_ipsos[[i]])) {
    return(i)
  }
}) %>% unlist()

ipsos <- bind_rows(
  lapply(datasets_ipsos[waves_inc], function(df) df %>% select(WT_FINAL, WAVE, Q107_1))
) %>%
  left_join(dates_ipsos, by = "WAVE")


### calculate vax rate

ipsos_vax <- calculate_ipsos_vax(ipsos)



### National trends and diff ####


csp_vax <- calculate_csp_vax(datasets_csp, "weight") 


national_trends <- filter(cdc, Location == "US",
                          Date < "2023-02-01") %>%
  select(Date, pct_vax) %>%
  rename(vax = pct_vax) %>%
  mutate(`Data Source` = "CDC",
         moe = 0,
         vax = vax) %>%
  bind_rows(mutate(ipsos_vax, `Data Source` = "Ipsos"),
            mutate(csp_vax, `Data Source` = "CSP")) 


national_diff <- vax_national_diff(cdc, csp_vax, ipsos_vax)
  

plot_national_trends <- ggplot(national_trends) + 
  geom_pointrange(aes(x=Date, y=vax,
                      ymin = vax - moe,
                      ymax = vax + moe,
                      color = `Data Source`, 
                      size = `Data Source`)) +
  geom_line(aes(x = Date,
                y = vax,
                color = `Data Source`),
            lwd = 0.75) + 
  scale_color_manual(values=c('black', 'blue', 'darkred'))+
  scale_size_manual(values = c(0.01, 0.3, 0.3), guide = "none") +
  scale_y_continuous(breaks = c(0,20,40,60,80)) +
  geom_vline(xintercept = as.Date("2021-09-22"),
            linetype = "dashed",
            color = "darkolivegreen4") + 
  annotate("text", label =  "FDA approves booster shots", x = as.Date("2021-09-27"), y = 30, hjust=0) + 
  labs(x = "Date", y = "Percentage vaccinated") + 
  scale_x_date(date_breaks = "2 months",
               date_labels = "%b%y",
               limits = as.Date(c("2020-12-01","2023-01-31"))) +
  theme_ipsum()+
  theme(legend.position.inside=c(0.12, 0.85), 
    legend.justification='center',
    legend.direction='horizontal',
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
    axis.title.y = element_text(color='black',face='bold', size=14,
                                margin = margin(0,8,0,0)),
    axis.title.x = element_text(color='black',face='bold', size=14,
                                margin = margin(5,0,0,0)),
    legend.text = element_text(color='black', size=15),
    legend.title = element_blank(),
    legend.spacing.x = unit(0.6, 'cm'),
    panel.border = element_blank(),
    plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")) +
    guides(color = guide_legend(label.position = "bottom")) 


plot_national_diff <- ggplot(national_diff, aes(x = Date, y = diff, colour = Survey)) + 
  geom_pointrange(aes(ymin = (diff - moe), 
                      ymax = (diff + moe),
                      size = Survey)) +
  scale_size_manual(values = c(0.3,0.3)) +
  scale_color_manual(values=c( 'blue', 'darkred'))+
  geom_hline(yintercept = 0) +
  scale_x_date(breaks = "2 months", 
               date_labels = "%b%y",
               limits = as.Date(c("2020-12-01","2023-01-31"))) +
  geom_vline(xintercept = as.Date("2021-09-22"),
             linetype = "dashed",
             color = "darkolivegreen4") + 
  # annotate("text", label =  "FDA approves booster shots", 
  #          x = as.Date("2021-09-27"), 
  #          y = 5, hjust=0) +
  scale_y_continuous(minor_breaks = c(-25, -15, -5)) +
  xlab("Date") + 
  ylab("Difference between survey estimate \n and closest CDC estimate")+
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


## Comparing Ipsos and CSP estimates: SI Table

tab_comp <- national_trends %>% filter(`Data Source` == "CSP") %>% 
  left_join(filter(national_trends, `Data Source` == "Ipsos"), join_by(closest(Date <= Date))) %>%
  left_join(filter(national_trends, `Data Source` == "Ipsos"), join_by(closest(Date.x >= Date))) %>%
  mutate(days_diff1 = abs(Date.y - Date.x),
         days_diff2 = abs(Date - Date.x),
         days_diff = pmin(days_diff1, days_diff2, na.rm = T),
         vax_diff = case_when(days_diff1 >= days_diff2 | is.na(days_diff1) ~ vax - vax.x,
                              (days_diff2 > days_diff1)  ~ vax.y - vax.x),
         moe_ipsos =  case_when(days_diff1 >= days_diff2 | is.na(days_diff1) ~ moe,
                                days_diff2 > days_diff1  ~ moe.y),
         moe_total = moe_ipsos + moe.x) %>%
  rename(Date_csp = Date.x) %>%
  select(Date_csp, days_diff, vax_diff, moe.x, moe_ipsos, moe_total)

tab_comp %>% summarise(sum(abs(vax_diff)<=moe_total))

tab_comp %>% filter(Date_csp > "2021-09-07") %>%
  summarise(mean(abs(vax_diff)))


## Text numbers: comparing to CDC estimates

cdc_comp1 <- filter(national_diff, Date < as.Date("2021-09-22")) %>%
  group_by(Survey) %>%
  summarise(mean_diff = mean(abs(diff), na.rm = T),
            n_sig = sum(abs(diff) < moe, na.rm = T),
            n = sum(!is.na(diff))) 

cdc_comp2 <- filter(national_diff, Date >= as.Date("2021-09-22")) %>%
  group_by(Survey) %>%
  summarise(mean_diff = mean(abs(diff), na.rm = T),
            n_sig = sum(abs(diff) < moe, na.rm = T),
            n = sum(!is.na(diff))) 

## Text numbers: average estimates

national_trends %>% filter(`Data Source` == "CSP", Date > "2022-01-01") %>% 
  summarise(csp_mean = mean(vax),
            csp_sd = sd(vax))

national_trends %>% filter(`Data Source` == "Ipsos", Date > "2022-01-01") %>% 
  summarise(ipsos_mean = mean(vax),
            ipsos_sd = sd(vax))



### Text numbers: % with a second shot

lapply(datasets_csp[3:8], function(d) d %>%
         as_survey_design(weight = weight) %>%
         summarise(survey_mean(vac_2nd != 1, na.rm = T),
                   wave = first(wave)))

### Text numbers: ipsos differences above CI

ipsos_vax %>% 
  mutate(vax_l = vax - lag(vax)) %>%
  filter(Date > as.Date("2021-12-31"))
    



##### State analysis ####


csp_vax_state <- calculate_csp_vax_state(datasets_csp)
  

cdc_flt <- filter(cdc, Location != "US", 
                  Date < "2023-02-01") %>%
  select(Date, Location, pct_vax) %>%
  rename(cdc_vax = pct_vax) %>%
  mutate(cdc_vax = cdc_vax/100) %>%
  left_join(select(csp_vax_state, state_code, State) %>%
              unique, 
            by = c("Location" = "state_code"))


df_state_diff <- cdc_flt %>%
  full_join(csp_vax_state,
            by = c("Date" = "Date", 
                   "State" = "State")) %>%
  group_by(Location) %>%
  complete(Date = seq(min(Date), max(Date), by = "day")) %>%
  mutate(cdc_vax = zoo::na.fill(cdc_vax, "extend")) %>%
  ungroup() %>%
  mutate(diff =  csp_vax - cdc_vax)
  
df_state_preboost <- df_state_diff %>%
  filter(!is.na(diff), wave <= 19) %>%
  select(Location, Date, wave, moe, diff, State, n) %>% 
  mutate(`Difference within 95% CI` = case_when(moe >= abs(diff) ~ "Yes",
                                             moe < abs(diff) ~ "No"))

df_state_all <- df_state_diff %>%
  filter(!is.na(diff)) %>%
  select(Location, Date, wave, moe, diff, State, n) %>% 
  mutate(`Difference within 95% CI` = case_when(moe >= abs(diff) ~ "Yes",
                                             moe < abs(diff) ~ "No"),
         Period = case_when(wave <=19 ~ "Pre-boosters",
                            wave > 19 ~ "Post-boosters"))

  

### histogram


state_histogram <- ggplot(df_state_all, aes(x = 100*diff, y = 100*after_stat(density),
                         fill = Period, color = Period)) +
  geom_histogram(position = "identity", alpha = 0.4, 
                 key_glyph = "rect") +
  geom_vline(data = df_state_all %>% group_by(Period) %>% 
               summarise(median = median(100*diff)),
             aes(xintercept = median, color = Period),
             linetype = "dashed", 
             linewidth = 1) + 
  scale_color_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2") +
  xlab("Difference between CSP and CDC state estimate of percentage vaccinated") +
  ylab("Percentage of state-wave pairs") +
  theme_ipsum() +
  theme(
    legend.position=c(0.1, 0.7), 
    legend.justification='left',
    legend.direction='vertical',
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(color='black', size=12),
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
    axis.title.y = element_text(color='black',face='bold', size=15,
                                margin = margin(0,10,0,0)),
    axis.title.x = element_text(color='black',face='bold', size=15),
    panel.border = element_blank(),
    plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))


# MERGE PLOTS
full_plot <- plot_national_trends / plot_national_diff / state_histogram + 
  plot_annotation(tag_levels = list(c('a)', 'b)', 'c)')),
  ) &
  theme(plot.title = element_text(size = 16),
        plot.tag = element_text(face = 'bold',
                                #size = 8,
                                hjust = -5,
                                vjust = 5))

ggsave("Plots/vax_full_plot.png",
       full_plot, 
       height = 14, width = 10,
       bg = "white", dpi = 300)


## Text numbers:

df_state_all %>% 
  filter(Period == "Post-boosters") %>%
  summarise(sum(diff > -0.1),
            n())

df_state_all %>% 
  filter(Period == "Pre-boosters") %>%
  summarise(sum(diff > -0.05),
            n())


### SI state plot 1


ggplot(df_state_preboost,
       aes(x = Date, y = diff)) + 
  geom_pointrange(aes(ymin = diff - moe, 
                      ymax = diff + moe,
                      colour = `Difference within 95% CI`),
                   size = 0.2) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray")+
  facet_wrap(~State, nrow = 9) +
  scale_y_continuous(labels = scales::percent, 
                     limits = c(-0.2, 0.1),
                     oob = scales::squish) +
  scale_color_manual(values=c( 'red', 'black')) + 
  ylab("Difference between survey estimate and closest CDC estimate") + 
  xlab("Median date of survey wave") +
  theme_ipsum()+
  theme(legend.position='top', 
        legend.justification='center',
        legend.direction='horizontal',
        legend.box.background = element_rect(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(color='black',face='bold', size=13),
        axis.title.x = element_text(color='black',face='bold', size=13),
        axis.title.y.right = element_text(color='blue',face='bold', size=15),
        legend.text = element_text(color='black', size=12),
        panel.border = element_blank(),
        panel.spacing = unit(1, "lines")) +
  guides(color = guide_legend(nrow = 1))


ggsave("Plots/cdc_diff_state_panels.png",
       height = 14, width = 12,
       bg = "white", dpi = 300)


### state plot 2

ggplot() + 
  geom_pointrange(aes(x=Date, y=csp_vax,
                      ymin = csp_vax - moe,
                      ymax = csp_vax + moe,
                      color = 'CSP',
                      size = 'CSP'),
                  data = csp_vax_state) +
  geom_line(aes(x = Date,
                y = csp_vax, 
                color = 'CSP'),
            data = csp_vax_state,
            lwd = 0.75) + 
  geom_point(aes(x=Date, 
                 y=cdc_vax,
                 color = 'CDC',
                 size = 'CDC'),
             data = cdc_flt) +
  scale_color_manual(values=c('black','darkred'), guide = "none")+
  scale_size_manual(values = c(0.005, 0.15), guide = "none") +
  scale_y_continuous(labels = scales::percent, 
                     limits = c(0, 1.3),
                     oob = scales::squish) +
  facet_wrap(~State, nrow = 9) +
  geom_vline(xintercept = as.Date("2021-09-22"),
             linetype = "dashed",
             color = "darkolivegreen4") + 
  labs(x = "Date", y = "Percentage vaccinated") + 
  scale_x_date(date_breaks = "4 months",
               date_labels = "%b%y",
               limits = as.Date(c("2020-12-01","2023-01-31"))) +
  theme_ipsum()+
  theme(legend.position='top', 
        legend.justification='center',
        legend.direction='horizontal',
        legend.box.background = element_rect(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(color='black',face='bold', size=13),
        axis.title.x = element_text(color='black',face='bold', size=13),
        axis.title.y.right = element_text(color='blue',face='bold', size=15),
        legend.text = element_text(color='black', size=12),
        panel.border = element_blank(),
        panel.spacing = unit(1, "lines")) +
    guides(color = guide_legend(nrow = 1, title = NULL))

ggsave("Plots/state_panels.png",
       height = 14, width = 12,
       bg = "white", dpi = 300)


### Text numbers: summary tables

sum_tab <- df_state_preboost %>% group_by(wave) %>%
  summarise(`Number of states` = n(),
            `Differences within margin of error` = sum(moe >= abs(diff)),
            `Estimates within 10% of CDC estimate` = sum(abs(diff) < 0.1),
            `Estimates within 5% of CDC estimate` = sum(abs(diff) < 0.05),
            `Estimates within 2.5% of CDC estimate` = sum(abs(diff) < 0.025),
            `Estimates within 1% of CDC estimate` = sum(abs(diff) < 0.01))
colSums(sum_tab)


sum_tab_state <- df_state_preboost %>% 
  group_by(State) %>%
  summarise(`Differences within margin of error` = sum(moe >= abs(diff)),
            `Estimates within 10% of CDC estimate` = sum(abs(diff) < 0.1),
            `Estimates within 5% of CDC estimate` = sum(abs(diff) < 0.05),
            `Estimates within 2.5% of CDC estimate` = sum(abs(diff) < 0.025),
            `Estimates within 1% of CDC estimate` = sum(abs(diff) < 0.01))
colSums(sum_tab)

bad_states <- sum_tab_state %>% filter(`Estimates within 5% of CDC estimate` < 3)

state_samples <- lapply(datasets_csp[1:5], function(d) d %>% group_by(state) %>% 
                summarise(n = n())) %>% 
  bind_rows() %>%
  group_by(state) %>%
  summarise(n = mean(n))

sum_tab_state %>% filter(`Estimates within 5% of CDC estimate` < 3) %>% 
  left_join(state_samples, by = c("State" = "state")) %>%
  summarise(mean(n))

sum_tab_state %>% filter(`Estimates within 5% of CDC estimate` >= 3) %>% 
  left_join(state_samples, by = c("State" = "state")) %>%
  summarise(mean(n))

df_state_preboost %>% filter(wave > 14, `Difference within MOE` == "No") %>% 
  filter(diff<0)
  


##### SI table 2 with unadjusted CDC values

cdc_un <- read_csv("../COVID-19_Vaccinations_in_the_United_States_Jurisdiction.csv") %>%
  select(Date, Location, Administered_Dose1_Recip_18PlusPop_Pct) %>%
  filter(Date == "01/04/2023") %>%
  right_join(filter(cdc, Date == "2023-01-02") %>%
               select(Location, pct_vax) %>%
               rename(cdc_unadjusted = pct_vax)) %>%
  right_join(filter(csp_vax_state, wave == 26), by = c("Location" = "state_code")) %>%
  select(State, Administered_Dose1_Recip_18PlusPop_Pct, cdc_unadjusted, csp_vax, moe) %>%
  mutate(csp_vax = round(100*csp_vax),
         moe = round(100*moe),
         cdc_unadjusted = round(cdc_unadjusted),
         Administered_Dose1_Recip_18PlusPop_Pct = round(Administered_Dose1_Recip_18PlusPop_Pct)) %>%
  arrange(State)



##### Vax by Trust: Extended Data Figure 9 ####

vax_trust_cdc <- lapply(
  datasets_csp[1:10], function(d)  d %>%
    mutate(EndDate = as.Date(EndDate),
           StartDate = as.Date(StartDate),
           days = EndDate - min(StartDate),
           Date = min(StartDate) + 
             weightedMedian(days, weight, ties = "weighted")) %>%
    as_survey_design(weight = weight) %>%
    group_by(cov_trust_cdc) %>%
    summarise(vax = 100*survey_mean(vaccine_1, na.rm = T), 
              Date = first(Date), 
              Wave = first(wave)) %>%
    mutate(moe = vax_se*qnorm(.025, lower.tail = F)) %>%
    select(-vax_se)
) %>%
  bind_rows() %>%
  filter(!is.na(cov_trust_cdc)) %>%
  mutate(cov_trust_cdc = as_factor(cov_trust_cdc))

vax_trust_science <- lapply(
  datasets_csp[1:10], function(d)  d %>%
    mutate(EndDate = as.Date(EndDate),
           StartDate = as.Date(StartDate),
           days = EndDate - min(StartDate),
           Date = min(StartDate) + 
             weightedMedian(days, weight, ties = "weighted")) %>%
    as_survey_design(weight = weight) %>%
    group_by(cov_trust_science) %>%
    summarise(vax = 100*survey_mean(vaccine_1, na.rm = T), 
              Date = first(Date), 
              Wave = first(wave)) %>%
    mutate(moe = vax_se*qnorm(.025, lower.tail = F)) %>%
    select(-vax_se)
) %>%
  bind_rows() %>%
  filter(!is.na(cov_trust_science)) %>%
  mutate(cov_trust_science = as_factor(cov_trust_science))


trust_science <- ggplot(vax_trust_science, 
                    aes(x = Date,
                        y = vax,
                        group = cov_trust_science)) +
  geom_pointrange(aes(ymin = (vax - moe), 
                      ymax = (vax + moe),
                      color = cov_trust_science),
                  position = position_dodge2(width = 10)) +
  xlab("Median date of survey wave") + 
  ylab("Vaccination rate") +
  scale_color_viridis_d(name = "Trust in scientists to handle \n the COVID-19 pandemic") +
  scale_y_continuous(breaks = c(0,20,40,60,80)) + 
  theme_bw() +
  theme(legend.position = "bottom")
  
  
trust_cdc <- ggplot(vax_trust_cdc, 
                          aes(x = Date,
                              y = vax,
                              group = cov_trust_cdc)) +
    geom_pointrange(aes(ymin = (vax - moe), 
                        ymax = (vax + moe),
                        color = cov_trust_cdc),
                    position = position_dodge2(width = 10)) +
    xlab("Median date of survey wave") + 
    ylab("Vaccination rate") +
    scale_color_viridis_d(name = "Trust in the CDC to handle \n the COVID-19 pandemic") +
    scale_y_continuous(breaks = c(0,20,40,60,80)) + 
    theme_bw() +
    theme(legend.position = "bottom")
  
trust_bw_plot <- trust_cdc + trust_science


ggsave("PROJECTS/Validation/paper/Plots/trust_bw_vax.png",
       trust_bw_plot,
       height = 7, width = 12,
       bg = "white", dpi = 300)


