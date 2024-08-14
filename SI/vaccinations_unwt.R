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

#### set weights to 1 for all responses (equivalent to unweighted analysis)
datasets_csp <- read_CSP_data(waves) %>% 
  lapply(., function(d) d %>% mutate(weight = 1,
                                     weight_state = 1))



##### RUN CODE AS IN vaccinations.R #####

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

ggsave("Plots/vax_full_plot_unwt.png",
       full_plot, 
       height = 14, width = 10,
       bg = "white", dpi = 300)

