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

datasets_csp <- read_CSP_data(waves)

datasets_csp_ps <- lapply(datasets_csp, function(d) filter(d, source == "PS"))
datasets_csp_vs <- lapply(datasets_csp, function(d) filter(d, source == "VS"))


lapply(datasets_csp, function(d) table(d$source)[2]) %>% unlist()
table(datasets_csp[[1]]$source)

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


csp_vax_ps <- calculate_csp_vax(datasets_csp_ps, "weight") 
csp_vax_vs <- calculate_csp_vax(datasets_csp_vs, "weight") 


national_trends <- filter(cdc, Location == "US",
                          Date < "2023-02-01") %>%
  select(Date, pct_vax) %>%
  rename(vax = pct_vax) %>%
  mutate(`Data Source` = "CDC",
         moe = 0,
         vax = vax) %>%
  bind_rows(mutate(ipsos_vax, `Data Source` = "Ipsos"),
            mutate(csp_vax_ps, `Data Source` = "CSP Pure Spectrum"),
            mutate(csp_vax_vs, `Data Source` = "CSP Facebook")) 


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
  scale_color_manual(values=c('black', 'steelblue1', 'darkorchid1', 'darkred'))+
  scale_size_manual(values = c(0.01, 0.3, 0.3, 0.3), guide = "none") +
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
  theme(legend.position=c(0.6, 0.2), 
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


ggsave("PROJECTS/Validation/paper/Plots/vaccination_rate_VS_PS.png", 
       plot_national_trends, 
       height = 6, width = 9,
       bg = "white", dpi = 300)



#### regressions ####

library(marginaleffects)

trust_reg <- function(dat) {
  dat <- mutate(dat,
                cov_trust_science = zap_missing(cov_trust_science))
  mod <- svyglm(cov_trust_science ~ source + factor(education) + factor(race) + factor(age_cat_6) +
                  factor(urbanicity) + factor(gender) + factor(income_cat_10),
                as_survey_design(dat, weights = weight), family = gaussian())
  return(list(mod, mod$coefficients, confint(mod)))
}

pid_reg <- function(dat) {
  mod <- svyglm(party7 ~ source + factor(education) + factor(race) + factor(age_cat_6) +
                  factor(urbanicity) + factor(gender) + factor(income_cat_10),
                as_survey_design(dat, weights = weight), family = gaussian())
  return(list(mod, mod$coefficients, confint(mod)))
}

dem_reg <- function(dat) {
  mod <- svyglm(party3 == "Democrat" ~ source + factor(education) + factor(race) + factor(age_cat_6) +
                  factor(urbanicity) + factor(gender) + factor(income_cat_10),
                as_survey_design(dat, weights = weight), family = gaussian())
  return(list(mod, mod$coefficients, confint(mod)))
}

trust_reg_no_control <- function(dat) {
  dat <- mutate(dat,
                trust_science = case_when(cov_trust_science == 4 ~ 1, 
                                          cov_trust_science %in% c(1,2,3) ~ 0, 
                                          T ~ NA_real_),
                cov_trust_science = zap_missing(cov_trust_science))
  mod <- svyglm(cov_trust_science ~ source,
                as_survey_design(dat, weights = weight), family = gaussian())
  return(c(mod$coefficients[[2]], confint(mod)[2,]))
}

trust_cdc_reg_no_control <- function(dat) {
  mod <- svyglm(cov_trust_cdc ~ source,
                as_survey_design(dat, weights = weight), family = gaussian(),
                cov_trust_cdc = zap_missing(cov_trust_cdc))
  return(c(mod$coefficients[[2]], confint(mod)[2,]))
}

vax_reg <- function(dat) {
  dat <- mutate(dat, 
                days = as.Date(StartDate) - min(as.Date(StartDate)))
  mod <- svyglm(vaccine_1 ~ source + days + factor(education) +  factor(race) + factor(age_cat_8) +
                  factor(urbanicity) + factor(gender) + factor(income_cat_10), 
                as_survey_design(dat, weights = weight), 
                family = quasibinomial())
  return(list(mod, mod$coefficients, confint(mod)))
}

vax_reg_trust <- function(dat) {
  dat <- mutate(dat, 
                days = as.Date(StartDate) - min(as.Date(StartDate)),
                cov_trust_science = zap_missing(cov_trust_science))
  mod <- svyglm(vaccine_1 ~ source + days + factor(education) + factor(race) + factor(age_cat_8) +
                  factor(urbanicity) + factor(cov_trust_science) + factor(gender) + factor(income_cat_10), 
                as_survey_design(dat, weights = weight), 
                family = quasibinomial())
  return(list(mod, mod$coefficients, confint(mod)))
}

vax_reg_pid <- function(dat) {
  dat <- mutate(dat, 
                days = as.Date(StartDate) - min(as.Date(StartDate)),
                cov_trust_science = zap_missing(cov_trust_science))
  mod <- svyglm(vaccine_1 ~ source + days + factor(education) + factor(race) + factor(age_cat_8) +
                  factor(urbanicity) + factor(party3) + factor(gender) + factor(income_cat_10), 
                as_survey_design(dat, weights = weight), 
                family = quasibinomial())
  return(list(mod, mod$coefficients, confint(mod)))
}

vax_reg_pid_trust <- function(dat) {
  dat <- mutate(dat, 
                days = as.Date(StartDate) - min(as.Date(StartDate)))
  mod <- svyglm(vaccine_1 ~ source + days + factor(education) +  factor(race) + factor(age_cat_8) +
                  factor(urbanicity) + factor(gender) + factor(income_cat_10) + factor(party3) + 
                  factor(cov_trust_science), 
                as_survey_design(dat, weights = weight), 
                family = quasibinomial())
  return(list(mod, mod$coefficients, confint(mod)))
}



vax_reg_no_controls <- function(dat) {
  mod <- svyglm(vaccine_1 ~ source, 
                as_survey_design(dat, weights = weight), 
                family = quasibinomial())
  return(c(mod$coefficients[[2]], confint(mod)[2,]))
}

vax_reg_no_controls_trust <- function(dat) {
  mod <- svyglm(vaccine_1 ~ source + factor(cov_trust_science), 
                as_survey_design(dat, weights = weight), 
                family = quasibinomial())
  return(c(mod$coefficients[[2]], confint(mod)[2,]))
}

regs_to_df <- function(regs, idxs = c(1:10)) {
  df <- reduce(regs, rbind) %>%
    as_tibble() %>%
    mutate(across(everything(), ~ round(., 2)),
           Wave = sapply(datasets_csp[idxs], function(d) d$wave[[1]]), 
           Estimate = str_c(V1, " (", `2.5 %`, " - ", `97.5 %`, ")")) %>%
    select(Wave, Estimate)

  return(df)
}



trust_regs <- lapply(datasets_csp[1:10], trust_reg)
trust_regs_df <- lapply(trust_regs, function(x) c(x[[2]][[2]], x[[3]][2,])) %>%
  regs_to_df()

pid_regs <- lapply(datasets_csp[1:10], pid_reg) 
pid_regs_df <- lapply(pid_regs, function(x) c(x[[2]][[2]], x[[3]][2,])) %>%
regs_to_df()

#   
# trust_df_no_control <- lapply(datasets_csp[1:10], trust_reg_no_control) %>%
#   regs_to_df(c(1:10))
# 
# trust_cdc_df_no_control <- lapply(datasets_csp[1:10], trust_cdc_reg_no_control) %>%
#   regs_to_df(c(1:10))


# vax 
vax_regs <- lapply(datasets_csp[2:10], vax_reg) 
vax_df <- lapply(vax_regs, function(x) c(x[[2]][[2]], x[[3]][2,])) %>%
  regs_to_df(c(2:10))

# counterfactual_dat <- datagrid(newdata = datasets_csp[[2]],
#                                source = labelled_spss(c("PS", "VS")),
#                                grid_type = "counterfactual")

preds_counterfactual <- 
  lapply(vax_regs, function(x) {
    predictions(x[[1]],
                type = "response",
                by = "source",
                newdata = datagrid(source = labelled_spss(c("PS", "VS")),
                                   grid_type = "counterfactual"))
  }
  )

pred_vax <- lapply(preds_counterfactual, function(d) d %>% select(estimate) %>% t()) %>%
  reduce(rbind)


# trust #
vax_regs_trust <- lapply(datasets_csp[2:10], vax_reg_trust) 
trust_coefs <- lapply(vax_regs_trust, function(reg) reg[[2]][27:29]) %>%
  bind_rows()
trust_CI <- lapply(vax_regs_trust, function(reg) reg[[3]][27:29,]) 

vax_df_trust <- lapply(vax_regs_trust, function(x) c(x[[2]][[2]], x[[3]][2,])) %>%
                         regs_to_df(c(2:10))

# preds_counterfactual_trust <- 
#   lapply(vax_regs_trust, function(x)
#     predictions(x[[1]],
#             type = "response",
#             by = "source",
#             newdata = datagrid(source = labelled_spss(c("PS", "VS")),
#                                grid_type = "counterfactual"))
#   )
# 
# pred_vax_trust <- lapply(preds_counterfactual_trust,
#                          function(d) d %>% select(estimate) %>% t()) %>%
#   reduce(rbind)

## trust + pid 

vax_regs_pid <- lapply(datasets_csp[2:10], vax_reg_pid) 
vax_df_pid<- lapply(vax_regs_pid, function(x) c(x[[2]][[2]], x[[3]][2,])) %>%
  regs_to_df(c(2:10))

vax_regs_pid_trust <- lapply(datasets_csp[2:10], vax_reg_pid_trust) 
vax_df_pid_trust <- lapply(vax_regs_pid_trust, function(x) c(x[[2]][[2]], x[[3]][2,])) %>%
  regs_to_df(c(2:10))


vax_df_no_controls <- lapply(datasets_csp[2:10], vax_reg_no_controls) %>%
  regs_to_df(c(2:10))

vax_df_no_controls_trust <- lapply(datasets_csp[2:10], vax_reg_no_controls_trust) %>%
  regs_to_df(c(2:10))


vs_regs <- lapply(datasets_csp[2:10], vs_reg_pid_trust) 
