library(tidyverse)
library(srvyr)
library(survey)
library(matrixStats)
library(data.table)


calc_err <- function(pct, wt) {
  if(is.list(wt)) wt <- unlist(wt)
  n <- length(wt)
  deff <- n*(sum(wt^2)/sum(wt)^2)
  sqrt(deff)*sqrt((pct*(1-pct))/(n-1))
}

calculate_ipsos_vax <- function(ipsos_data) {
  
  ipsos_vax <- ipsos_data %>% group_by(WAVE) %>%
    group_map(~ as_survey_design(., weight = WT_FINAL) %>%
                summarise(vax = survey_mean(Q107_1 == 1, na.rm = T),
                          Date = first(Date),
                          Wave = first(WAVE),
                          wt = list(WT_FINAL)),
              .keep = T) %>%
    bind_rows() %>%
    rowwise() %>%
    mutate(moe = calc_err(vax, wt)*100*qnorm(.025, lower.tail = F),
           vax = 100*vax) %>%
    select(-wt, -Wave, -vax_se) %>%
    ungroup()
  
  return(ipsos_vax)
}


calculate_csp_vax <- function(datasets, weight_var = "weight") {
  lapply(
    datasets, function(d)  d %>%
      mutate(EndDate = as.Date(EndDate),
             StartDate = as.Date(StartDate),
             days = EndDate - min(StartDate),
             Date = min(StartDate) + 
               weightedMedian(days, !!sym(weight_var), ties = "weighted")) %>%
      as_survey_design(weight = !!sym(weight_var)) %>%
      summarise(vax = 100*survey_mean(vaccine_1, na.rm = T), 
                Date = first(Date), 
                Wave = first(wave)) %>%
      mutate(moe = vax_se*qnorm(.025, lower.tail = F)) %>%
      select(-vax_se)
  ) %>%
    bind_rows()
}

vax_national_diff <- function(cdc, csp_vax, ipsos_vax) {
  
  filter(cdc, Location == "US", 
                          Date < "2023-02-01") %>%
    select(Date, pct_vax) %>%
    rename(cdc_vax = pct_vax) %>%
    full_join(csp_vax %>% 
                rename(csp_moe = moe,
                       csp_vax = vax), 
              by = "Date") %>%
    full_join(ipsos_vax %>% 
                rename(ipsos_vax = vax,
                       ipsos_moe = moe), 
              by = "Date") %>%
    complete(Date = seq(min(Date), max(Date), by = "day")) %>%
    mutate(cdc_vax = zoo::na.fill(cdc_vax, "extend")) %>%
    mutate(csp_diff =  csp_vax - cdc_vax,
           ipsos_diff = ipsos_vax - cdc_vax) %>%
    filter(!is.na(ipsos_diff) | !is.na(csp_diff)) %>%
    select(-cdc_vax, -ipsos_vax, -csp_vax, -Wave) %>%
    pivot_longer(cols = -Date, 
                 names_sep = "_",
                 names_to = c("Survey", "val")) %>%
    pivot_wider(names_from = val, 
                values_from = value) %>%
    mutate(Survey = ifelse(Survey == "csp", "CSP", "Ipsos"))
}

vax_diff_downsample <- function(cdc, csp_vax) {
  
  filter(cdc, Location == "US", 
         Date < "2023-02-01") %>%
    select(Date, pct_vax) %>%
    rename(cdc_vax = pct_vax) %>%
    full_join(csp_vax, 
              by = "Date") %>%
    complete(Date = seq(min(Date), max(Date), by = "day")) %>%
    mutate(cdc_vax = zoo::na.fill(cdc_vax, "extend")) %>%
    mutate(csp_diff =  vax - cdc_vax) %>%
    filter(!is.na(csp_diff)) 
}

calculate_csp_vax_state <- function(datasets, weight_var = "weight_state") {
  lapply(
    datasets, function(d)  d %>%
      mutate(EndDate = as.Date(EndDate),
             StartDate = as.Date(StartDate)) %>%
      as_survey_design(weight = !!sym(weight_var)) %>%
      group_by(state_code) %>%
      summarise(csp_vax = survey_mean(vaccine_1, na.rm = T),
                Date = min(StartDate) + 
                  weightedMedian(EndDate - min(StartDate), !!sym(weight_var), ties = "weighted"),
                wave = first(wave),
                State = first(state),
                n = n())
  ) %>%
    bind_rows() %>%
    mutate(moe = csp_vax_se*qnorm(.025, lower.tail = F)) %>%
    select(-csp_vax_se) 
}


calculate_ipsos_infections <- function(datasets_ipsos_n_inf) {
  
  datasets_ipsos_n_inf %>% 
    group_by(WAVE) %>%
    group_map( ~ as_survey_design(., weight = WT_FINAL) %>%
                 summarise(avg_infections = survey_mean(n_infections, na.rm = T),
                           Date = first(Date), 
                           Wave = first(WAVE),
                           wt = list(WT_FINAL)),
               .keep = T) %>%
    bind_rows()  %>%
    rowwise() %>%
    mutate(std_error = calc_err(avg_infections, wt)*100,
           avg_infections = 100*avg_infections) %>%
    select(-wt, -avg_infections_se) %>%
    ungroup()
  
}


read_wave_month_table <- function(path) {
  wave_month_match_table <- fread(str_c(path, "infections_months.csv"))
  wave_month_match_table$months <- as.Date(wave_month_match_table$months , format = "%m/%d/%y")
  wave_month_match_table$wave_date <- as.Date(wave_month_match_table$wave_date , format = "%m/%d/%y")
  wave_month_match_table <- data.frame(wave_month_match_table)
  
  months_per_wave <- wave_month_match_table %>% 
    group_by(wave_x) %>%
    summarise(month_min = ymd(min(months)),
              month_max = ymd(max(months))) %>%
    rename(wave = wave_x)
  
  return(list(wave_month_match_table, months_per_wave))
}

## given a wave dataset and a list of cov_month variables as strings, generate the corresponding
## eligible variable that's 1 when the respondent got covid in any of these months. 
generate_infection_var <- function(dat, vars) {
  
  ### create the case_when argument as string from the provided cov_month vars
  vars <- paste0("as.integer(", vars, ")")
  pattern <- "(cov_test == 1) & ("
  pattern <- paste0(pattern, 
                    paste(vars, collapse = " | "),
                    ") ~ 1"
  )
  
  ### generate the eligible variable by evaluating the string within case_when
  dat <- mutate(dat, eligible = case_when(eval(rlang::parse_expr(pattern)),
                                          T ~ 0))
  return(dat)
}

get_datasets_method1 <- function(datasets_csp, wave_month_match_table) {
  
  ### use the wave_month_match_table to create a list with the cov_month variables
  ### corresponding to each wave. 
  month_var_blocks <- lapply(datasets_csp, function(d) 
    filter(wave_month_match_table, wave_x == d$wave[1]) %>%
      pull(variable))
  
  ### iterate over the survey datasets and run 'generate_infection_var' on each of them
  ### when they have at least one cov_month var (some waves are not used)
  datasets_csp_method1 <- lapply(seq_along(datasets_csp), function(i) {
    month_vars <- month_var_blocks[[i]]
    
    if(length(month_vars > 0)) {
      return(generate_infection_var(datasets_csp[[i]],
                                    month_vars))
    }
    
    else(return(NULL))
    
  })
  
  
  ### remove the datasets that are not used from the list of datasets. 
  datasets_csp_method1 <- Filter(Negate(is.null), datasets_csp_method1)
  
  return(datasets_csp_method1)
  
}

### function to calculate the percentage of infected within the block of months
### corresponding to each dataset
calculate_infected_method1 <- function(dat, weight_var = "weight"){
  dat %>% as_survey_design(ids = 1, weights = !!sym(weight_var)) %>%
    summarize(wave= first(wave),
              wave_end_date = first(wave_date),
              estimate = 100*survey_mean(eligible, na.rm = T)) %>%
    rename(std_error = estimate_se)
} 

### function to calculate the percentage of infected within the block of months
### corresponding to each dataset BY STATE
calculate_infected_method1_states <- function(dat){
  dat %>% as_survey_design(ids = 1, weights = weight_state) %>%
    group_by(state_code) %>%
    summarize(state = first(state),
              wave = first(wave),
              wave_end_date = first(wave_date),
              estimate = 100*survey_mean(eligible, na.rm = T)) %>%
    rename(std_error = estimate_se)
} 

generate_var_method2 <- function(datasets_csp) {
  
  ### function to count number of unique non-consecutive infections from the months
  ### a respondent was infected as a string of 0's and 1's
  how_many_blocks2 <- function(string){
    count = 0
    if(substr(string, 1,1) == '1'){count = count+1}
    for(i in 1:(nchar(string)-2)){ ### minus 1 to include the last cov_month, minus 2 to discard the last one
      if(substr(string, i, i+1) == '01'){count = count+1}
    }
    return(count)
  }
  
  how_many_blocks2 <- Vectorize(how_many_blocks2, USE.NAMES = F)
  
  ### apply 'how_many_blocks2' to generate a new variable with the number of infections
  ### of each respondent. Only count cases when cov_test is equal to 1.
  #datasets_csp_backup <- datasets_csp
  datasets_csp <- lapply(datasets_csp, function (dat) {
    dat %>%
      mutate(across(starts_with("cov_month"), ~replace_na(., 0))) %>%
      unite(block_infections, starts_with("cov_month"), sep = "") %>%
      mutate(count_cases = case_when(
        cov_test != 1 ~ 0, 
        cov_test == 1 ~ how_many_blocks2(block_infections)
        #,T ~ NA_integer_
      )
      )
  })
  return(datasets_csp)
}

calculate_infected_method2 <- function(dat, weight_var="weight"){
  dat %>% as_survey_design(ids = 1, weights = !!sym(weight_var)) %>%
    summarize(wave = first(wave),
              weighted_median_date = survey_median(EndDate_numeric, na.rm = T),
              wave_end_date = as.Date(first(wave_date)),
              estimate = 100*survey_mean(count_cases, na.rm = T)) %>%
    rename(std_error = estimate_se)
} 

calculate_infected_method2_states <- function(dat){
  dat %>% as_survey_design(ids = 1, weights = weight) %>%
    group_by(state_code) %>%
    summarize(state = first(state),
              wave = first(wave),
              wave_end_date = as.Date(first(wave_date)),
              weighted_median_date = survey_median(EndDate_numeric, na.rm = T),
              estimate = 100*survey_mean(count_cases, na.rm = T),
              n_respondents = n()) %>%
    rename(std_error = estimate_se)
} 

calculate_infections_diff <- function(us_cumulative_nyt, method1_results,
                                      method2_results, n_infections_ipsos) {
  us_cumulative_nyt %>%
    rename(nyt_pct_cumulative_cases = pct_cumulative_cases) %>%
    full_join(select(method1_results, date_plot, cdf_estimate, cdf_std) %>%
                rename(method1_pct_cumulative_cases = cdf_estimate,
                       method1_std = cdf_std,
                       date = date_plot), 
              by = "date")  %>%
    full_join(method2_results %>%
                rename(method2_pct_cumulative_cases = estimate,
                       method2_std = std_error,
                       date = weighted_median_date), 
              by = "date") %>% 
    full_join(n_infections_ipsos %>% 
                rename(ipsos_pct_cumulative_cases = avg_infections,
                       ipsos_std = std_error, 
                       date = Date), 
              by = "date") %>%
    mutate(method1_diff =  method1_pct_cumulative_cases - nyt_pct_cumulative_cases,
           method2_diff = method2_pct_cumulative_cases - nyt_pct_cumulative_cases,
           ipsos_diff = ipsos_pct_cumulative_cases - nyt_pct_cumulative_cases) %>%
    filter(!is.na(method1_diff) | !is.na(method2_diff) | !is.na(ipsos_diff)) %>%
    select(date, method1_diff, method2_diff, ipsos_diff, method1_std, method2_std, ipsos_std)
}

downsample_csp_dataset_random <- function(dat, prop_keep) {
  # n_keep <- dat %>%
  #   group_by(cov_trust_cdc) %>%
  #   summarise(n = n()) %>%
  #   left_join(fraction_keep) %>%
  #   mutate(n_keep = fraction_keep*n) %>%
  #   summarise(sum = sum(n_keep)) %>%
  #   pull(sum) 
  # 
  dat_random <- slice_sample(dat, prop = prop_keep)
  return(dat_random)
}

downsample_csp_dataset_trust <- function(dat, fraction_keep ){

  dat_trust <- dat %>%
    group_by(cov_trust_cdc) %>%
    group_modify(~ slice_sample(.x, 
                                prop = fraction_keep %>% 
                                  filter(cov_trust_cdc == unique(.x$cov_trust_cdc)) %>%
                                  pull(fraction_keep)) %>%
                   select(-cov_trust_cdc),
                 .keep = T) %>%
    ungroup
  

    return(dat_trust)
}

downsample_csp_dataset_state <- function(dat, population_perc ){
  
  sampling <- dat %>%
    group_by(state_code) %>% 
    summarise(n= n()) %>%
    mutate(freq = n/sum(n)) %>%
    select(-n) %>%
    left_join(population_perc) %>%
    mutate(repr = Population_perc/freq,
           undersampling = repr/max(repr))
  
  dat <- dat %>%
    group_by(state_code) %>%
    group_map(~ slice_sample(.x, 
                                prop = sampling %>% 
                                  filter(state_code == unique(.x$state_code)) %>%
                                  pull(undersampling)),
                 .keep = T) %>%
    bind_rows()
  
  
  return(dat)
}

downsample_csp_trust <- function(datasets, fraction_keep = tibble(cov_trust_cdc = c(1,2,3,4, -99), 
                                                           fraction_keep=c(0.7,0.8,0.9,1, 1))) {
  
  datasets_trust <- lapply(datasets, downsample_csp_dataset_trust, fraction_keep = fraction_keep)
  
  fraction_kept <- sapply(seq_along(datasets_trust), 
                                      function(i) nrow(datasets_trust[[i]])/nrow(datasets[[i]]))
  
  datasets_trust <- lapply(datasets_trust, function(d) mutate(d, weight = weight_csp(dat = d)))
  
  return(list(datasets_trust, fraction_kept))
  
}

downsample_csp_rand <- function(datasets, prop_keep) {
  
  datasets_rand <- lapply(datasets, downsample_csp_dataset_random, prop_keep)
  
  datasets_rand <- lapply(datasets_rand, function(d) mutate(d, weight = weight_csp(dat = d)))
  
  return(datasets_rand)
  
}

downsample_csp_state <- function(datasets, population_perc) {
  
  datasets_state <- lapply(datasets, downsample_csp_dataset_state, population_perc)
  
  fraction_kept <- sapply(seq_along(datasets_state), 
                          function(i) nrow(datasets_state[[i]])/nrow(datasets[[i]]))
  
  datasets_state <- lapply(datasets_state, function(d) mutate(d, weight = weight_csp(dat = d)))
  
  return(list(datasets_state, fraction_kept))
  
}

# --------------------- Generate weights ------------------

# Source this R file and use the 'weight_csp()' function to generate weights.
# Weighted data has to include the standard COVID States Project demographic variables.
# Variables used for weighing cannot contain missing data. Drop the NAs or impute first.

# PARAMETERS for weight_cps():
#
# dat:          Data used to create the weights.
# wave:         Survey wave. If 'dat' is NULL, the function will attempt to read a wave file
#                 using the COVID States naming convention from the directory 'data_folder'.
# data_folder:  The directory where the wave data is expected to be stored.
# state:        If FALSE, generate national weights. If TRUE, generate state weights.
# interlock:    If TRUE and "state" is FALSE, use interlocking variables for national weights.
# trim:         If TRUE the weights get trimmed at .01 and .99 percentile.
# party:        If FALSE, do not use party when generating weights.
#               If party==1, use Pew 2021 NPORS Dataset to weigh by party (only national weights)
#               If party==2, use Pew 2021 NPORS Dataset to weigh by party x gender (only national)
#               If party==3, use FEC election data to weigh by turnout and 2020 vote (national or state)
#               (If party is anything but FALSE and state is TRUE, 2020 election weights are used.)


weight_csp <- function(wave=24, dat=NULL, state=FALSE, interlock=TRUE, 
                       party=FALSE, trim=TRUE, data_folder="../DATA/CLEAN/RDS/") {
  
  library(survey)
  
  # Function testing variables:
  # wave <- '23'
  # data_folder="../DATA/CLEAN/RDS/"
  # dat=NULL
  # state=FALSE
  # interlock=TRUE
  # party=3
  # 
  
  # ~~~~~~~~~~~~~ Check data and variables ~~~~~~~~~~~~~~~
  
  # Missing dataset or wave:
  if(is.null(dat)) {
    
    if(is.na(wave)) {
      stop("Bad news: Missing both wave number and data object.") }
    
    file_name <- paste0(data_folder, "CSP_W", wave, ".RDS")
    
    if(!file.exists(file_name)) {
      stop(paste("Bad news: Missing data file", file_name)) } 
    
    dat <- readRDS(file_name)  
  }
  
  # Key demographics
  demo_vars <- c("gender", "race", "age", "education", 
                 "region", "state", "state_code", "urbanicity")
  
  # Missing demographic, party, or vote variables:
  
  if (!all(demo_vars %in% names(dat))) {
    stop("Bad news: The data is missing key demographic variables.") }
  
  if (party %in% 1:2 & !"party" %in% names(dat))   {
    stop("Bad news: The data is missing a 'party' variable.") }
  
  if (party %in% 3   & !"voted20" %in% names(dat)){
    stop("Bad news: The data is missing a 2020 election variable.") }
  
  
  keep_vars <- c(demo_vars, "state_code", "party", "voted20")
  keep_vars <- intersect(keep_vars, names(dat))
  dat <- dat[, keep_vars]  
  
  
  # ~~~~~~~~~~~~~ Get population & sample data ~~~~~~~~~~~~~~~
  
  # Get population and sample data:
  sample_n <- nrow(dat) # Sample size
  
  # Population percentage data:
  pop <- weight_pop()   
  
  # Generate weighting variables:
  dat <- weight_vars(dat, party=party) 
  
  # Switch the population percentages to expected Ns for our sample:
  pop <- lapply(pop, function(x) {x$Freq <- sample_n * x$Freq; x})
  
  
  # ~~~~~~~~~~~~~ Check interlocking categories ~~~~~~~~~~~~~~~
  
  # Remove any interlocking categories missing from our data. 
  # Ideally there should be none of these!
  # NOTE! if the dataset is too small, doing this will throw off the interlocking categories.
  
  include_interlock     <- pop$race_sex_age_dist$race_sex_age %in% dat$race_sex_age
  pop$race_sex_age_dist <- pop$race_sex_age_dist[include_interlock,]
  
  include_interlock2 <- pop$edu_age_dist$edu_age %in% dat$edu_age
  pop$edu_age_dist   <- pop$edu_age_dist[include_interlock2,]
  
  # Print out a message warning about missing categories
  if((any(!include_interlock) | any(!include_interlock2))
     & interlock %in% TRUE & state %in% FALSE) {
    message( paste( "Some interlocking categories have 0 cases. \n See "),
             paste(pop$race_sex_age_dist$race_sex_age[!include_interlock], collapse=" "),
             " ", paste( pop$edu_age_dist$edu_age[!include_interlock2]), collapse=" ",
             "\n Proceed with caution! Consider non-interlocked weights.")    }
  
  
  # ~~~~~~~~~~~~~ party variable options ~~~~~~~~~~~~~~~
  
  # Pew 2021 NPORS Dataset - Party
  if(party==1) { 
    dat$party_w    <- dat$party_pid 
    pop$party_dist <- pop$party_pid_dist }
  
  # Pew 2021 NPORS Dataset - Party by gender
  if(party==2) { 
    dat$party_w    <- dat$party_gender
    pop$party_dist <- pop$party_gender_dist }
  
  # Voting 2020 from FEC.gov data
  if(party==3) { 
    dat$party_w    <- dat$party_vote20
    pop$party_dist <- pop$party_vote20_dist }
  
  
  # ~~~~~~~~~~~~~ National weights ~~~~~~~~~~~~~~~
  
  if(!state & interlock & sample_n<5000) {
    warning("NOTE: Sample may be too small for interlocking weights!") }
  
  if(!state & !interlock & !party) {
    
    sample_margins <- list("gender", "region", "race_cat", 
                           "age_cat", "edu_cat", "urban_cat")
    pop_margins    <- pop[c("gender_dist", "region_dist", "race_dist", 
                            "age_dist", "edu_dist", "urban_dist")]  
  }
  
  if(!state & !interlock & party) {
    
    sample_margins <- list("gender", "region", "race_cat", "age_cat", 
                           "edu_cat", "urban_cat", "party_w")
    pop_margins    <- pop[c("gender_dist", "region_dist", "race_dist",  "age_dist",
                            "edu_dist", "urban_dist", "party_dist")]  
  }
  
  if(!state & interlock & !party)  {
    
    sample_margins <- list("race_sex_age", "edu_age",
                           "region", "urban_cat")
    pop_margins   <- pop[c("race_sex_age_dist", "edu_age_dist", 
                           "region_dist", "urban_dist")]  
  }
  
  if(!state & interlock & party) {
    
    sample_margins <- c("race_sex_age", "edu_age",
                        "region", "urban_cat", "party_w")
    pop_margins   <- pop[c("race_sex_age_dist", "edu_age_dist", 
                           "region_dist", "urban_dist", "party_dist")]  
  }
  
  # NATIONAL WEIGHTS:
  
  if(!state) {
    
    sample_margins <- paste0("~", sample_margins)
    sample_margins <- lapply(sample_margins, as.formula)
    
    dat_survey   <- svydesign(ids=~1, weights=~1, data=dat)
    dat_survey_w <-  rake(design = dat_survey,
                          sample.margins = sample_margins,
                          population.margins = pop_margins)
    weights(dat_survey_w)
    summary(weights(dat_survey_w))
    
    # Trimming the weights (top and bottom 1%):
    if(trim) { 
      
      low_trim  <- quantile(weights(dat_survey_w), .01) 
      high_trim <- quantile(weights(dat_survey_w), .99) 
      
      dat_survey_w <- trimWeights(dat_survey_w, strict=FALSE,
                                  lower=low_trim, upper=high_trim) 
      summary(weights(dat_survey_w)) }
    
    return(weights(dat_survey_w))
  }
  
  
  # ~~~~~~~~~~~~~ State weights ~~~~~~~~~~~~~~~
  
  if(state & sample_n<50) {
    stop("Bad news: You need a larger sample for proper state weights.") } 
  
  state_file <-  "../DATA/AUX DATA/state_data.csv"
  if(!file.exists(state_file)) state_file <- "http://www.kateto.net/covid19/state_data.csv"
  
  state_data <- tryCatch( 
    {  state_data <- read.csv(state_file, header=T, as.is=T)  },
    error = function(e) {
      stop("Unable to read the state data file from '../DATA/AUX DATA/state_data.csv'") } )
  
  if(party) dat$party_w <- dat$party_vote20
  dat$weight_st <- NA
  
  # Generate weights for each state
  for (selected_state in unique(dat$state_code)) {
    
    pop <- weight_states(state_data, selected_state, dat)
    ind_st <- which(dat$state_code %in% selected_state)
    dat_st <- dat[ind_st,] 
    
    sample_margins <- list("gender", "race_cat", "age_cat", "edu_cat", 
                           "urban_cat", "party_w")
    
    pop_margins <- pop[c("gender_di_st", "race_di_st", "age_di_st", "edu_di_st", 
                         "urban_di_st", "party_di_st")] 
    
    # For some states, there may be 0 people living in certain rural/urban areas:
    urban_st <- sum(pop$urban_di_st$Freq==0)
    
    margin_include <- 1:6
    if(!party) margin_include <- 1:5 
    
    sample_margins <- sample_margins[margin_include]  
    pop_margins    <- pop_margins[margin_include]  
    
    sample_margins <- paste0("~", sample_margins)
    sample_margins <- lapply(sample_margins, as.formula)
    
    dat_state_survey   <- svydesign(ids=~1,  weights=~1, data=dat_st) 
    dat_state_survey_w <- rake(design = dat_state_survey,
                               sample.margins = sample_margins,
                               population.margins = pop_margins)
    weights(dat_state_survey_w)
    summary(weights(dat_state_survey_w))
    
    if(trim) {
      low_trim <-  max(quantile(weights(dat_state_survey_w), .01), .1)
      high_trim <- min(quantile(weights(dat_state_survey_w), .99), 10)
      dat_state_survey_w <- trimWeights(dat_state_survey_w, lower=low_trim, upper=high_trim,  strict=F)
    }
    
    dat$weight_st[ind_st] <- weights(dat_state_survey_w)
  }
  
  summary(dat$weight_st)
  return(dat$weight_st)
}  


# Generate variables used for weighting:
weight_vars <- function(dat, party=FALSE) {
  
  library(dplyr)
  
  dat$party_vote20 <- NA
  
  if ("voted20" %in% colnames(dat)) {
    dat <- mutate(dat, 
                  party_vote20 = case_when(voted20 %in% 1   ~ "biden",
                                           voted20 %in% 2   ~ "trump",
                                           voted20 %in% 3   ~ "other",
                                           voted20 %in% 4:5 ~ "none",
                                           voted20 %in% NA  ~ "none",
                                           voted20 %in% -99 ~  NA_character_ )) 
  }
  table(dat$party_vote20, useNA = "always")
  
  
  # There may be a handful of missing cases in party, vote, or urban type.
  # We can impute them though this becomes a problem if we have more than a few.
  demo_vars <- c("gender", "race", "age", "education",  "state", "urbanicity")
  
  demo_miss  <-  sapply(dat[,demo_vars], anyNA)
  party_miss <- "party"   %in% colnames(dat) &  any(dat$party  %in% c(NA, -99))
  vote_miss  <- "voted20" %in% colnames(dat) &  any(dat$party_vote20  %in% NA) 
  
  imp_vars <- c(demo_vars, "party", "party_vote20") 
  imp_vars <- intersect(imp_vars, colnames(dat)) 
  imp_vars <- dat[,imp_vars]
  imp_vars <- mutate(imp_vars, across(-age, factor))
  imp_vars$age <- as.numeric(imp_vars$age)
  
  if ( any(demo_miss) | (party_miss & (party %in% 1:2)) | (vote_miss & (party %in% 3))) {
    
    message("Imputing missing data -- this may take a second.")
    mice_installed <- suppressPackageStartupMessages(require(mice))
    if(!mice_installed) stop("Please install the 'mice' library for imputation.")
    
    imp <- mice(imp_vars, m=1, maxit=1)
    imp <- complete(imp)
    imp <- mutate(imp, across(-age, as.character))
    imp$education <- as.numeric(imp$education)
    imp$urbanicity <- as.numeric(imp$urbanicity)
    
    table(dat[,demo_vars] == imp[,demo_vars], useNA="always")
    
    if (any(demo_miss))  dat[,demo_vars[demo_miss]] <- imp[,demo_vars[demo_miss]] 
    if (vote_miss)  dat$party_vote20 <- imp$party_vote20 
    if (party_miss) dat$party <- imp$party
    
    detach(package:mice)
  } 
  
  # Recode party variables for weighing
  if ("party" %in% colnames(dat)) {
    
    dat <- mutate(dat, 
                  party_pid = case_when(party %in% "Republican"  ~ "Rep",
                                        party %in% "Democrat"    ~ "Dem",
                                        party %in% "Independent" ~ "Ind",
                                        TRUE ~ "Oth")) 
    table(dat$party_pid, useNA = "always")
    
    dat$party_gender <- dat$party_pid    
    dat <- mutate(dat, 
                  party_gender = case_when(gender %in% "Male"   ~ paste0("M_", party_gender),
                                           gender %in% "Female" ~ paste0("F_", party_gender) ))
    table(dat$party_gender, useNA = "always")
    
  }
  
  # Recode age, race, and education variables for weighing  
  dat <- mutate(dat, 
                age_cat = case_when( age <= 24             ~ "18-24",
                                     age >= 25 & age <= 44 ~ "25-44",
                                     age >= 45 & age <= 64 ~ "45-64",
                                     age >= 65             ~ "65+",
                                     TRUE ~ NA_character_) )
  table(dat$age_cat,useNA = "always")
  
  
  dat$edu_cat <- c("No HS",  "HS", "Some Col", "Col", "Grad")[dat$education] 
  table(dat$edu_cat,useNA = "always")
  
  dat$race_cat <- dat$race
  dat$race_cat[dat$race %in% c('Pacific Islander', 'Native American', "Other race")] <- 'Other'
  table(dat$race_cat, useNA = "always")
  
  dat$urban_cat <- paste0("nchs_",  dat$urbanicity)
  
  # Create interlocking demogrpahics for weighting
  
  dat <- mutate(dat, 
                edu_age = case_when( age < 25 ~ "age_18_24",
                                     age > 24 & dat$age < 45 ~ "age_25_44",   
                                     age > 44 & dat$age < 65 ~ "age_45_64", 
                                     age > 64 ~ "age_65plus" ))
  
  dat <- mutate(dat, 
                edu_age = case_when( education %in% 1 ~ paste0(edu_age, "_no_hs"),
                                     education %in% 2 ~ paste0(edu_age, "_hs"),
                                     education %in% 3 ~ paste0(edu_age, "_some_col"),
                                     education %in% 4 ~ paste0(edu_age, "_ba"), 
                                     education %in% 5 ~ paste0(edu_age, "_grad") ))
  
  table(dat$edu_age, useNA = "always")               
  
  
  
  dat <- mutate(dat, 
                race_sex_age = case_when( age < 25 ~ "18_24",
                                          age > 24 & dat$age < 35 ~ "25_34",
                                          age > 34 & dat$age < 45 ~ "35_44",  
                                          age > 44 & dat$age < 55 ~ "45_54",
                                          age > 54 & dat$age < 65 ~ "55_64",
                                          age > 64 & dat$age < 75 ~ "65_74",
                                          age > 74 ~ "75plus" ))
  dat <- mutate(dat, 
                race_sex_age = case_when(gender %in% "Male"   ~ paste0("male_", race_sex_age),
                                         gender %in% "Female" ~ paste0("female_", race_sex_age)))
  
  dat <- mutate(dat, 
                race_sex_age = case_when(race_cat=="White"            ~ paste0("white_", race_sex_age), 
                                         race_cat=="Hispanic"         ~ paste0("hisp_",  race_sex_age), 
                                         race_cat=="African American" ~ paste0("black_", race_sex_age),
                                         race_cat=="Asian American"   ~ paste0("asian_", race_sex_age),
                                         race_cat=="Other"            ~ paste0("oth_",   race_sex_age)  ))
  table(dat$race_sex_age,useNA = "always")
  
  return(dat)
  
}



# Get the population demographics:
weight_pop <- function() {
  
  out <- list()
  
  # Pew 2021 NPORS Dataset - Party
  out$party_pid_dist <- data.frame(party_w = c("Rep", "Dem", "Ind", "Oth"),
                                   Freq   =  c(.278, .325, .261, .136) )
  
  # Pew 2021 NPORS Dataset - Party by Gender
  out$party_gender_dist <- data.frame(party_w = c("M_Rep", "M_Dem", "M_Ind", "M_Oth",
                                                  "F_Rep", "F_Dem", "F_Ind", "F_Oth"),
                                      Freq =  c(.145, .129, .138, .061, 
                                                .135, .199, .124, .070))
  # Voting 2020 from FEC.gov data
  out$party_vote20_dist <- data.frame(party_w = c("none", "biden", "trump", "other"),
                                      Freq = c(.522, .245, .224, .009) )
  
  # Census/ACS USA demographics
  out$gender_dist <- data.frame(gender = c("Male", "Female"),
                                Freq = c(.484, .516))
  
  out$region_dist <- data.frame(region = c("Northeast", "Midwest", "South", "West"),
                                Freq = c(.18, .21, .38, .24) )
  
  out$race_dist <- data.frame(race_cat = c("White","Hispanic", "African American",
                                           "Asian American", "Other"),
                              Freq = c(.629, .166, .12, .063, .022) )
  
  out$age_dist <- data.frame(age_cat = c("18-24", "25-44", "45-64", "65+"),
                             Freq = c(.12, .34, .33, .20) )
  
  out$edu_dist <- data.frame(edu_cat =  c("No HS", "HS", "Some Col", "Col", "Grad"),
                             Freq =  c(.15, .27, .18, .29, .11) )
  
  out$urban_dist <- data.frame(urban_cat = c("nchs_1", "nchs_2", "nchs_3",
                                             "nchs_4", "nchs_5", "nchs_6"),
                               Freq = c(.305, .247, .209, .092,  .087,  0.061) )
  
  
  edu_age_pop <- c( age_18_24_no_hs    = 0.0169,
                    age_18_24_hs       = 0.03628,
                    age_18_24_some_col = 0.04599,
                    age_18_24_ba       = 0.01362,
                    age_18_24_grad     = 0.001,
                    age_25_44_no_hs    = 0.02618,
                    age_25_44_hs       = 0.08671,
                    age_25_44_some_col = 0.08776,
                    age_25_44_ba       = 0.094,
                    age_25_44_grad     = 0.04979,
                    age_45_64_no_hs    = 0.02905,
                    age_45_64_hs       = 0.09134,
                    age_45_64_some_col = 0.08161,
                    age_45_64_ba       = 0.07336,
                    age_45_64_grad     = 0.04607,
                    age_65plus_no_hs   = 0.02391,
                    age_65plus_hs      = 0.06877,
                    age_65plus_some_col = 0.05532,
                    age_65plus_ba       = 0.04101,
                    age_65plus_grad     = 0.03132)
  
  out$edu_age_dist <- data.frame(edu_age = names(edu_age_pop),
                                 Freq = edu_age_pop)
  
  race_sex_age_pop <- c( white_male_18_24 = 0.0329,
                         white_male_25_34 = 0.05098,
                         white_male_35_44 = 0.04722,
                         white_male_45_54 = 0.05,
                         white_male_55_64 = 0.05738,
                         white_male_65_74 = 0.04456,
                         white_male_75plus = 0.02953,
                         white_female_18_24 = 0.0312,
                         white_female_25_34 = 0.04918,
                         white_female_35_44 = 0.04661,
                         white_female_45_54 = 0.05017,
                         white_female_55_64 = 0.06007,
                         white_female_65_74 = 0.04929,
                         white_female_75plus = 0.0405,
                         black_male_18_24 = 0.0086,
                         black_male_25_34 = 0.01292,
                         black_male_35_44 = 0.01006,
                         black_male_45_54 = 0.00951,
                         black_male_55_64 = 0.00903,
                         black_male_65_74 = 0.00529,
                         black_male_75plus = 0.00273,
                         black_female_18_24 = 0.00838,
                         black_female_25_34 = 0.01304,
                         black_female_35_44 = 0.01111,
                         black_female_45_54 = 0.01081,
                         black_female_55_64 = 0.01069,
                         black_female_65_74 = 0.00705,
                         black_female_75plus = 0.00475,
                         asian_male_18_24 = 0.00344,
                         asian_male_25_34 = 0.00632,
                         asian_male_35_44 = 0.00569,
                         asian_male_45_54 = 0.00485,
                         asian_male_55_64 = 0.00384,
                         asian_male_65_74 = 0.00261,
                         asian_male_75plus = 0.00167,
                         asian_female_18_24 = 0.00339,
                         asian_female_25_34 = 0.00657,
                         asian_female_35_44 = 0.00644,
                         asian_female_45_54 = 0.00556,
                         asian_female_55_64 = 0.00461,
                         asian_female_65_74 = 0.00333,
                         asian_female_75plus = 0.00232,
                         hisp_male_18_24 = 0.01399,
                         hisp_male_25_34 = 0.01982,
                         hisp_male_35_44 = 0.01776,
                         hisp_male_45_54 = 0.0143,
                         hisp_male_55_64 = 0.00987,
                         hisp_male_65_74 = 0.00515,
                         hisp_male_75plus = 0.00291,
                         hisp_female_18_24 = 0.01332,
                         hisp_female_25_34 = 0.01811,
                         hisp_female_35_44 = 0.01671,
                         hisp_female_45_54 = 0.01408,
                         hisp_female_55_64 = 0.01035,
                         hisp_female_65_74 = 0.00609,
                         hisp_female_75plus = 0.00432,
                         oth_male_18_24 = 0.00064,
                         oth_male_25_34 = 0.00097,
                         oth_male_35_44 = 0.00078,
                         oth_male_45_54 = 0.0007,
                         oth_male_55_64 = 0.00068,
                         oth_male_65_74 = 0.00043,
                         oth_male_75plus = 0.00022,
                         oth_female_18_24 = 0.00061,
                         oth_female_25_34 = 0.00094,
                         oth_female_35_44 = 0.00078,
                         oth_female_45_54 = 0.00073,
                         oth_female_55_64 = 0.00076,
                         oth_female_65_74 = 0.00049,
                         oth_female_75plus  = 0.0003 )
  
  
  out$race_sex_age_dist <- data.frame(race_sex_age = names(race_sex_age_pop),
                                      Freq = race_sex_age_pop)
  
  return(out)
}


weight_states <- function(state_data=NULL, selected_state="CA", dat ) {
  
  state_data$oth <- state_data$other + state_data$natam
  
  ind_st <- which(dat$state_code==selected_state)
  dat_st <- dat[ind_st,]
  dist_st   <- state_data[state_data$code==selected_state,]
  
  gender_st <- dist_st[c("male", "female")]
  race_st   <- dist_st[c("white", "latino", "black", "asian", "oth")] 
  age_st    <- dist_st[c("age_18_24", "age_25_44", "age_45_64", "age_65" )]
  edu_st    <- dist_st[c("edu_some_hs", "edu_hs", "edu_some_col", "edu_ba", "edu_grad" )]
  urban_st  <- dist_st[c("nchs_1", "nchs_2", "nchs_3",
                         "nchs_4", "nchs_5", "nchs_6")]
  
  race_inc  <-  c( 'White', 'Hispanic', 'African American', 
                   'Asian American', 'Other') %in% dat_st$race_cat
  age_inc   <-  c("18-24", "25-44", "45-64", "65+") %in% dat_st$age_cat
  edu_inc   <-  1:5 %in% dat_st$education
  urban_inc <-  paste0("nchs_",1:6) %in% dat_st$urban_cat
  
  party_st <-  dist_st[c("perc_biden_20", "perc_trump_20", "perc_other_20", "perc_vote_20")]
  names(party_st) <- c("biden", "trump", "other", "none")
  party_st["none"] <- 1-party_st["none"]
  party_inc <-   names(party_st) %in% dat_st$party_vote20
  
  excluded <- sapply(list(race_inc, age_inc, edu_inc, party_inc), function(x) any(!x))
  if(!"voted20" %in% colnames(dat)) excluded[4] <- FALSE
  if (any(excluded)) message("NOTE: Missing strata ", 
                             paste(c("race_inc", 'age_inc', "edu_inc", "party_inc")[excluded], collapse=" "),
                             " in the great state of ", selected_state)
  out <- list()
  
  sample_st <- nrow(dat_st)
  
  out$gender_di_st <- data.frame(gender = c("Male", "Female"),
                                 Freq = sample_st * unlist(gender_st) )
  
  out$race_di_st <- data.frame(race_cat = c("White","Hispanic", "African American",
                                            "Asian American", "Other")[race_inc],
                               Freq = sample_st * unlist(race_st)[race_inc] )
  
  out$age_di_st <- data.frame(age_cat = c("18-24", "25-44", "45-64", "65+")[age_inc], 
                              Freq = sample_st * unlist(age_st)[age_inc] )
  
  out$edu_di_st <- data.frame(edu_cat =  c("No HS", "HS", "Some Col", "Col", "Grad")[edu_inc],
                              Freq = sample_st * unlist(edu_st)[edu_inc])
  
  out$urban_di_st <- data.frame(urban_cat = c("nchs_1", "nchs_2", "nchs_3",
                                              "nchs_4", "nchs_5", "nchs_6")[urban_inc],
                                Freq = sample_st * unlist(urban_st)[urban_inc] ) 
  
  if("voted20" %in% colnames(dat)) {
    out$party_di_st <- data.frame(party_w = c("biden", "trump", "other", "none")[party_inc],
                                  Freq = sample_st * unlist(party_st)[party_inc] ) }
  
  # State satistics are approximate, so 0 strata means very few cases rounded down
  out <- lapply(out, function(x) {x$Freq[x$Freq==0] <- 1; x})
  
  return(out)
  
}
