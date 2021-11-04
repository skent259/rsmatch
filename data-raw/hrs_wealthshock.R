#' Purpose: code to prepare `hrs_wealthshock` dataset
#'
#' Due to restrictions on sharing the data, it is not included in this package
#' distribution.  However, the processing steps for the data are provided here
#' to generate an identical analysis data set from which matching can be done.
#'
#' Instructions for data download:
#'
#' 1. Download the Version 1992-2016 v1: released May 2019 RAND HRS Longitudinal
#' Data in STATA format from
#' https://hrsdata.isr.umich.edu/data-products/rand-hrs-archived-data-products.
#' This requires making an account to access the data.
#' 2. Unzip the .zip file, and move  randhrs1992_2016v1.dta to an appropriate
#' directory that matches the `fname` object below.
#' 3. Depending on where you run this file, you may need to move
#' ./data-raw/bls_cpi.xlsx to the appropriate directory and change `fname`.
#' This gives the consumer price index for each year and month between 01-1992
#' and 10-2019.
#'

library(tidyverse)
library(lubridate)
library(here) # use common sense file paths
library(haven) # `read_dta()`
library(readxl) # For CPI data from BLS
library(rlang)

## pull in raw data
fname <- here("data-raw/randhrs1992_2016v1.dta")
d <- read_dta(file = fname)

## Getting variables needed for analysis
d_analysis <-
  d %>%
  filter(hacohort == 3) %>%
  select(
    hhidpn, # unique identifier
    matches("r[0-9]*wthh$"), # household weights
    radyear, radmonth, raddate, radsrc, # death information
    # Baseline covariates
    r1agem_b, # age at enrollment
    ragender, # self-reported sex
    raracem, # self-reported race
    raedyrs, # educational Attainment
    matches("h[0-9]*atotb$"), # household net worth (including 2nd residence)
    matches("h[0-9]*atota$"), # household net worth (not including 2nd residence)
    matches("r[0-9]*smokev$"), # smoking status
    matches("r[0-9]*drink$"), # alcohol consumption
    matches("r[0-9]*vigact$"), # physical activity
    matches("r[0-9]*bmi$"), # BMI
    rahispan, # self-reported ethnicity
    # Time varying covariates
    matches("h[0-9]*itot$"), # household total income
    matches("r[0-9]*mstat$"), # marital status
    matches("r[0-9]*lbrf$"), # labor-force status
    #   health insurance status
    matches("r[0-9]*higov$"), # covered by federal govnt health insurance
    matches("r[0-9]*covr$"), # covered by employer sponsored health insurance
    matches("r[0-9]*hiothp$"), # covered by other health insurance
    #
    matches("r[0-9]*shlt$"), # self-rated health
    matches("r[0-9]*hlthlm$"), # health limitation on work ability
    matches("r[0-9]*hosp$"), # hospitalization
    matches("r[0-9]*oopmd$"), # out-of-pocket HC costs
    #   history of 8 conditions
    matches("r[0-9]*hibpe$"), # high-blood pressure or hypertension
    matches("r[0-9]*diabe$"), # diabetes or high blood sugar
    matches("r[0-9]*cancre$"), # cancer
    matches("r[0-9]*lunge$"), # chronic lung disease
    matches("r[0-9]*hearte$"), # heart attack
    matches("r[0-9]*stroke$"), # stroke
    matches("r[0-9]*psyche$"), # emotional, nervous or psychiatric problems
    matches("r[0-9]*arthre$"), # arthritis or rheumatism
    matches("r[0-9]*conde$"), # count of chronic conditions
    #   limitations of 5 daily activities (ADLs)
    matches("r[0-9]*walkr$"), # walking across a room
    matches("r[0-9]*dress$"), # dressing
    matches("r[0-9]*bath$"), # bathing
    matches("r[0-9]*eat$"), # eating
    matches("r[0-9]*bed$"), # getting in and out of bed
    #
    matches("r[0-9]*risk$"), # financial risk aversion
    r1beqlrg, # expectation of leaving a large bequest at death
    hhid, hacohort # other identifiers and variables
  ) %>%
  mutate_at( # >= 2 chronic conditions
    .vars = vars(matches("r[0-9]*conde$")),
    .funs = list(ge2 = ~1*(. >= 2))
  ) %>%
  select(-matches("r[0-9]*conde$"))

#' Create health insured variable `r_w_hiany` as being covered by govnt, employer
#' sponsored, or other health insurance for each wave

for (w in 1:13) {
  rwhiany <- paste0("r", w, "hiany")
  rwhigov <- paste0("r", w, "higov")
  rwcovr <- paste0("r", w, "covr")
  rwhiothp <- paste0("r", w, "hiothp")

  expr <- paste(c(rwhigov,rwcovr,rwhiothp), collapse = "|") # e.g. r1higov | r1covr | r1hiothp
  d_analysis <-
    d_analysis %>%
    mutate(!!rwhiany := as.factor(1*(.data[[rwhigov]] + .data[[rwcovr]] + .data[[rwhiothp]] > 0)))
}



#' Create race/ethnicity variable to match Pool et al. paper
#' levels are non-Hispanic black, non-Hispanic white, Hispanic, or other race
#' mapped to levels 1,2,3,4 respectively

d_analysis <-
  d_analysis %>%
  mutate(raethrace = as.factor(case_when(
    rahispan == 0 & raracem == 2 ~ 1, # non-Hispanic black
    rahispan == 0 & raracem == 1 ~ 2, # non-Hispanic white
    rahispan == 1                ~ 3, # Hispanic
    TRUE                         ~ 4 # other race
  )))

#' Out of pocket medical expenses are all NA for `r2oopmd`, and don't exists
#' for `r1oopmd`.  Best thing to do is set these to be `r3oopmd` for balance
#' and matching purposes

d_analysis <-
  d_analysis %>%
  mutate(r1oopmd = r3oopmd,
         r2oopmd = r3oopmd)

#' Adjust household wealth and total income variables
#' to be in 2014 dollars

# CPI data from https://data.bls.gov/pdq/SurveyOutputServlet
fname <- here("data-raw/bls_cpi.xlsx")
cpi <- read_excel(fname, skip = 11) %>%
  select(Year, Annual)

# Make 100 = 2014
index2014 <- cpi %>% filter(Year == 2014) %>% pull(Annual)
cpi <-
  cpi %>%
  mutate(cpi2014 = Annual / index2014)

cols_to_add <-
  d_analysis %>%
  select(hhidpn, matches("h[0-9]*atota$"), matches("h[0-9]*itot$")) %>%
  pivot_longer(-hhidpn,
               values_to = "value",
               names_to = "wave_var",
  ) %>%
  mutate(wave = str_extract(wave_var, "[0-9]+")) %>%
  mutate(year = 1990 + as.numeric(wave)*2) %>%
  left_join(cpi,
            by = c("year" = "Year")) %>%
  mutate(value_adj = value / cpi2014) %>%
  mutate(wave_var = str_c(wave_var,"_cpi14")) %>%
  select(hhidpn, wave_var, value_adj) %>%
  pivot_wider(hhidpn,
              values_from = value_adj,
              names_from = wave_var
  )

d_analysis <-
  d_analysis %>%
  left_join(cols_to_add, by = "hhidpn")


#' Calculate wealth shock as change of wealth
#' NOTE: hwwlthshck = 2 when current wealth is negative, in which case
#' they should have a past wealth shock (or be classified in asset
#' poverty at baseline).  Value of 2 is primarily to make sure this
#' variable doesn't get used in the wrong way
d_analysis <-
  d_analysis %>%
  mutate(
    h2wlthshck = ifelse(h1atota_cpi14 < 0, 2, 1*(h2atota_cpi14 / h1atota_cpi14 < 0.25)),
    h3wlthshck = ifelse(h2atota_cpi14 < 0, 2, 1*(h3atota_cpi14 / h2atota_cpi14 < 0.25)),
    h4wlthshck = ifelse(h3atota_cpi14 < 0, 2, 1*(h4atota_cpi14 / h3atota_cpi14 < 0.25)),
    h5wlthshck = ifelse(h4atota_cpi14 < 0, 2, 1*(h5atota_cpi14 / h4atota_cpi14 < 0.25)),
    h6wlthshck = ifelse(h5atota_cpi14 < 0, 2, 1*(h6atota_cpi14 / h5atota_cpi14 < 0.25)),
    h7wlthshck = ifelse(h6atota_cpi14 < 0, 2, 1*(h7atota_cpi14 / h6atota_cpi14 < 0.25)),
    h8wlthshck = ifelse(h7atota_cpi14 < 0, 2, 1*(h8atota_cpi14 / h7atota_cpi14 < 0.25)),
    h9wlthshck = ifelse(h8atota_cpi14 < 0, 2, 1*(h9atota_cpi14 / h8atota_cpi14 < 0.25)),
    h10wlthshck = ifelse(h9atota_cpi14 < 0, 2, 1*(h10atota_cpi14 / h9atota_cpi14 < 0.25)),
    h11wlthshck = ifelse(h10atota_cpi14 < 0, 2, 1*(h11atota_cpi14 / h10atota_cpi14 < 0.25)),
    h12wlthshck = ifelse(h11atota_cpi14 < 0, 2, 1*(h12atota_cpi14 / h11atota_cpi14 < 0.25)),
    h13wlthshck = ifelse(h12atota_cpi14 < 0, 2, 1*(h13atota_cpi14 / h12atota_cpi14 < 0.25)),
  )


#' Add column for date of death in lubridate format
#' "raddate" is time from 1960, january 1
d_analysis <- d_analysis %>%
  mutate(raddate_lub = lubridate::ymd("19600101") + lubridate::ddays(raddate))

#' Calculate the first wealthshock time and whether an individual experiences
#' a wealthshock, has asset poverty at baseline, or never experiences a
#' wealth shock

calc_wlthshck <- function(waveshock){
  if (all(is.na(1 == waveshock))) NA
  else any(1 == waveshock, na.rm = TRUE)
}

first_wealthshock <-
  select(d_analysis, hhidpn, contains("wlthshck")) %>%
  pivot_longer(-hhidpn, values_to = "value", names_to = "wave_var") %>%
  mutate(wave = as.numeric(str_extract(wave_var, "[0-9]+"))) %>%
  mutate(wlthshck_wave = ifelse(value == 1,wave, NA)) %>%
  group_by(hhidpn) %>%
  summarize(first_wlthshck = min(wlthshck_wave, na.rm = TRUE))

# Create treatment variable and time of treatment
d_analysis <-
  d_analysis %>%
  mutate(
    treatment = case_when(
      h1atota_cpi14 < 0 ~ -1, # Asset poverty as baseline
      apply(select(., contains("wlthshck")), 1, calc_wlthshck) ~ 1, # wealth shock
      is.na(apply(select(., contains("wlthshck")), 1, calc_wlthshck)) ~ -1000, # NA
      TRUE ~ 0
    )
  ) %>%
  left_join(first_wealthshock, by = "hhidpn") %>%
  mutate(treatment_time = ifelse(treatment == 1, first_wlthshck, NA))


# Change column types
d_analysis <-
  d_analysis %>%
  mutate(raedyrs = as.numeric(raedyrs)) %>% # remove haven labeled class (this is really not a factor)
  mutate_at(vars(contains("conde_ge2")), as_factor) %>%
  mutate_if(haven::is.labelled, as_factor) # take all labeled classes and make them into factors

# Keep only wealth shock or positive wealth without shock
d_analysis <-
  d_analysis %>%
  filter(treatment != -1000) %>%
  filter(treatment != -1)

#' Take in baseline variables `bl_var` and timevarying variables `t_var` and
#' a dataframe with these variables, along with hhidpn as unique identifier
#'
#' @Return a dataset (unique at wave x hhidpn pair) which has bl_var, t_var for
#' each hhidpn and wave
#' Note: assumes 13 waves of data
generate_wave_based_data <- function(df, bl_var, t_var) {

  out <- list()
  for (wave in 1:13) {
    out_i <-
      df[, c("hhidpn", bl_var,str_replace(t_var, "_w_", as.character(wave)))] %>%
      mutate(wave = wave)

    colnames(out_i) <- c("hhidpn", bl_var,t_var, "wave")
    out[[wave]] <- out_i
  }
  do.call(rbind, out)
}

baseline_vars <- c(
  "r1agem_b", "ragender", "raethrace",
  "raedyrs", "h1atotb", "r1smokev",
  "r1drink", "r1vigact", "r1bmi"
)

timevarying_vars <- c(
  "h_w_atota_cpi14", "r_w_mstat", "r_w_lbrf",
  "r_w_hiany", "r_w_shlt", "r_w_hlthlm", "r_w_oopmd",
  "r_w_hosp", "r_w_conde_ge2",
  "r_w_hibpe","r_w_diabe","r_w_cancre","r_w_lunge",
  "r_w_hearte","r_w_stroke","r_w_psyche","r_w_arthre",
  "r_w_walkr","r_w_dress",
  "r_w_bath","r_w_eat","r_w_bed"
)


hrs_wealthshock <-
  d_analysis %>%
  generate_wave_based_data(baseline_vars, timevarying_vars) %>%
  na.omit() %>%
  left_join(select(d_analysis, hhidpn, treatment, treatment_time),
            by = "hhidpn") %>%
  select(hhidpn, wave, treatment_time, treatment, everything())

hrs_wealthshock <- hrs_wealthshock %>%
  mutate(across(c("r_w_walkr","r_w_dress", "r_w_bath","r_w_eat","r_w_bed"),
                ~ !(.x %in% c("1.not at all diff", "0.no", "9.don't do")),
                .names = "{col}_flag")) %>%
  mutate(r_w_adllim = as.factor(1*(r_w_walkr_flag + r_w_dress_flag + r_w_bath_flag + r_w_eat_flag + r_w_bed_flag >= 1))) %>%
  select(-c("r_w_walkr","r_w_dress", "r_w_bath","r_w_eat","r_w_bed"), -ends_with("flag"))


save(hrs_wealthshock, file = here("data/hrs_wealthshock.rda"))
# usethis::use_data(hrs_wealthshock, overwrite = TRUE)
