## code to prepare `oasis` dataset goes here

# downloaded from https://www.kaggle.com/jboysen/mri-and-alzheimers?select=oasis_longitudinal.csv
df <- read.csv("data-raw/oasis_longitudinal.csv")

# select participants to include
df1 <- df %>%
  dplyr::group_by(Subject.ID) %>%
  dplyr::filter(
    identical(unique(CDR), c(0.5, 1)) | # MCI to AD
      identical(unique(CDR), 0.5) # MCI throughout
    # excludes AD in baseline, CN throughout
  ) %>%
  dplyr::ungroup()

# calculate treatment time, treatment, and select appropriate variables
oasis <- df1 %>%
  dplyr::mutate(trt_time = ifelse(CDR == 1, Visit, NA)) %>%
  dplyr::group_by(Subject.ID) %>%
  dplyr::mutate(
    trt_time = min(trt_time, na.rm = TRUE),
    trt_time = ifelse(trt_time == Inf, NA, trt_time),
    trt = ifelse(is.na(trt_time), 0, 1),
    SES = ifelse(is.na(SES), -1, SES),
    SES = as.factor(SES)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    Subject.ID, Visit, trt_time, # core variables
    M.F, EDUC, SES, # baseline covariates
    Age, MR.Delay, eTIV, nWBV, ASF
  )

# clean up the names
colnames(oasis) <- c("subject_id", "visit", "time_of_ad",
                   "m_f", "educ", "ses",
                   "age", "mr_delay", "e_tiv", "n_wbv", "asf")


usethis::use_data(oasis, overwrite = TRUE)
