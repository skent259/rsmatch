context("Test the functions used in Balanced Risk-set Matching brsmatching.R")
library(dplyr)

test_that("enumerate_edges_and_compute_distances has correct output", {
  wave_data <- data.frame(
    hhidpn = rep(1:3, each = 3),
    wave = rep(1:3, 3),
    X1 = c(2,2,2,3,3,3,9,9,9),
    X2 = rep(c("a","a","b"), each = 3),
    X3 = c(9,4,5,6,7,2,3,4,8),
    X4 = c(8,9,4,5,6,7,2,3,4)
  )
  treated_data <- data.frame(
    hhidpn = 1:2,
    treatment = c(1,1),
    treatment_time = c(2,3)
  )

  df <- wave_data %>%
    left_join(treated_data, by = "hhidpn") %>%
    select(-treatment)
  # id <- "hhidpn"
  # time <- "wave"
  # trt_time <- "treatment_time"
  edges <- enumerate_edges_and_compute_distances(df, "hhidpn", "wave", "treatment_time")

  expect_equal(edges$trt_id, c(1,1,2))
  expect_equal(edges$all_id, c(2,3,3))
  expect_equal(edges$trt_time, c(2,2,3))
  expect_equal(round(edges$dist,3),
               c(7.647, 11.309, 7.068))


})
