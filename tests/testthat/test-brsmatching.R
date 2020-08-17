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
  edges <- enumerate_edges_and_compute_distances(df, "hhidpn", "wave", "treatment_time")

  expect_equal(edges$trt_id, c(1,1,2))
  expect_equal(edges$all_id, c(2,3,3))
  expect_equal(edges$trt_time, c(2,2,3))
  expect_equal(round(edges$dist,3), c(7.647, 11.309, 7.068))
})



test_that("balance_columns has correct output", {
  df <- data.frame(
    hhidpn = rep(1:3, each = 3),
    wave = rep(1:3, 3),
    treatment_time = rep(c(2,3,NA), each = 3),
    X1 = c(2,2,2,3,3,3,9,9,9),
    X2 = rep(c("a","a","b"), each = 3),
    X3 = c(9,4,5,6,7,2,3,4,8),
    X4 = c(8,9,4,5,6,7,2,3,4)
  )

  balance_covariates <- c("X1", "X2", "X3", "X4")
  bal <- balance_columns(df, "hhidpn", "wave", "treatment_time", balance_covariates)

  expect_equal(length(colnames(bal)), 10)
  expect_equivalent(bal[, "X1.q1"], c(rep(0,3), rep(1,6)))
  expect_equivalent(bal[, "X1.q2"], c(rep(0,3), rep(1,6)))
  expect_equivalent(bal[, "X2a"], c(rep(1,6), rep(0,3)))
  expect_equivalent(bal[, "X2b"], c(rep(0,6), rep(1,3)))
  expect_equivalent(bal[, "X3.q1"], c(1,1,1,1,1,0,1,1,1))
  expect_equivalent(bal[, "X3.q2"], c(1,1,1,1,1,0,0,1,1))
  expect_equivalent(bal[, "X4.q1"], c(1,1,0,0,0,0,0,0,0))
  expect_equivalent(bal[, "X4.q2"], c(0,1,0,0,0,0,0,0,0))

})
