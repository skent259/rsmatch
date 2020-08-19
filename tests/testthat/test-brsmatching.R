context("Test the functions used in Balanced Risk-set Matching brsmatching.R")

test_that("enumerate_edges_and_compute_distances has correct output", {
  library(dplyr)
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


test_that("rsm_optimization_model has correct output", {
  df <- data.frame(
    hhidpn = rep(1:3, each = 3),
    wave = rep(1:3, 3),
    treatment_time = rep(c(2,3,NA), each = 3),
    X1 = c(2,2,2,3,3,3,9,9,9),
    X2 = rep(c("a","a","b"), each = 3),
    X3 = c(9,4,5,6,7,2,3,4,8),
    X4 = c(8,9,4,5,6,7,2,3,4)
  )
  edges <- enumerate_edges_and_compute_distances(df, "hhidpn", "wave", "treatment_time")
  bal <- balance_columns(df, "hhidpn", "wave", "treatment_time")
  n_unique_id <- length(unique(df$hhidpn))

  ## Gurobi, balanced
  model <- rsm_optimization_model(1, edges, bal, optimizer = "gurobi", balance = TRUE)
  expect_equal(names(model), c("modelsense", "obj", "varnames", "A", "sense", "rhs", "vtype"))
  expect_equal(model$obj[1:nrow(edges)], edges$dist)
  expect_equivalent(as.matrix(model$A[3:(2+n_unique_id), 1:nrow(edges)]),
                    matrix(c(1,0,1, 0,1,1, 1,1,0), nrow = n_unique_id, byrow = TRUE))
  ## Gurobi, unbalanced
  model <- rsm_optimization_model(1, edges, bal, optimizer = "gurobi", balance = FALSE)
  expect_equal(nrow(model$A), 2 + n_unique_id)
  expect_equal(unique(model$vtype), "B")
  ## GLPK, balanced
  model <- rsm_optimization_model(1, edges, bal, optimizer = "glpk", balance = TRUE)
  expect_equal(names(model), c("max", "obj", "varnames", "mat", "dir", "rhs", "types"))
})


test_that("brsmatch has correct output", {
  df <- data.frame(
    hhidpn = rep(1:3, each = 3),
    wave = rep(1:3, 3),
    treatment_time = rep(c(2,3,NA), each = 3),
    X1 = c(2,2,2,3,3,3,9,9,9),
    X2 = rep(c("a","a","b"), each = 3),
    X3 = c(9,4,5,6,7,2,3,4,8),
    X4 = c(8,9,4,5,6,7,2,3,4)
  )

  brsmatch(n_pairs = 1, df = df, id = "hhidpn", time = "wave", trt_time = "treatment_time",
           optimizer = "glpk")


  brsmatch(n_pairs = 1, df = df, id = "hhidpn", time = "wave", trt_time = "treatment_time",
           optimizer = "glpk", balance = FALSE)

  # TODO: come up with a few tests here, check gurobi works
})

