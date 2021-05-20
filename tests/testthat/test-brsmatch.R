context("Test the functions used in Balanced Risk-set Matching brsmatching.R")

check_for_gurobi <- function() {
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    skip("The package gurobi is not available.")
  }
}

check_for_glpk <- function() {
  if (!requireNamespace("Rglpk", quietly = TRUE)) {
    skip("The package Rglpk is not available.")
  }
}

test_that("compute_distances has correct output", {
  suppressWarnings(library(dplyr))
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
  edges <- .compute_distances(df, "hhidpn", "wave", "treatment_time")

  df$treatment_time <- df$treatment_time - 1
  edges <- .compute_distances(df, "hhidpn", "wave", "treatment_time", options = list(time_lag = TRUE))

  expect_equal(edges$trt_id, c(1,1,2))
  expect_equal(edges$all_id, c(2,3,3))
  expect_equal(edges$trt_time, c(2,2,3)-1)
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

  # mix of character and numeric variables
  balance_covariates <- c("X1", "X2", "X3", "X4")
  df$treatment_time <- df$treatment_time - 1
  bal <- .balance_columns(df, "hhidpn", "wave", "treatment_time", balance_covariates)

  expect_equal(length(colnames(bal)), 10)
  expect_equivalent(bal[, "X1.q1"], c(rep(0,3), rep(1,6)))
  expect_equivalent(bal[, "X1.q2"], c(rep(0,3), rep(1,6)))
  expect_equivalent(bal[, "X2a"], c(rep(1,6), rep(0,3)))
  expect_equivalent(bal[, "X2b"], c(rep(0,6), rep(1,3)))
  expect_equivalent(bal[, "X3.q1"], c(1,0,0,0,0,0,0,0,1))
  expect_equivalent(bal[, "X3.q2"], c(1,0,0,0,0,0,0,0,0))
  expect_equivalent(bal[, "X4.q1"], c(1,1,0,0,0,1,0,0,0))
  expect_equivalent(bal[, "X4.q2"], c(1,1,0,0,0,0,0,0,0))

  # single character balance variable
  balance_covariates <- c("X2")
  bal <- .balance_columns(df, "hhidpn", "wave", "treatment_time", balance_covariates)

  expect_equal(length(colnames(bal)), 4)
  expect_equivalent(bal[, "X2a"], c(rep(1,6), rep(0,3)))
  expect_equivalent(bal[, "X2b"], c(rep(0,6), rep(1,3)))

  # single numeric balance variable
  balance_covariates <- c("X1")
  bal <- .balance_columns(df, "hhidpn", "wave", "treatment_time", balance_covariates)

  expect_equal(length(colnames(bal)), 4)
  expect_equivalent(bal[, "X1.q1"], c(rep(0,3), rep(1,6)))
  expect_equivalent(bal[, "X1.q2"], c(rep(0,3), rep(1,6)))

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
  df$treatment_time <- df$treatment_time - 1
  edges <- .compute_distances(df, "hhidpn", "wave", "treatment_time")
  bal <- .balance_columns(df, "hhidpn", "wave", "treatment_time")
  n_unique_id <- length(unique(df$hhidpn))

  ## GLPK, balanced
  model <- .rsm_optimization_model(1, edges, bal, optimizer = "glpk", balance = TRUE)
  expect_equal(names(model), c("max", "obj", "varnames", "mat", "dir", "rhs", "types"))

  ## Gurobi, balanced
  model <- .rsm_optimization_model(1, edges, bal, optimizer = "gurobi", balance = TRUE)
  expect_equal(names(model), c("modelsense", "obj", "varnames", "A", "sense", "rhs", "vtype"))
  expect_equal(model$obj[1:nrow(edges)], edges$dist)
  expect_equivalent(as.matrix(model$A[3:(2+n_unique_id), 1:nrow(edges)]),
                    matrix(c(1,0,1, 0,1,1, 1,1,0), nrow = n_unique_id, byrow = TRUE))
  ## Gurobi, unbalanced
  model <- .rsm_optimization_model(1, edges, bal, optimizer = "gurobi", balance = FALSE)
  expect_equal(nrow(model$A), 2 + n_unique_id)
  expect_equal(unique(model$vtype), "B")
})

test_that("rsm_optimization_model() doesn't re-order the edges", {
  edges <- data.frame(
    trt_id = c(1, 1, 1, 2, 2, 2, 3, 4, 3),
    all_id = c(2, 3, 4, 3, 4, 5, 2, 3, 5),
    trt_time = c(2, 2, 2, 3, 3, 3, 2, 2, 2),
    dist = 0.1 * 3:11
  )

  bal <- data.frame(
    id = rep(1:5, each = 3),
    time = rep(1:3, 5),
    X2a = c(rep(1, 10), rep(0, 5)),
    X1.q1 = c(1, rep(0, 13), 1)
  )
  verbose <- interactive()

  model <- .rsm_optimization_model(2, edges, bal, optimizer = "gurobi", balance = TRUE)
  expect_equal(edges$dist, model$obj[1:nrow(edges)])

  model <- .rsm_optimization_model(2, edges, bal, optimizer = "glpk", balance = TRUE)
  check_for_glpk()
  res <- with(model, Rglpk::Rglpk_solve_LP(obj, mat, dir, rhs, types = types, max = max,
                                           control = list(verbose = verbose, presolve = TRUE)))
  matches <- res$solution[grepl("f", model$varnames)]
  matched_ids <- edges[matches == 1, c("trt_id", "all_id")]
  .output_pairs(matched_ids, "id", id_list = sample(unique(bal$id)))

  model <- .rsm_optimization_model(2, edges, bal, optimizer = "gurobi", balance = FALSE)
  expect_equal(edges$dist, model$obj[1:nrow(edges)])

})

test_that("output_pairs() has correct output.", {
  matched_ids <- data.frame(
    trt_id = c(2, 5, 7),
    all_id = c(3, 1, 4)
  )

  out <- .output_pairs(matched_ids)
  expect_equal(colnames(out), c("id", "pair_id", "type"))
  expect_setequal(out$id, c(1,2,3,4,5,7))
  expect_equal(unique(as.vector(table(out$pair_id))), 2)
  expect_equal(unique(as.vector(table(out$type))), nrow(matched_ids))

  out <- .output_pairs(matched_ids, id_list = 1:10)
  expect_setequal(out$id, 1:10)
  expect_equal(unique(as.vector(table(out$pair_id))), 2)
  expect_equal(unique(as.vector(table(out$type))), nrow(matched_ids))
  expect(is.na(unique(out[out$id %in% c(6, 8, 9, 10), "pair_id"])), "pair_id should have NA values at id 6, 8, 9, and 10")
  expect(is.na(unique(out[out$id %in% c(6, 8, 9, 10), "type"])), "type should have NA values at id 6, 8, 9, and 10")

  out <- .output_pairs(matched_ids, id = "hhidpn")
  expect_equal(colnames(out), c("hhidpn", "pair_id", "type"))

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

  check_for_glpk()
  pairs <- brsmatch(n_pairs = 1, data = df, id = "hhidpn", time = "wave", trt_time = "treatment_time",
                    options = list(time_lag = TRUE, optimizer = "glpk"))
  expect_equal(colnames(pairs), c("hhidpn", "pair_id", "type"))
  expect_equal(length(unique(na.omit(pairs$pair_id))), 1)
  expect_equal(unique(pairs$hhidpn), unique(df$hhidpn))
  expect_equal(pairs$hhidpn[which(pairs$pair_id == 1)], c(2,3))

  # check runs properly with other arguments
  brsmatch(n_pairs = 1, data = df, id = "hhidpn", time = "wave", trt_time = "treatment_time",
           balance = FALSE)

  check_for_gurobi()
  pairs <- brsmatch(n_pairs = 1, data = df, id = "hhidpn", time = "wave", trt_time = "treatment_time",
                    options = list(time_lag = TRUE, optimizer = "gurobi"))
})


test_that("options 'between period treatment' works with dead individuals", {
  df <- data.frame(
    hhidpn = rep(1:4, each = 3),
    wave = rep(1:3, 4),
    treatment_time = rep(c(2,3,NA,NA), each = 3),
    X1 = c(2,2,2,3,3,3,9,9,9,2,2,2),
    X2 = rep(c("a","a","b","a"), each = 3),
    X3 = c(9,4,5,6,7,2,3,4,8,9,3,5),
    X4 = c(8,9,4,5,6,7,2,3,4,9,8,7)
  )
  df <- df[-12, ] # hhidpn=4 dies at wave=3

  check_for_glpk()
  pairs <- brsmatch(n_pairs = 1, data = df, id = "hhidpn", time = "wave", trt_time = "treatment_time",
                    options = list(time_lag = TRUE, optimizer = "glpk"))
  expect_equal(pairs %>% filter(!is.na(pair_id)) %>% pull(hhidpn),
               c(1,4))

  df$treatment_time <- df$treatment_time - 1
  edges <- .compute_distances(df, "hhidpn", "wave", "treatment_time", options = list(time_lag = TRUE))
  # expect that the possible pairs should be 1,2 1,3 1,4 or
  expect_equal(edges$trt_id, c(1,1,1,2))
  expect_equal(edges$all_id, c(2,3,4,3))
  expect_equal(round(edges$dist, 4), c(6.6990, 9.8932, 0.4854, 8.2104))

})


test_that("`brsmatch()` works when 'id' is a character vector", {
  df <- data.frame(
    hhidpn = rep(1:3, each = 3),
    wave = rep(1:3, 3),
    treatment_time = rep(c(2,3,NA), each = 3),
    X1 = c(2,2,2,3,3,3,9,9,9),
    X2 = rep(c("a","a","b"), each = 3),
    X3 = c(9,4,5,6,7,2,3,4,8),
    X4 = c(8,9,4,5,6,7,2,3,4)
  )

  check_for_glpk()
  pairs1 <- brsmatch(n_pairs = 1, data = df, id = "hhidpn", time = "wave", trt_time = "treatment_time",
                     options = list(time_lag = TRUE, optimizer = "glpk"))

  df$hhidpn <- as.character(df$hhidpn)

  pairs2 <- brsmatch(n_pairs = 1, data = df, id = "hhidpn", time = "wave", trt_time = "treatment_time",
                     options = list(time_lag = TRUE, optimizer = "glpk"))

  expect_equivalent(pairs1[, 2:3], pairs2[, 2:3])
  expect_equivalent(as.character(pairs1$hhidpn), pairs2$hhidpn)
})

test_that("`brsmatch()` returns warning when 'trt_time' is not numeric", {
  df <- data.frame(
    hhidpn = rep(1:3, each = 3),
    wave = rep(1:3, 3),
    treatment_time = rep(c(2,3,NA), each = 3),
    X1 = c(2,2,2,3,3,3,9,9,9),
    X2 = rep(c("a","a","b"), each = 3),
    X3 = c(9,4,5,6,7,2,3,4,8),
    X4 = c(8,9,4,5,6,7,2,3,4)
  )

  check_for_glpk()
  pairs1 <- brsmatch(n_pairs = 1, data = df, id = "hhidpn", time = "wave", trt_time = "treatment_time",
                     options = list(time_lag = TRUE, optimizer = "glpk"))

  df$treatment_time <- as.character(df$treatment_time)

  expect_warning({
    pairs2 <- brsmatch(n_pairs = 1, data = df, id = "hhidpn", time = "wave", trt_time = "treatment_time",
                       options = list(time_lag = TRUE, optimizer = "glpk"))
  })

  expect_equivalent(pairs1, pairs2)
})


test_that("`brsmatch()` works when there are no never-treated individuals", {
  df1 <- data.frame(
    hhidpn = rep(1:5, each = 7),
    wave = rep(1:7, 5),
    treatment_time = rep(c(2,3,3,4,7), each = 7),
    X1 = c(2,2,4,5,5,5,4,
           9,9,10,10,10,7,7,
           2,3,4,5,6,6,7,
           4,5,6,6,6,5,1,
           3,5,6,6,7,5,6),
    X2 = rep(c("a","a","b","c","d"), each = 7),
    X3 = c(9,4,5,6,7,2,3,
           4,8,5,7,8,5,8,
           7,4,5,6,7,7,8,
           4,5,6,7,8,9,7,
           5,6,7,5,6,5,5),
    X4 = c(8,9,4,5,6,7,2,
           3,4,6,4,2,5,7,
           3,3,4,6,2,4,5,
           3,5,6,3,4,3,3,
           3,2,3,3,5,6,3)
  )
  df2 <- df1
  df2[df2$treatment_time == 7, "treatment_time"] <- NA

  dist1 <- .compute_distances(df1, id = "hhidpn", time = "wave", trt_time = "treatment_time")
  dist2 <- .compute_distances(df2, id = "hhidpn", time = "wave", trt_time = "treatment_time")
  expect_equal(dist1, dist2)

  check_for_glpk()
  pairs <- brsmatch(n_pairs = 2, data = df1,
                    id = "hhidpn", time = "wave", trt_time = "treatment_time")
})

test_that("`brsmatch() works with exact_match covariates", {

  check_exact_match <- function(pairs, n_pairs) {
    tmp <- oasis %>%
      left_join(pairs, by = "subject_id") %>%
      distinct(subject_id, m_f, pair_id)

    tbl <- table(tmp$m_f, tmp$pair_id)
    expect_setequal(unique(as.numeric(tbl)), c(0, 2))
    expect_equal(sum(rowSums(tbl) / 2), n_pairs)
    expect_equal(unique(colSums(tbl)), 2)
  }

  pairs <- brsmatch(n_pairs = 5, oasis,
                    id = "subject_id", time = "visit", trt_time = "time_of_ad",
                    balance = FALSE,
                    exact_match = c("m_f"))

  check_exact_match(pairs, 5)

  pairs <- brsmatch(n_pairs = 1, oasis,
                    id = "subject_id", time = "visit", trt_time = "time_of_ad",
                    balance = FALSE,
                    exact_match = c("m_f"))

  check_exact_match(pairs, 1)

  pairs <- brsmatch(n_pairs = 10, oasis,
                    id = "subject_id", time = "visit", trt_time = "time_of_ad",
                    balance = FALSE,
                    exact_match = c("m_f"))

  check_exact_match(pairs, 10)

})
