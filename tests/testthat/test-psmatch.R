context("Test the functions used in psmatch.R")

df <- data.frame(
  hhidpn = rep(1:3, each = 3),
  wave = rep(1:3, 3),
  treatment_time = rep(c(2,3,NA), each = 3),
  X1 = c(2,2,2,3,3,3,9,9,9),
  X2 = rep(c("a","a","b"), each = 3),
  X3 = c(9,4,5,6,7,2,3,4,8),
  X4 = c(8,9,4,5,6,7,2,3,4)
)

test_that("`coxph_match()` has correct output", {

  expect_warning({
    pairs <- coxph_match(n_pairs = 1, df = df,
                         id = "hhidpn",
                         time = "wave",
                         trt_time = "treatment_time")
  })
  expect_equal(colnames(pairs), c("hhidpn", "pair_id", "type"))
  expect_equal(length(unique(na.omit(pairs$pair_id))), 1)
  expect_equal(unique(pairs$hhidpn), unique(df$hhidpn))
  expect_equal(pairs$hhidpn[which(pairs$pair_id == 1)], c(2,3))

  # check runs properly with other arguments
  expect_warning({
    coxph_match(n_pairs = 1, df = df, id = "hhidpn", time = "wave", trt_time = "treatment_time")
  })

})


test_that("`coxph_match()` works when 'id' is a character vector", {
  expect_warning({
    pairs1 <- coxph_match(n_pairs = 1, df = df, id = "hhidpn",
                        time = "wave", trt_time = "treatment_time")
  })

  df$hhidpn <- as.character(df$hhidpn)

  expect_warning({
    pairs2 <- coxph_match(n_pairs = 1, df = df, id = "hhidpn",
                          time = "wave", trt_time = "treatment_time")
  })

  expect_equivalent(pairs1[, 2:3], pairs2[, 2:3])
  expect_equivalent(as.character(pairs1$hhidpn), pairs2$hhidpn)

  df$hhidpn <- as.numeric(df$hhidpn)
})

test_that("`coxph_match()` returns warning when 'trt_time' is not numeric", {

  expect_warning({
    pairs1 <- coxph_match(n_pairs = 1, df = df, id = "hhidpn",
                          time = "wave", trt_time = "treatment_time")
  }, "ghost value")

  df$treatment_time <- as.character(df$treatment_time)

  expect_warning({
    pairs2 <- coxph_match(n_pairs = 1, df = df, id = "hhidpn",
                          time = "wave", trt_time = "treatment_time")
  }, "should be numeric")

  expect_equivalent(pairs1, pairs2)
})


test_that("`coxph_match()` works when there are no never-treated individuals", {
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

  expect_warning({
    pairs <- coxph_match(n_pairs = 2, df = df1, id = "hhidpn", time = "wave",
                         trt_time = "treatment_time")
  }, "ghost value")

})


