#' Compute distance on valid matches in Risk Set Matching.
#'
#' TODO: Come up with a better description here.  Mention the calculation heuristics from the Li et al. (2001) Balanced Risk Set Matching
#' paper.
#'
#' @param df data frame containing columns matching the \code{id, time, trt_time} arguments, and covariates.
#'   This data frame is expected to be in tidy, long format, so that the \code{id}, \code{trt_time}, and
#'   baseline  variables may be repeated for different values of \code{time}.  Data frame should be unique
#'   at \code{id} and \code{time}.
#' @param id optional parameter to specify the name of the id column in \code{df}.
#' @param time optional parameter to specify the name of the time column in \code{df}.
#' @param trt_time optional parameter to specify the name of the treatment time column in \code{df}.
#' @param covariates optional parameter to specify the names of covariates to be used. If \code{NULL}, will
#'   default to all columns except those named by \code{id, time, trt_time}.
#'
#' @return a data frame with valid risk set matching pairs and their corresponding distance.  This data frame
#'   will have four columns.  \code{trt_id} refers to the treated id, \code{all_id} refers to the id that is
#'   a possible control, \code{trt_time} refers to the treatment time of \code{trt_id}, and \code{dist} refers to
#'   the Mahalanobis distance
#'
#' @examples
#' df <- data.frame(
#'   hhidpn = rep(1:3, each = 3),
#'   wave = rep(1:3, 3),
#'   treatment_time = rep(c(2,3,NA), each = 3),
#'   X1 = c(2,2,2,3,3,3,9,9,9),
#'   X2 = rep(c("a","a","b"), each = 3),
#'   X3 = c(9,4,5,6,7,2,3,4,8),
#'   X4 = c(8,9,4,5,6,7,2,3,4)
#' )
#'
#' enumerate_edges_and_compute_distances(df, "hhidpn", "wave", "treatment_time")
#'
#' @export
enumerate_edges_and_compute_distances <- function(df, id = "id", time = "time", trt_time = "trt_time", covariates = NULL) {
  # TODO: check that the all units treatment time is not needed in this (it shows up in edges in the code)

  if (is.null(covariates)) {
    covariates <- setdiff(colnames(df), c(id, time, trt_time))
  }
  df[[trt_time]][which(is.na(df[[trt_time]]))] <- 0
  # TODO: check treatment times are either 0, NA, or match times

  # pre-compute covariance matrix
  cov_mat <- cov(model.matrix(~ 0 + ., data = df[, covariates]))
  cov_mat_inv <- MASS::ginv(cov_mat)

  split_id <- split(df, factor(df[[id]]))
  ids <- unlist(lapply(split_id, function(x) unique(x[[id]])))
  trt_times <- unlist(lapply(split_id, function(x) unique(x[[trt_time]])))
  treated_ids <- ids[trt_times > 0]

  # Iterate over treated ids, calculate Mahalanobis dist, and put into a data.frame
  out <- list()
  for (i in treated_ids) {
    trt_time_i <- trt_times[[which(ids == i)]]
    df_at_trt <- df[df[[time]] == trt_time_i - 1, ] # TODO: rename to reflect this is pre-treatment data

    if(i %in% df_at_trt[[id]]) {
      covariates_at_trt <- model.matrix(~ 0 + ., data = df_at_trt[, covariates])
      dists <- stats::mahalanobis(covariates_at_trt,
                                  covariates_at_trt[which(df_at_trt[[id]] == i),],
                                  cov = cov_mat_inv,
                                  inverted = TRUE)
      valid_match <- df_at_trt[[id]] != i & # can't match control with itself
        (df_at_trt[[trt_time]] > trt_time_i | df_at_trt[[trt_time]] == 0) # control receives treatment later, or not at all

      out[[i]] <- data.frame(
        trt_id = i,
        all_id = df_at_trt[[id]][valid_match],
        trt_time = trt_time_i,
        dist = dists[valid_match]
      )
    }
  }
  do.call(rbind, out)
}

#' Compute balance covariate indicators
#'
#' TODO: add description
#'
#' @param df data frame containing columns matching the \code{id, time, trt_time} arguments, and covariates.
#'   This data frame is expected to be in tidy, long format, so that the \code{id}, \code{trt_time}, and
#'   baseline  variables may be repeated for different values of \code{time}.  Data frame should be unique
#'   at \code{id} and \code{time}.
#' @param id optional parameter to specify the name of the id column in \code{df}.
#' @param time optional parameter to specify the name of the time column in \code{df}.
#' @param trt_time optional parameter to specify the name of the treatment time column in \code{df}.
#' @param balance_covariates optional parameter to specify the names of covariates to be used in balancing.
#'   If \code{NULL}, will default to all columns except those named by \code{id, time, trt_time}.
#'
#' @return TODO: add return description
#'
#' @examples
#' df <- data.frame(
#    hhidpn = rep(1:3, each = 3),
#    wave = rep(1:3, 3),
#'   treatment_time = rep(c(2,3,NA), each = 3),
#'   X1 = c(2,2,2,3,3,3,9,9,9),
#'   X2 = rep(c("a","a","b"), each = 3),
#'   X3 = c(9,4,5,6,7,2,3,4,8),
#'   X4 = c(8,9,4,5,6,7,2,3,4)
#' )
#'
#' balance_covariates <- c("X1", "X2", "X3", "X4")
#' bal <- balance_columns(df, "hhidpn", "wave", "treatment_time", balance_covariates = balance_covariates)
#'
#' @export
balance_columns <- function(df, id = "id", time = "time", trt_time = "trt_time", balance_covariates = NULL) {
  if (is.null(balance_covariates)) {
    balance_covariates <- setdiff(colnames(df), c(id, time, trt_time))
  }
  df[[trt_time]][which(is.na(df[[trt_time]]))] <- 0

  # split balance_covariates into types
  bal_cov_types <- sapply(df[, balance_covariates], class)
  numeric_cov <- balance_covariates[which(bal_cov_types %in% c("numeric", "integer"))]
  factor_cov <- balance_covariates[which(bal_cov_types %in% c("factor", "character"))]
  other_cov <- setdiff(balance_covariates, union(numeric_cov, factor_cov))
  if (length(other_cov) > 0) {
    warning(paste0("Columns ", other_cov, " are of unusable type and will be omitted.  Useable types are numeric, integer, factor, and character"))
  }

  # calculate quantiles based on treated ids at treatment time
  trt_at_trt_time_df <- df[df[[time]] == df[[trt_time]], ]
  quantiles <- apply(trt_at_trt_time_df[, numeric_cov],
                     MARGIN = 2,
                     quantile,
                     probs = c(1/3, 2/3))
  # balance on numeric columns
  out <- list()
  for (col in colnames(quantiles)) {
    for (row in 1:nrow(quantiles)) {
      out[[paste0(col, ".q",row)]] <- 1*(df[[col]] > quantiles[row, col])
    }
  }
  bal_numeric <- do.call(cbind, out)
  # balance on character and factor columns
  bal_factor <- model.matrix(~ 0 + ., data = df[,factor_cov, drop = FALSE])

  cbind(
    id = df[[id]],
    time = df[[time]],
    bal_factor,
    bal_numeric
  )
}









