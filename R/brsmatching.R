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

