#' Not in, as opposed to %in%
#' @noRd
`%ni%` <- Negate(`%in%`)

#' Output pairs to new format
#'
#' Takes a data frame with each row as a pair and returns output in long, tidy
#' format that indicates the matched pairs
#'
#' @inheritParams brsmatch
#' @param matched_ids A data.frame with two columns: `trt_id, all_id`.  Each row
#'   consists of a matched pair, where "trt_id" provides the id of the treated
#'   and "all_id" provides the id of the control.
#' @param id_list An optional vector of ids to include in the output.
#'
#' @return a data frame with columns "id", "pair_id", and "type".  "id" refers
#'   to the individual ids, "pair_id" is a unique identifier for each pair, and
#'   "type" indicates whether the id is from treatment ("trt") or control ("all")
#'
#' @examples
#' matched_ids <- data.frame(
#'   trt_id = c(2, 5, 7),
#'   all_id = c(3, 1, 4)
#' )
#'
#' output_pairs(matched_ids)
#'
#' @noRd
.output_pairs <- function(matched_ids, id = "id", id_list = NULL) {
  if (is.null(id_list)) id_list <- unlist(matched_ids, use.names = FALSE)
  pairs <- data.frame(id = id_list,
                      pair_id = NA,
                      type = NA)
  for (rowid in 1:nrow(matched_ids)) {
    match_ind <- pairs$id %in% matched_ids[rowid, ]
    pairs$pair_id[match_ind] <- rowid
  }
  trt_ind <- pairs$id %in% matched_ids[, "trt_id"]
  all_ind <- pairs$id %in% matched_ids[, "all_id"]
  pairs$type[trt_ind] <- "trt"
  pairs$type[all_ind] <- "all"
  names(pairs)[1] <- id
  return(pairs)
}



#' Split an integer according to weights
#'
#' Splits `n` according to `weights` in the most even way possible.  By this, we
#' mean that numbers closest to the 'next digit' when considering `n * weights`
#' will be rounded up or down first.  The procedure is kind of complicated, but
#' I'm fairly sure that it works.
#'
#' @param n An integer.
#' @param weights A set of weights which determine how many groups to split `n`
#'   into, and their relative strengths
#'
#' @return A vector of size `length(weights)` which should sum to `n`.
#'
#' @noRd
.weighted_split <- function(n, weights) {
  weights <- weights / sum(weights)

  x <- n*weights
  out <- round(x)
  diff <- x - round(x)
  n_left <- n - sum(out)

  if (n_left == 0) {
    return(out)
  } else if (n_left < 0) {
    # overestimated, need to subtract
    diff[diff < 0] <- diff[diff < 0] + 1
    for (i in 1:-n_left) {
      add_to <- which.min(diff)
      out[add_to] <- out[add_to] - 1
      diff[add_to] <- diff[add_to] + 1 # ensure it won't be the min again
    }
  } else {
    # underestimated, need to add
    diff[diff > 0] <- diff[diff > 0] - 1
    for (i in 1:n_left) {
      add_to <- which.max(diff)
      out[add_to] <- out[add_to] + 1
      diff[add_to] <- diff[add_to] - 1 # ensure it won't be the min again
    }
  }
  stopifnot(sum(out) == n)
  return(out)
}
