#' Output pairs to new format
#'
#' Takes a data frame with each row as a pair and returns output in long, tidy
#' format that indicates the matched pairs
#'
#' @inheritParams brsmatch
#' @param matched_ids data frame with two columns: trt_id, all_id.  Each row
#'   consists of a matched pair, where "trt_id" provides the id of the treated
#'   and "all_id" provides the id of the control.
#' @param id_list optional vector of ids to include in the output.
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
