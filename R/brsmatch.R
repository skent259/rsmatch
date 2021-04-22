#' Balanced Risk Set Matching
#'
#' Perform balanced risk set matching as described in Li et al. (2001) "Balanced
#' Risk Set Matching".  Given a longitudinal data frame with covariate
#' information, along with treatment time, build a MIP problem that matches
#' treated individuals to those that haven't been treated yet (or are never
#' treated) based on minimizing the Mahalanobis distance between covariates. If
#' balancing is desired, the model will try to minimize the imbalance in terms
#' of specified balancing covariates in the final pair output.  Each treated
#' individual is matched to one other individual.
#'
#' @param n_pairs number of pairs desired from matching
#' @param df data frame containing columns matching the \code{id, time,
#'   trt_time} arguments, and covariates. This data frame is expected to be in
#'   tidy, long format, so that the \code{id}, \code{trt_time}, and baseline
#'   variables may be repeated for different values of \code{time}.  Data frame
#'   should be unique at \code{id} and \code{time}.
#' @param id optional parameter to specify the name of the id column in
#'   \code{df}.
#' @param time optional parameter to specify the name of the time column in
#'   \code{df}.
#' @param trt_time optional parameter to specify the name of the treatment time
#'   column in \code{df}.
#' @param covariates optional parameter to specify the names of covariates to be
#'   used. If \code{NULL}, will default to all columns except those named by
#'   \code{id, time, trt_time}.
#' @param balance_covariates optional parameter to specify the names of
#'   covariates to be used in balancing. If \code{NULL}, will default to all
#'   columns except those named by \code{id, time, trt_time}.
#' @param optimizer "gurobi" or "glpk". Specifies which optimizer output to use;
#'   defaults to "gurobi".
#' @param verbose logical; if TRUE, will print some useful information for
#'   potentially long model calls; defaults to FALSE
#' @param balance When TRUE, a balanced risk set matching model will be built.
#'   When FALSE, or when bal_all = NULL, balancing constraints will not be
#'   included.
#' @param options additional arguments for rare scenarios, default = "none".  If
#'   options = 'between period treatment' then matches must be made on the time
#'   period preceding treatment and some special care must be taken.
#'
#' @return a data frame containing the pair information.  The data frame has
#'   columns \code{id}, "pair_id", and "type". \code{id} matches the input
#'   parameter and will contain all ids from the input data frame.  "pair_id"
#'   refers to the id of the computed pairs; NA values indicate unmatched
#'   individuals.  "type" indicates whether the individual in the pair is
#'   considered as treatment ("trt") or control ("all") in that pair.
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
#' brsmatch(n_pairs = 1, df = df, id = "hhidpn", time = "wave",
#'          trt_time = "treatment_time", optimizer = "glpk")
#'
#' @export
#' @author Sean Kent
brsmatch <- function(n_pairs,
                     df,
                     id = "id", time = "time", trt_time = "trt_time",
                     covariates = NULL, balance_covariates = NULL,
                     optimizer = "gurobi", verbose = FALSE, balance = TRUE,
                     options = c("none", "between period treatment")) {

  options <- match.arg(options)
  optimizer <- match.arg(optimizer, c("gurobi", "glpk"))
  if (optimizer == "gurobi" & !requireNamespace("gurobi", quietly = TRUE)) {
    rlang::abort(c(
      "Package 'gurobi' required when `optimizer == 'gurobi'`.",
      i = "This package requires gurobi to be installed on your computer.",
      i = "If you have gurobi installed, see https://www.gurobi.com/documentation/9.1/refman/ins_the_r_package.html for package installation. "
    ))
  } else if (optimizer == "glpk" & !requireNamespace("Rglpk", quietly = TRUE)) {
    rlang::abort(c(
      "Package 'Rglpk' required when `optimizer == 'glpk'`.",
      i = "Please install the package and retry the funciton."
    ))
  }
  if (!is.numeric(df[[trt_time]])) {
    rlang::warn(c(
      paste0("Treatment time `", trt_time, "` should be numeric."),
      i = "Converting to a numeric column."
    ))
    df[[trt_time]] <- as.numeric(df[[trt_time]])
  }
  if (options == "between period treatment") {
    # need to match on time just before treatment
    df[[trt_time]] <- df[[trt_time]] - 1
  }
  if (verbose) message("Computing distances from df...")
  edges <- compute_distances(df, id, time, trt_time, covariates, options)
  bal <- NULL
  if (balance) {
    if (verbose) message("Building balance columns from df...")
    bal <- balance_columns(df, id, time, trt_time, balance_covariates)
  }
  if (verbose) message("Constructing optimization model...")
  model <- rsm_optimization_model(n_pairs, edges, bal, optimizer, verbose, balance)

  if (verbose) message("Preparing to run optimization model")
  if (optimizer == "gurobi") {
    res <- gurobi::gurobi(model)
    matches <- res$x[grepl("f", model$varnames)]
  } else if (optimizer == "glpk") {
    res <- with(model, Rglpk::Rglpk_solve_LP(obj, mat, dir, rhs, types = types, max = max,
                                             control = list(verbose = verbose, presolve = TRUE)))
    matches <- res$solution[grepl("f", model$varnames)]
  }

  matched_ids <- edges[matches == 1, c("trt_id", "all_id")]
  output_pairs(matched_ids, id = id, id_list = unique(df[[id]]))
}

#' Compute distance on valid matches in Risk Set Matching.
#'
#' The \code{compute_distances} function takes in longitudinal data and
#' computes the Mahalanobis distance between eligible edges according to the
#' procedure of risk set matching. In risk set matching, each individual that
#' is treated at time t can be matched to someone who hasn't been treated yet
#' as of time t.  The Mahalanobis distance is computed for each treated
#' id against possible matches based on baseline and timevarying covariates at
#' time t. See Li et al. (2001) "Balanced Risk Set Matching" for additional
#' details.
#'
#' @inheritParams brsmatch
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
#' compute_distances(df, "hhidpn", "wave", "treatment_time")
#'
#' @export
compute_distances <- function(df, id = "id", time = "time", trt_time = "trt_time", covariates = NULL, options = "none") {
  if (is.null(covariates)) {
    covariates <- setdiff(colnames(df), c(id, time, trt_time))
  }
  df[[trt_time]][which(is.na(df[[trt_time]]))] <- 0
  # TODO: check treatment times are either 0, NA, or match times

  # pre-compute covariance matrix
  cov_mat <- stats::cov(stats::model.matrix(~ 0 + ., data = df[, covariates]))
  cov_mat_inv <- MASS::ginv(cov_mat)

  split_id <- split(df, factor(df[[id]]))
  ids <- unlist(lapply(split_id, function(x) unique(x[[id]])))
  trt_times <- unlist(lapply(split_id, function(x) unique(x[[trt_time]])))
  treated_ids <- ids[trt_times > 0]

  # Iterate over treated ids, calculate Mahalanobis dist, and put into a data.frame
  out <- lapply(1:length(treated_ids), FUN = function(j) {
    i <- treated_ids[[j]]
    trt_time_i <- trt_times[[which(ids == i)]]
    df_at_trt <- df[df[[time]] == trt_time_i, ]

    if (i %in% df_at_trt[[id]]) {
      covariates_at_trt <- stats::model.matrix(~ 0 + ., data = df_at_trt[, covariates])
      dists <- stats::mahalanobis(covariates_at_trt,
                                  covariates_at_trt[which(df_at_trt[[id]] == i),],
                                  cov = cov_mat_inv,
                                  inverted = TRUE)
      valid_match <- df_at_trt[[id]] != i & # can't match control with itself
        (df_at_trt[[trt_time]] > trt_time_i | df_at_trt[[trt_time]] == 0) # control receives treatment later, or not at all
      if (options == "between period treatment") {
        # potential matches must also exist at trt_time (unless trt_time is the last period)
        exist_after_trt <- df_at_trt[[id]] %in% df[df[[time]] == trt_time_i + 1, ][[id]]
        valid_match <- valid_match & (exist_after_trt | trt_time_i == max(df[[time]]) )
      }

      out_j <- data.frame(
        trt_id = i,
        all_id = df_at_trt[[id]],
        trt_time = trt_time_i,
        dist = dists
      )
      return(out_j[valid_match, , drop = FALSE])
    } else {
      return(NULL)
    }
  })
  do.call(rbind, out)
}

#' Compute balance covariate indicators
#'
#' The \code{balance_columns} function takes longitudinal data input and returns
#' corresponding balance columns on desired covariates according to the process
#' described in Li et al. (2001) "Balanced Risk Set Matching". Each balance
#' column is an indicator variable.
#'
#' For numeric covariates, quantiles are computed for treated individuals at
#' their treated time, and indicator variables for whether each covariate
#' exceeds the quantiles are returned.  Default quantiles used are 1/3 and 2/3.
#'
#' For factor and character covariates, indicators for each level (except one)
#' are returned as would be done in a call to \code{stats::model.matrix}.
#'
#' @inheritParams brsmatch
#'
#' @return a matrix with the same number of rows as the input dataframe.  Each
#'   column is an indicator variable corresponds to a balancing criteria.  The
#'   number of balance columns returned is 2*n_numeric_cols + n_factor_cols.
#'
#' @examples
#' df <- data.frame(
#'    hhidpn = rep(1:3, each = 3),
#'    wave = rep(1:3, 3),
#'   treatment_time = rep(c(2,3,NA), each = 3),
#'   X1 = c(2,2,2,3,3,3,9,9,9),
#'   X2 = rep(c("a","a","b"), each = 3),
#'   X3 = c(9,4,5,6,7,2,3,4,8),
#'   X4 = c(8,9,4,5,6,7,2,3,4)
#' )
#'
#' balance_covariates <- c("X1", "X2", "X3", "X4")
#' bal <- balance_columns(df, "hhidpn", "wave", "treatment_time",
#'                        balance_covariates = balance_covariates)
#'
#' @export
balance_columns <- function(df, id = "id", time = "time", trt_time = "trt_time", balance_covariates = NULL) {
  if (is.null(balance_covariates)) {
    balance_covariates <- setdiff(colnames(df), c(id, time, trt_time))
  }
  df[[trt_time]][which(is.na(df[[trt_time]]))] <- 0

  # split balance_covariates into types
  bal_cov_types <- sapply(df[, balance_covariates, drop = FALSE], class)
  numeric_cov <- balance_covariates[which(bal_cov_types %in% c("numeric", "integer"))]
  factor_cov <- balance_covariates[which(bal_cov_types %in% c("factor", "character"))]
  other_cov <- setdiff(balance_covariates, union(numeric_cov, factor_cov))
  if (length(other_cov) > 0) {
    warning(paste0("Columns ", other_cov, " are of unusable type and will be omitted.  Useable types are numeric, integer, factor, and character"))
  }

  empty_df <- matrix(nrow = nrow(df), ncol = 0)
  # calculate quantiles based on treated ids at treatment time
  if (length(numeric_cov) > 0) {
    trt_at_trt_time_df <- df[df[[time]] == df[[trt_time]], ]
    quantiles <- apply(trt_at_trt_time_df[, numeric_cov, drop = FALSE],
                       MARGIN = 2,
                       stats::quantile,
                       probs = c(1/3, 2/3))
    # balance on numeric columns
    out <- list()
    for (col in colnames(quantiles)) {
      for (row in 1:nrow(quantiles)) {
        out[[paste0(col, ".q",row)]] <- 1*(df[[col]] > quantiles[row, col])
      }
    }
    bal_numeric <- do.call(cbind, out)
  } else {
    bal_numeric <- empty_df
  }
  # balance on character and factor columns
  if (length(factor_cov) > 0) {
    bal_factor <- stats::model.matrix(~ 0 + ., data = df[,factor_cov, drop = FALSE])
  } else {
    bal_factor <- empty_df
  }

  res <- cbind(
    df[, c(id, time)],
    bal_factor,
    bal_numeric
  )
  colnames(res)[1:2] <- c("id", "time")
  return(res)
}


#' Build optimization model for Balanced Risk Set Matching
#'
#' The \code{rsm_optimization_model} function takes in a distance data-frame
#' (tidy-format) and optional balancing columns and returns an optimization
#' model in the corresponding optimizer format.  The model is build to generate
#' \code{n_pairs} pairs that minimize the distance while ensuring that each id
#' gets used in at most one pair.  If \code{balance} is TRUE, the optimization
#' model will also add constraint variables that heavily penalize violations
#' of the balance criterion.  This function is rarely useful on its own, it
#' is preferred to make a function call to \code{brsmatch}.
#'
#' @inheritParams brsmatch
#' @param edges data frame with columns "trt_id", "all_id", "trt_time", "dist";
#'   for example, the output from a call to \code{compute_distances()}
#' @param bal_all matrix with columns "id", "time", and additional balance
#'   columns; for example, the output from a call to \code{balance_columns()};
#'   defaults to NULL, indicating balance is not used.
#'
#' @return an optimization model that can be readily passed to the optimizer
#'   parameter.  Defines the mixed integer programming problem for risk set
#'   matching in terms of specified data.
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
#' edges <- compute_distances(df, "hhidpn", "wave", "treatment_time")
#' bal <- balance_columns(df, "hhidpn", "wave", "treatment_time")
#' n_unique_id <- length(unique(df$hhidpn))
#'
#' model <- rsm_optimization_model(1, edges, bal, optimizer = "gurobi", balance = TRUE)
#'
#' @export
rsm_optimization_model <- function(n_pairs,
                                   edges,
                                   bal_all = NULL,
                                   optimizer = "gurobi",
                                   verbose = FALSE,
                                   balance = TRUE) {

  if (is.null(bal_all)) balance <- FALSE
  if (balance) {
    # TODO: check that no columns have the .trt or .all name in them already. This is unlikely
    bal_all <- as.data.frame(bal_all)
    # bal_all$time <- bal_all$time + 1 # want to match when time == trt_time - 1
    edges$.rowid <- 1:nrow(edges)
    edges <- merge(edges, bal_all, by.x = c("trt_id", "trt_time"), by.y = c("id", "time"))
    edges <- merge(edges, bal_all, by.x = c("all_id", "trt_time"), by.y = c("id", "time"),
                   suffixes = c(".trt", ".all"))
    edges <- edges[order(edges$.rowid), ]
    edges$.rowid <- NULL
  }

  S <- n_pairs # number of pairs
  K <- ifelse(balance, ncol(bal_all)-2, 0) # number of balancing constraints TODO: maybe rename this as n_gp, n_gm
  E <- nrow(edges) # number of edges  TODO: maybe rename this as n_f
  n_vars <- E + 2*K

  delta <- edges$dist
  lambda_k <- sum(delta) + 100

  ## A.2 - There are S matched sets
  A.2 <- c(rep(1, E), rep(0, 2*K))
  A.2 <- rbind(A.2, -A.2)

  ## A.3 - each unit has at most one edge
  # all_ids <- unique(edges$all_id)
  all_ids <- unique(c(edges$all_id, edges$trt_id))
  edge_ids <- edges[, c("trt_id", "all_id")]

  j_ind <- lapply(all_ids, function(id) {
    # TODO: might be able to vectorize this?
    # if (verbose & (id %% 50 == 0)) { cat("  i:", id, "/", length(all_ids), "\n") }
    # id <- all_ids[x]
    # pairs_with_hhidpn_i <- union(
    #   which(tmp$treated_units == hhidpn_i),
    #   which(tmp$all_units == hhidpn_i)
    # )
    which(rowSums(edge_ids == id) == 1) # rows where either treated == hhidpn or all == hhidpn
    # TODO: use microbenchmark to try this a few different ways.
  })
  i_ind <- mapply(
    function(x,y) {
      rep.int(y, times = length(x))
    },
    j_ind,
    1:length(j_ind),
    SIMPLIFY = FALSE
  )
  j_ind <- do.call(c, j_ind)
  i_ind <- do.call(c, i_ind)

  A.3 <- Matrix::sparseMatrix(i_ind, j_ind,
                              dims = c(length(all_ids), n_vars))

  # # I think this is faster....
  # TODO: check if this is faster
  # sparse_inds <- lapply(1:length(all_ids), function(i) {
  #   # TODO: might be able to vectorize this?
  #   if (verbose & (i %% 50 == 0)) { cat("  i:", i, "/", length(all_ids), "\n") }
  #   id <- all_ids[i]
  #   # pairs_with_hhidpn_i <- union(
  #   #   which(tmp$treated_units == hhidpn_i),
  #   #   which(tmp$all_units == hhidpn_i)
  #   # )
  #   active_rows <- which(rowSums(edge_ids == id) == 1) # rows where either treated == hhidpn or all == hhidpn
  #   rbind(i = rep(i, length(active_rows)),
  #         j = active_rows)
  #
  #   # TODO: use microbenchmark to try this a few different ways.
  # })
  # sparse_inds <- do.call(cbind, sparse_inds)
  # A.3 <- Matrix::sparseMatrix(sparse_inds["i",], sparse_inds["j",],
  #                             # x = 1,
  #                             dims = c(length(all_ids), n_vars))
  #

  ## A.4 and A.5 - balance constraints
  if (balance) {
    B_e <- as.matrix(edges[, grepl(".trt", names(edges))])
    B_p <- as.matrix(edges[, grepl(".all", names(edges))])
    f_e_conditions <- t(B_p) - t(B_e)
    I_minus_k <- diag(rep(-1, K)) # I_{k \by k}
    zero_k <- matrix(0, nrow = K, ncol = K) # 0_{k \by k}

    A.45 <- rbind(
      cbind(f_e_conditions, I_minus_k, zero_k),
      cbind(-f_e_conditions, zero_k, I_minus_k)
    )
  } else {
    A.45 <- NULL
  }

  model <- list()
  if (balance) {
    ## Objective
    model$modelsense <- "min"
    model$obj <- c(delta, rep(lambda_k, 2*K))
    ## Constraints
    model$varnames <- c(paste0("f",1:E), paste0("gp",1:K), paste0("gm",1:K))
    model$A <- rbind(A.2, A.3, A.45)
    model$sense <- rep("<=", nrow(model$A))
    model$rhs <- c(S, -S, rep(1, length(all_ids)), rep(0, 2*K))
    model$vtype <- c(rep("B", E), rep("C", 2*K))
  } else {
    ## Objective
    model$modelsense <- "min"
    model$obj <- c(delta)
    ## Constraints
    model$varnames <- c(paste0("f",1:E))
    model$A <- rbind(A.2, A.3)
    model$sense <- rep("<=", nrow(model$A))
    model$rhs <- c(S, -S, rep(1, length(all_ids)))
    model$vtype <- c(rep("B", E))
  }

  if (optimizer == "gurobi") {
    return(model)
  } else  if (optimizer == "glpk") {
    names(model) <- c("max", "obj", "varnames", "mat", "dir", "rhs", "types")
    model$max <- ifelse(model$max == "min", FALSE, NA)
    return(model)
  } else {
    stop(paste0("Model optimizer of type ", optimizer, "is not supported.  Use 'gurobi' or 'glpk'."))
  }
  return(model)
}

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
#' @export
output_pairs <- function(matched_ids, id = "id", id_list = NULL) {
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








