#' Propensity Score Matching with Time-Dependent Covariates
#'
#' Perform propensity score matching as described in Lu (2005) "Propensity Score
#' Matching with Time-Dependent Covariates".  Given a longitudinal data frame
#' with covariate information, along with treatment time, build a MIP problem
#' that matches treated individuals to those that haven't been treated yet (or
#' are never treated) based on time-dependent propensity scores from a Cox
#' proportional hazards model. If balancing is desired, the model will try to
#' minimize the imbalance in terms of specified balancing covariates in the
#' final pair output.  Each treated individual is matched to one other
#' individual.
#'
#' @inheritParams brsmatch
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
#' coxph_match(n_pairs = 1, df = df, id = "hhidpn", time = "wave", trt_time = "treatment_time")
#'
#' @export
#' @author Mitchell Paukner
coxph_match <- function(n_pairs = 10^10,
                        df,
                        id = "id",
                        time = "time",
                        trt_time = "trt_time",
                        covariates = NULL,
                        balance_covariates = NULL){

  if (!is.null(balance_covariates)) {
    balance_split <- split(df, df[, balance_covariates])
    matches <- NULL

    for (i in 1:length(balance_split)) {
      matches <- rbind(matches,coxph_match1(balance_split[[i]],id, time,
                                            trt_time, trt, covariates))
    }
  } else {
    matches <- coxph_match1(df, id, time, trt_time, covariates)
  }

  matches <- matches[order(matches$distance), ]

  if(length(matches$trt.id) > n_pairs) {
    matches <- matches[1:n_pairs,]
  }

  matches$pair_id <- c(1:dim(matches)[1])

  matches <- tidyr::pivot_longer(matches, c("trt.id","con.id"),
                                 names_to = "type", values_to = "id")

  final <- data.frame("id" = matches$id, "pair_id" = matches$pair_id, "type" = matches$type)
  final$type <- ifelse(final$type == "trt.id", "trt", "all")

  return(final)
}

#' Propensity Score Matching with Time-Dependent Covariates
#'
#' Perform propensity score matching as described in Lu (2005) "Propensity Score
#' Matching with Time-Dependent Covariates".  Given a longitudinal data frame
#' with covariate information, along with treatment time, build a MIP problem
#' that matches treated individuals to those that haven't been treated yet (or
#' are never treated) based on time-dependent propensity scores from a Cox
#' proportional hazards model. If balancing is desired, the model will try to
#' minimize the imbalance in terms of specified balancing covariates in the
#' final pair output.  Each treated individual is matched to one other
#' individual.
#'
#' @inheritParams brsmatch
#' @return a data frame with matched pairs based on propensity scores generated
#'   using a Cox proportional hazards model. This data frame will have three
#'   columns.  \code{trt.id} refers to the treated id, \code{con.id} refers to
#'   the id that is a possible control, \code{distance} refers to the scaled
#'   euclidean distance generated based on the propensity scores.
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
#' coxph_match1(df = df, id = "hhidpn", time = "wave", trt_time = "treatment_time")
#'
#' @export
coxph_match1 <- function(df,
                         id = "id",
                         time = "time",
                         trt_time = "trt_time",
                         covariates = NULL) {
  if(is.null(covariates)){
    data <- df[,-which(names(df) %in% c(id, time, trt_time))]
  } else{
    data <- df[, c(covariates)]
  }
  data$id <- as.factor(unlist(df[, id]))
  data$time <- unlist(df[, time])
  data$trt_time <- unlist(df[, trt_time])
  data$trt <- ifelse(is.na(data$trt_time), 0, 1)
  data[, covariates] <-
    data[, covariates] %>%
    dplyr::mutate_if(is.numeric, scale)
  data$trt_time <- ifelse(is.na(data$trt_time), 13, data$trt_time)

  time.max <- max(data$trt_time[!is.na(data$trt_time)])
  data <- transform(data, iid = as.numeric(id))

  data.cox <-
    data[which(data$time != 1 & data$time != time.max &
                 data$time <= data$trt_time),]
  data.cox$start <- data.cox$time - rep(1, length(data.cox$time))

  form <- survival::Surv(start, time, trt) ~ . - id - trt_time - iid
  model <- survival::coxph(form, data = data.cox)
  data.cox$p <- predict(model)
  nbp.distance <-
    matrix(rep(NA, length(unique(data$iid)) * length(unique(data$iid))),
           ncol = length(unique(data$iid)))
  nbp.t <-
    matrix(rep(NA, length(unique(data$iid)) * length(unique(data$iid))),
           ncol = length(unique(data$iid)))
  for (i in 1:length(unique(data$iid))) {
    for (j in 1:length(unique(data$iid))) {
      if (i == j) {
        nbp.distance[i, j] <- 999
      }
      if (i < j) {
        a <- unique(data[which(data$iid == i),]$trt_time)
        b <- unique(data[which(data$iid == j),]$trt_time)
        if (a == b) {
          nbp.distance[i, j] <- 999
        }
        else {
          trt_time <- min(a, b)
          nbp.t[i, j] <- trt_time
          if (length(data.cox[which(data.cox$iid == i &
                                    data.cox$time == trt_time - 1),]$p) != 0 &
              length(data.cox[which(data.cox$iid == j &
                                    data.cox$time == trt_time - 1),]$p) != 0 &
              length(data.cox[which(data.cox$iid == i &
                                    data.cox$time == trt_time),]$p) != 0 &
              length(data.cox[which(data.cox$iid == j &
                                    data.cox$time == trt_time),]$p) != 0) {
            nbp.distance[i, j] <- (data.cox[which(data.cox$iid == i &
                                                    data.cox$time == trt_time - 1),]$p -
                                     data.cox[which(data.cox$iid == j &
                                                      data.cox$time == trt_time - 1),]$p) ^ 2
          } else {
            nbp.distance[i, j] <- 999
          }
        }
      }
      if (i > j) {
        nbp.distance[i, j] <- nbp.distance[j, i]
      }
    }
  }
  nbp.distance1 <- nbp.distance * 1000
  nbp.distance2 <- nbpMatching::distancematrix(nbp.distance1)
  nbp <- nbpMatching::nonbimatch(nbp.distance2)
  nbpp <- nbp$halves[which(nbp$halves$Distance != 999000),]
  if (length(unique(data$id)) %% 2 != 0) {
    nbpp <- nbpp[-which(nbpp$Group1.ID == "ghost" | nbpp$Group2.ID == "ghost"),]
  }
  nbp.id1 <- nbpp$Group1.Row
  nbp.id2 <- nbpp$Group2.Row

  trt.id <- NULL
  con.id <- NULL

  for (i in 1:length(nbp.id1)) {
    if (data$trt_time[which(data$iid == nbp.id1[i])[1]] <
        data$trt_time[which(data$iid == nbp.id2[i])[1]]) {
      trt.id[i] = as.numeric(as.character(nbp.id1[i]))
      con.id[i] = as.numeric(as.character(nbp.id2[i]))
    } else {
      trt.id[i] = as.numeric(as.character(nbp.id2[i]))
      con.id[i] = as.numeric(as.character(nbp.id1[i]))
    }
  }

  for (i in 1:length(nbp.id1)) {
    trt.id[i] <-
      as.numeric(as.character(data$id[which(data$iid == trt.id[i])[1]]))
    con.id[i] <-
      as.numeric(as.character(data$id[which(data$iid == con.id[i])[1]]))
  }

  coxph_matches <- data.frame("trt.id" = trt.id, "con.id" = con.id,
                              "distance" = nbpp$Distance)

  return(coxph_matches)
}


