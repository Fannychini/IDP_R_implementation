# ----------------------------------------- #
# R script for optimised IDP implementation #
# ----------------------------------------- #

## Library
library(tidyverse)
library(lubridate)
library(data.table)

## Reproducibility
set.seed(239)
sessionInfo()

# This is the optimised R implementation of IDP
# from https://github.com/c-bharat/IDP_exposure_method

# !! The implementation was tested on generated data only, validation is pending !!


## helper functions ----

### check the dates ----
#' validate date in a dataframe
#'
#' @param data df to validate
#' @param date_cols vector of column names that should be dates
#' @return input df
#'
validate_dates <- function(data, date_cols) {
  invalid_rows <- 0
  for (col in date_cols) {
    if (col %in% names(data)) {
      date_vals <- data[[col]]
      invalid <- !is.na(date_vals) & !inherits(date_vals, "Date")
      invalid_rows <- invalid_rows + sum(invalid)
    }
  }
  if (invalid_rows > 0) warning(paste0(invalid_rows, " rows with invalid date formats"))
  return(data)
}

### check required variables ----
#' validate required columns in a df
#'
#' @param data df to validate
#' @param required_cols vector of required column names
#' @return get TRUE if ok, otherwise stops with error
#'
validate_columns <- function(data, required_cols) {
  missing_cols <- required_cols[!required_cols %in% names(data)]
  if (length(missing_cols) > 0) {
    stop(paste("missing required columns:", paste(missing_cols, collapse=", ")))
  }
  return(TRUE)
}


## data preparation ----
# function to deal with necessary inputs and optional col


## e_pop_estimate ----
#' calculate population estimates for medication patterns
#'
#' @description
#' analyses dispensing data to estimate typical patterns in a population,
#' calculating percentiles of days/unit values
#'
#' @param drug_code the medication code to analyse
#' @param macro_d the input dispensing data with required columns:
#' PPN = unique person ids,
#' Group = col name of drug_code
#' Quantity = quantity of drug_code supplied,
#' DateSupplied = date of supply of drug_code
#'
#' @return get a df of percentile values for days per unit using type 2 quantiles
#' (P_20, P_50, P_60, P_70, P_80, P_85, P_90, Group)
#'
#' and extra values of P_80 testing for the different quantile types
#' (P_80_type1, P_80_type3, P_80_type4, P_80_type5, P_80_type6, P_80_type7, P_80_type8, P_80_type9)
#'
e_pop_estimate <- function(drug_code, macro_d) {
  # validate required columns with helper
  required_cols <- c("PPN", "DateSupplied", "Quantity", "Group")
  validate_columns(macro_d, required_cols)

  # validate dates with helper
  date_cols <- c("DateSupplied")
  macro_d <- validate_dates(macro_d, date_cols)

  # all dates are properly formatted
  macro_d <- macro_d %>%
    mutate(across(contains("date") | contains("Date"), as.Date))

  # filter and arrange relevant data
  filtered_data <- macro_d %>%
    filter(Group == drug_code) %>%
    arrange(PPN, DateSupplied)

  # convert to dt for faster processing
  dt <- as.data.table(filtered_data)

  # calculate intervals and identify new episodes
  dt[, t_nm1 := shift(DateSupplied), by = PPN]
  dt[, q_nm1 := shift(Quantity), by = PPN]
  # flag used to be as.integer
  dt[, new_episode := is.na(t_nm1) | as.integer(DateSupplied - t_nm1) >= 365]

  # generate e_pope (continuous episodes) and new_episode df
  e_pope <- dt[new_episode == FALSE]
  new_episode_data <- dt[new_episode == TRUE]

  # checks
  if (nrow(e_pope) > 0) {
    print("first 50 rows of e_pope:")
    print(head(e_pope, 50))
  }
  if (nrow(new_episode_data) > 0) {
    print("first 50 rows of new_episode:")
    print(head(new_episode_data, 50))
  }

  # calculate days per unit following a dispensing day
  if (nrow(e_pope) > 0) {
    # flag used to be as.numeric, updated to align best with SAS
    e_pope[, dif := as.numeric(DateSupplied - t_nm1)]
    e_pope[, e_pop_est := dif / q_nm1]

    # ! if zero or NA values
    e_pope[q_nm1 == 0 | is.na(q_nm1), e_pop_est := NA]

    e_pope[, first_interval := seq_len(.N) == 1, by = PPN]

    # checks
    print("first 50 rows of e_pop_est:")
    print(head(e_pope, 50))

    # calculate percentiles, using type 2 as this corresponds to SAS default as per
    # Wicklin, R. (2017) Sample quantiles: A comparison of 9 definitions; SAS Blog.
    # https://blogs.sas.com/content/iml/2017/05/24/definitions-sample-quantiles.html
    percentiles <- data.frame(
      N = nrow(e_pope),
      P_20 = quantile(e_pope$e_pop_est, 0.20, type=2, na.rm = TRUE),
      P_50 = quantile(e_pope$e_pop_est, 0.50, type=2, na.rm = TRUE),
      P_60 = quantile(e_pope$e_pop_est, 0.60, type=2, na.rm = TRUE),
      P_70 = quantile(e_pope$e_pop_est, 0.70, type=2, na.rm = TRUE),
      P_80 = quantile(e_pope$e_pop_est, 0.80, type=2, na.rm = TRUE),
      P_85 = quantile(e_pope$e_pop_est, 0.85, type=2, na.rm = TRUE),
      P_90 = quantile(e_pope$e_pop_est, 0.90, type=2, na.rm = TRUE),
      Group = drug_code,
      # included more pop estimates for further testing of P80
      P_80_type1 = quantile(e_pope$e_pop_est, 0.80, type=1, na.rm = TRUE),
      P_80_type3 = quantile(e_pope$e_pop_est, 0.80, type=3, na.rm = TRUE),
      P_80_type4 = quantile(e_pope$e_pop_est, 0.80, type=4, na.rm = TRUE),
      P_80_type5 = quantile(e_pope$e_pop_est, 0.80, type=5, na.rm = TRUE),
      P_80_type6 = quantile(e_pope$e_pop_est, 0.80, type=6, na.rm = TRUE),
      P_80_type7 = quantile(e_pope$e_pop_est, 0.80, type=7, na.rm = TRUE),
      P_80_type8 = quantile(e_pope$e_pop_est, 0.80, type=8, na.rm = TRUE),
      P_80_type9 = quantile(e_pope$e_pop_est, 0.80, type=9, na.rm = TRUE))

    # checks
    print("percentiles:")
    print(percentiles)
    return(percentiles)
  }
  else {
    warning("no continuous episodes found for this drug code")
    return(data.frame(N = 0, P_20 = NA, P_50 = NA, P_60 = NA, P_70 = NA,
                      P_80 = NA, P_85 = NA, P_90 = NA, Group = drug_code,
                      P_80_type1 = NA, P_80_type3 = NA, P_80_type4 = NA,
                      P_80_type5 = NA, P_80_type6 = NA, P_80_type7 = NA,
                      P_80_type8 = NA, P_80_type9 = NA))
  }
}


## exposure_by_drug ----
#' calculate medication exposure periods using IDP
#'
#' @description
#' process dispensing data to create exposure periods (current, recent, former)
#' for each PPN based on their dispensing pattern and population estimates.
#'
#' @param drug_code medication code
#' @param macro_d df with dispensing data
#' @param EndDate study end of follow-up/right censoring date
#' @param pop_estimates population estimates output from e_pop_estimate
#' @param new_episode_threshold days threshold for defining new episodes (default 365)
#' @param recent_exposure_window days for recent exposure window (default 7)
#' OPTIONAL PARAMETERS:
#' @param keep_tmp_variable keep temporary variables in final df (default FALSE)
#' @param pop_estimates if need another population estimates output from e_pop_estimate (default P_80)
#' @param DeathDate if not supplied, will create it and assign NA (default NULL)
#' @param DateSupplied_index index supply date for study period, if not supplied,
#' will assign earliest supply date of chosen drug_code for each person id
#' @param output_name name for the output df in global environment (default NULL)

#' @return get a df with exposure periods
#'
exposure_by_drug <- function(drug_code, macro_d, EndDate, pop_estimates,
                             new_episode_threshold = 365, recent_exposure_window = 7,
                             keep_tmp_variable = FALSE, percentile_col = "P_80",
                             DeathDate = NULL, DateSupplied_index = NULL, output_name = NULL) {

  # remove dt warnings
  datatable.verbose.orig <- getOption("datatable.verbose")
  options(datatable.verbose = FALSE)
  # on.exit to restore the original setting at the end
  on.exit(options(datatable.verbose = datatable.verbose.orig), add = TRUE)

  # validate required columns with helper
  required_cols <- c("PPN", "DateSupplied", "Quantity", "Group")
  validate_columns(macro_d, required_cols)

  # validate dates with helper
  date_cols <- c("DateSupplied")
  macro_d <- validate_dates(macro_d, date_cols)

  # handle DeathDate if not supplied
  if (is.null(DeathDate) || !(DeathDate %in% names(macro_d))) {
    # set col name if NULL
    DeathDate <- ifelse(is.null(DeathDate), "DeathDate", DeathDate)
    # set DeathDate with default value = NA
    macro_d[[DeathDate]] <- as.Date(rep(NA, nrow(macro_d)))
    message(sprintf("created '%s' column as Date type with NA values", DeathDate))
  }

  # handle DateSupplied_index if not supplied
  if (is.null(DateSupplied_index) || !(DateSupplied_index %in% names(macro_d))) {
    # set col name if NULL
    DateSupplied_index <- ifelse(is.null(DateSupplied_index), "DateSupplied_index", DateSupplied_index)
    # convert to dt for processing
    dt <- as.data.table(macro_d)
    # calculate min date by id using DateSupplied
    dt[, (DateSupplied_index) := min(DateSupplied, na.rm = TRUE), by = PPN]
    # back to df
    macro_d <- as.data.frame(dt)
    message(sprintf("created '%s' column with minimum date from DateSupplied for each id",
                    DateSupplied_index))
  }

  # re-validate dates with helper
  date_cols <- c("DateSupplied", "DeathDate", "DateSupplied_index")
  macro_d <- validate_dates(macro_d, date_cols)

  # all dates are properly formatted at the beginning
  macro_d <- macro_d %>%
    mutate(across(contains("date") | contains("Date"), as.Date))

  # get population estimate, updated to use another percentile if needed
  e_est_macro <- pop_estimates %>%
    filter(Group == drug_code) %>%
    select(Group, !!sym(percentile_col))

  # check that there exactly one value
  if(nrow(e_est_macro) != 1) {
    stop(paste0("there is more than 1 value for population quantile estimate for drug code ", drug_code))
  }

  # add to dispensing data directly, updated to use another percentile if needed
  P_value <- e_est_macro[[percentile_col]][1]
  macro_d_temp <- macro_d %>%
    filter(Group == drug_code & DateSupplied >= DateSupplied_index) %>%
    mutate(P_80 = P_80_value)

  ### calculate IDP exposure per PPN ----
  #' process data per PPN to calculate exposure periods
  #'
  #' @param patient_data data for a single patient
  #' @return processed data with exposure calculations
  calculate_patient_exposure <- function(patient_data) {
    # init necessary variables
    n_rows <- nrow(patient_data)

    # pre-allocate all results
    results <- patient_data
    results$ep_num <- integer(n_rows)
    results$episode_dispensing <- integer(n_rows)
    results$e_n <- numeric(n_rows)
    results$unique_ep_id <- integer(n_rows)
    results$recent_exp <- recent_exposure_window
    results$no_formerly_exposed <- 0
    # below for tracking t_ and q_ vars
    results$t_nm1 <- as.Date(rep(NA, n_rows))
    results$t_nm2 <- as.Date(rep(NA, n_rows))
    results$t_nm3 <- as.Date(rep(NA, n_rows))
    results$q_nm1 <- numeric(n_rows)
    results$q_nm2 <- numeric(n_rows)
    results$q_nm3 <- numeric(n_rows)
    # need to trackg for actual exposure days (pdays as in methods/SAS)
    results$actual_exposure_nm1 <- numeric(n_rows)
    results$actual_exposure_nm2 <- numeric(n_rows)
    results$actual_exposure_nm3 <- numeric(n_rows)

    # init tracking variables: track current episode numb / dispensing count within episode / unique episode id
    ep_num <- 0
    episode_dispensing <- 0
    unique_ep_id <- 0

    # row wise process
    for (i in 1:n_rows) {
      # curent values
      current_time <- as.Date(patient_data$DateSupplied[i])
      current_q <- patient_data$Quantity[i]

      # here want to determine if this is an index dispensing, i.e. as per methods
      # index dispensing = first dispensing OR first after former exposure gap
      is_index_dispensing <- FALSE

      if (i == 1) {
        is_index_dispensing <- TRUE
      } else {
        # check for episode break
        prev_time <- as.Date(patient_data$DateSupplied[i-1])
        prev_e_n <- results$e_n[i-1]
        gap_days <- as.numeric(current_time - prev_time)

        # episode breaks if gap > (e_n + recent_window)
        if (gap_days > (ceiling(prev_e_n) + recent_exposure_window)) {
          is_index_dispensing <- TRUE
        }
      }

      if (is_index_dispensing) {
        # new episode starts
        ep_num <- ep_num + 1
        unique_ep_id <- unique_ep_id + 1
        episode_dispensing <- 1

        # index dispensing using population estimate
        e_n <- current_q * patient_data$P_80[i]

        # adding printings for index dispensings
        if (i <= 2) {  # First two dispensings
          cat("\n check index dispensing", i-1, "calculatios \n")
          cat("date:", format(current_time, "%Y-%m-%d"), "\n")
          cat("new episode:", ep_num, "\n")
          cat("quantity:", current_q, "\n")
          cat("P80:", patient_data$P_80[i], "\n")
          cat("e_n:", e_n, "\n")
        }

        # flag if after former exposure
        if (ep_num > 1) {
          results$no_formerly_exposed[i] <- 1
        }

      } else {
        # then continue episode
        episode_dispensing <- episode_dispensing + 1

        # here now calculate e_n using actual exposures from previous dispensings
        # first find all previous dispensings in this episode
        episode_start_idx <- which(results$ep_num == ep_num)[1]
        episode_history_idx <- episode_start_idx:(i-1)

        ### flag: I misunderstood what the paper/SAS code used when describing
        ### 'actual exposure' = "define exposure as current until the next dispensing (dn + 1)"
        ### whicih in SAS is:
        ### `end_date = MIN(INTNX('DAY',date_of_supply,e_n,'END'), SEE1-1, DeathDate, &EndDate.);`
        ### BUT, in my previous code I thought if SAS caps exposure for the output intervals,
        ### it uses the same "actual" capped values when calculating the weighted average afterwards
        ### -> no, there are actually 2 types of dispensing intervals in the code
        ### 1) output exposure intervals where capped values are used (current exposure ends at next dispensing)
        ### 2) weighted average calculation where RAW !! intervals between dispensing dates are used
        ### SAS used `(Date_of_Supply-t_nm1)/q_nm1` but I used actual days of exposure instead

        ### -> final version of R implementation now use the raw intervals
        ### the use of the capped values cause issues when quantities varied because
        ### it was introducing different ratio in the weighted calculations

        # calculate actual exposures for each historical dispensing (e_n = 1)
        # using RAW intervals for weighted average and NOT capped actual exposures
        actual_exposures <- numeric(length(episode_history_idx))
        quantities <- numeric(length(episode_history_idx))

        for (j in seq_along(episode_history_idx)) {
          hist_idx <- episode_history_idx[j]

          # get dispensing data
          hist_date <- as.Date(patient_data$DateSupplied[hist_idx])
          hist_q <- patient_data$Quantity[hist_idx]

          # find when next dispensing occurred
          if (hist_idx < (i-1)) {
            next_date <- as.Date(patient_data$DateSupplied[hist_idx + 1])
          } else {
            # for the most recent historical dispensing
            next_date <- current_time
          }

          # need to use RAW interval between dispensings -- actually matching SAS
          actual_exposures[j] <- as.numeric(next_date - hist_date)
          quantities[j] <- hist_q
        }

        # get the most recent 3 values for weighted average
        n_hist <- length(actual_exposures)

        # calculate terms using actual exposures
        safe_divide <- function(num, den, default) {
          if (!is.na(den) && den > 0) num / den else default
        }

        ## WEIGHT CALCULATIONS ##
        # for most recent (weight 3/6)
        term1 <- if (n_hist >= 1) {
          safe_divide(actual_exposures[n_hist], quantities[n_hist], patient_data$P_80[i])
        } else {
          patient_data$P_80[i]
        }

        # second most recent (weight 2/6)
        term2 <- if (n_hist >= 2) {
          safe_divide(actual_exposures[n_hist-1], quantities[n_hist-1], patient_data$P_80[i])
        } else {
          patient_data$P_80[i]
        }

        # third most recent (weight 1/6)
        term3 <- if (n_hist >= 3) {
          safe_divide(actual_exposures[n_hist-2], quantities[n_hist-2], patient_data$P_80[i])
        } else {
          patient_data$P_80[i]
        }

        # now can calculate e_n
        e_n <- current_q * ((3/6) * term1 + (2/6) * term2 + (1/6) * term3)
      }

      # ## include printings -- something os weird when quantity varies on generated data
      # if (i <= 20 && !is_index_dispensing) {
      #   cat("\n check dispensing", i-1, "calculations \n")
      #   cat("episode history indices:", episode_history_idx, "\n")
      #   # I want each historical dispensing
      #   for (j in seq_along(episode_history_idx)) {
      #     idx <- episode_history_idx[j]
      #     cat(sprintf("history %d: date=%s, e_n=%.1f, actual_exp=%.0f, q=%d\n",
      #                 j,
      #                 format(patient_data$DateSupplied[idx], "%Y-%m-%d"),
      #                 results$e_n[idx],
      #                 actual_exposures[j],
      #                 quantities[j]))
      #   }
      #   # I want to see what I use for each calculation
      #   cat("\n am using for weighted average: \n")
      #   cat(sprintf("position %d (most recent): %.0f / %d = %.1f\n",
      #               n_hist, actual_exposures[n_hist], quantities[n_hist],
      #               actual_exposures[n_hist]/quantities[n_hist]))
      #   cat(sprintf("position %d: %.0f / %d = %.1f\n",
      #               n_hist-1, actual_exposures[n_hist-1], quantities[n_hist-1],
      #               actual_exposures[n_hist-1]/quantities[n_hist-1]))
      #   cat(sprintf("position %d: %.0f / %d = %.1f\n",
      #               n_hist-2, actual_exposures[n_hist-2], quantities[n_hist-2],
      #               actual_exposures[n_hist-2]/quantities[n_hist-2]))
      # }

      # store results
      results$ep_num[i] <- ep_num
      results$episode_dispensing[i] <- episode_dispensing
      results$e_n[i] <- e_n
      results$unique_ep_id[i] <- unique_ep_id
    }

    # add remaining tracking variables
    results$first_date <- as.Date(patient_data$DateSupplied[1])
    results$last_time <- results$DateSupplied
    results$last_q <- results$Quantity

    return(results)
  }

  ### process per PPN using helper function ----
  # convert to dt
  dt <- as.data.table(macro_d_temp)

  # get unique PPN
  patients <- unique(dt$PPN)

  # process each PPN and combine results
  processed_data <- lapply(patients, function(p) {
    patient_data <- dt[PPN == p][order(DateSupplied)]
    calculate_patient_exposure(patient_data)
  })

  # combine all results
  macro_d1 <- rbindlist(processed_data)

  # include death and censoring information
  # flag floor, used to be as.integer, now will use as.numeric
  macro_d1[, `:=`(death = as.integer(!is.na(DeathDate)),
                  cens_date = pmin(DeathDate, EndDate, na.rm = TRUE))]

  # look ahead to get date of next dispensing
  setorder(macro_d1, PPN, -DateSupplied)
  macro_d1[, `:=`(PPN1 = shift(PPN),
                  SEE1 = shift(DateSupplied),
                  EP1 = shift(ep_num)), by = PPN]

  # set values for first row (=last dispensing chronologically) of each PPN
  macro_d1[, row_num := seq_len(.N), by = PPN]
  macro_d1[row_num == 1, `:=`(SEE1 = as.Date("9999-12-31"),
                              last = 1L,
                              EP1 = NA_integer_)]
  # remove tmp column
  macro_d1[, row_num := NULL]

  # go back to chronological order
  setorder(macro_d1, PPN, DateSupplied)

  ### create exposure intervals ----
  #' create exposure intervals for a single dispensing record
  #'
  #' @param row a single dispensing record
  #' @return get a list of exposure intervals
  #'
  create_exposure_intervals <- function(row) {
    # handle if NAs
    death_end_date <- if(is.na(row$DeathDate)) EndDate else min(row$DeathDate, EndDate)

    # check if dispensing is after death/end date
    # return nothing if after deceased date
    if (row$DateSupplied > death_end_date) {
      return(NULL)
    }

    ## current exposure = from dispensing date - 1 in SAS
    current_start <- as.Date(row$DateSupplied) - 1

    # current exposure end is calculated from the dispensing date (as per methods/SAS)
    theoretical_end <- current_start + ceiling(row$e_n)
    next_limit <- if (!is.na(row$SEE1) && row$SEE1 < as.Date("9999-01-01")) {
      as.Date(row$SEE1) - 1
    } else {
      death_end_date
    }

    current_end <- min(theoretical_end, next_limit, death_end_date)

    # generate intervals
    intervals <- list()

    # always generate current exposure interval
    interval1 <- row
    interval1$es <- 1L
    interval1$start_date <- current_start
    interval1$end_date <- current_end
    intervals[[1]] <- interval1

    ## recent exposure = starts day after current end + 1 in SAS
    if (current_end < min(next_limit, death_end_date)) {
      recent_start <- current_end + 1

      # end of recent exposure:
      recent_theoretical_end <- theoretical_end + row$recent_exp
      recent_end <- min(recent_theoretical_end, next_limit, death_end_date)

      if (recent_start <= recent_end) {
        interval2 <- row
        interval2$es <- 2L
        interval2$start_date <- recent_start
        interval2$end_date <- recent_end
        intervals[[2]] <- interval2

        ## former exposure = from end of recent until next dispensing in metods
        if (recent_end < min(next_limit, death_end_date)) {
          former_start <- recent_end + 1
          former_end <- min(next_limit, death_end_date)

          if (former_start <= former_end) {
            interval3 <- row
            interval3$es <- 3L
            interval3$start_date <- former_start
            interval3$end_date <- former_end
            intervals[[3]] <- interval3
          }
        }
      }
    }

    return(intervals)
  }

  # get interval, row wise
  all_intervals <- lapply(1:nrow(macro_d1), function(i) {
    create_exposure_intervals(macro_d1[i])
  })

  # unlist and combine
  all_intervals_flat <- unlist(all_intervals, recursive = FALSE)
  all_intervals_flat <- all_intervals_flat[!sapply(all_intervals_flat, is.null)]

  macro_episodes <- rbindlist(all_intervals_flat)

  # resort episodes
  setorder(macro_episodes, PPN, start_date)

  ### final processing ----
  # clean up variables for non-active periods
  macro_episodes[es %in% c(2L, 3L), `:=`(date_of_supply = as.Date(NA),
                                         Quantity = NA_real_,
                                         e_n = NA_real_,
                                         episode_dispensing = NA_integer_)]

  # calculate start_t and end_t
  # used to be as.integer, flag as changed to floor
  macro_episodes[, `:=`(start_t = as.integer(start_date - (DateSupplied_index - 1)),
                        end_t = as.integer(end_date - (DateSupplied_index - 1)))]

  # check dates again because there is yet another date issue
  date_cols <- c("t_nm1", "t_nm2", "t_nm3", "last_time", "first_date")
  for (col in date_cols) {
    if (col %in% names(macro_episodes)) {
      macro_episodes[[col]] <- as.Date(macro_episodes[[col]])
    }
  }

  # refresh the dt reference to avoid the warning  Invalid .internal.selfref detected and
  # fixed by taking a (shallow) copy of the data.table so that := can add this new column by reference
  setDT(macro_episodes)
  # get episode start date
  macro_episodes[, ep_st_date := first(DateSupplied), by = .(PPN, ep_num)]

  # generate person-days
  # flag floor, used to be as.integer -- aliging with SAS
  macro_episodes[, pdays := as.integer(end_date - start_date) + 1]

  # get lag_es for tracking state changes
  macro_episodes[, lag_es := shift(es), by = PPN]

  # create custom drug-specific columns
  ep_num_col <- paste0("ep_num", drug_code)
  first_date_col <- paste0("first_", drug_code, "_date")

  # init columns
  macro_episodes[[ep_num_col]] <- NA_integer_
  macro_episodes[[first_date_col]] <- as.Date(NA)

  ### optimise episode numbering ----
  #' process episode numbering for a single PPN
  #'
  #' @param patient_data data for a single PPN
  #' @return get processed data with episode numbering
  #'
  process_patient_episodes <- function(patient_data) {
    n_rows <- nrow(patient_data)

    # init values for first row
    patient_data[[ep_num_col]][1] <- 1L
    patient_data$rec_num[1] <- 1L
    patient_data[[first_date_col]][1] <- as.Date(patient_data$DateSupplied[1])

    # process remaining rows
    if (n_rows > 1) {
      for (i in 2:n_rows) {
        # remember first date
        patient_data[[first_date_col]][i] <- patient_data[[first_date_col]][1]

        # set record number based on exposure state changes
        patient_data$rec_num[i] <- ifelse(
          patient_data$es[i] != patient_data$es[i-1],
          patient_data$rec_num[i-1] + 1L,
          patient_data$rec_num[i-1])

        # set episode number, increment if es=1 follows es=3
        patient_data[[ep_num_col]][i] <- ifelse(
          patient_data$es[i] == 1 && patient_data$lag_es[i] == 3,
          patient_data[[ep_num_col]][i-1] + 1L,
          patient_data[[ep_num_col]][i-1])
      }
    }

    return(patient_data)
  }

  # process per PPN
  patients <- unique(macro_episodes$PPN)
  processed_episodes <- lapply(patients, function(p) {
    patient_data <- macro_episodes[PPN == p][order(start_date)]
    process_patient_episodes(patient_data)
  })

  # combine results
  output_data <- rbindlist(processed_episodes)

  # remove ep_num to avoid confusion
  output_data[, ep_num := NULL]

  # convert to dt for consistency
  output_data <- as.data.frame(output_data)

  # get rid of intermediate columns, should set = TRUE if needed for dvp/testing
  if (!keep_tmp_variable) {
    # tmp variables to drop -- see end of script for what they mean
    intermediate_cols <- c("t_nm1", "t_nm2", "t_nm3", "q_nm1", "q_nm2", "q_nm3",
                           "last_time", "last_q", "no_formerly_exposed",
                           "lag_es", "PPN1")
    output_data <- output_data[, !names(output_data) %in% intermediate_cols]
  }

  # assign to global environment if output_name is provided
  if (!is.null(output_name)) {
    assign(output_name, output_data, envir = .GlobalEnv)
  }

  return(output_data)
}


## Column outputs detail ----
#
# -- Intermediate variables (removed unless keep_tmp_variable = TRUE) --
# historical dispensing columns used in weighted average calculation:
# t_nm1, t_nm2, t_nm3 - dates of previous dispensings
# q_nm1, q_nm2, q_nm3 - quantities from previous dispensings
# last_time - copy of DateSupplied for current row
# last_q - copy of Quantity for current row
# no_formerly_exposed - indicator (1) if this dispensing starts a new episode after former exposure
# lag_es - exposure status (es) from previous row, used to detect state transitions
# PPN1 - PPN from next row, only used in look-ahead
#
# -- Main output variables --
# 1) core dispensing columns:
# PPN - unique person identifier
# DateSupplied - date of dispensing
# Quantity - quantity dispensed on this date
# Group - medication/drug code
# P_80 - population-level 80th percentile estimate of days per unit
# DateSupplied_index - first dispensing date in study period for this ppn
# DeathDate - date of death (if applicable)
#
# 2) episode tracking columns:
# ep_num{drug_code} - episode number for this specific drug (increments after former exposure)
# episode_dispensing - dispensing count within current episode (1 for index, >=2 for subsequent)
# unique_ep_id - globally unique identifier for each medication episode
# e_n - estimated exposure days for this dispensing (quantity × weighted avg days/unit)
# first_{drug_code}_date - first dispensing date for this drug for this person
# ep_st_date - start date of current episode (DateSupplied of first dispensing in episode)
#
# 3) exposure interval columns:
# es - exposure status (1=current, 2=recent, 3=former)
# start_date - start date of this exposure interval (DateSupplied - 1 for current exposure)
# end_date - end date of this exposure interval
# pdays - person-days in this interval (end_date - start_date + 1)
# start_t - days from DateSupplied_index to start_date
# end_t - days from DateSupplied_index to end_date
#
# 4) next dispensing columns (look-ahead):
# SEE1 - date of next dispensing (set as 9999-12-31 for last dispensing)
# EP1 - episode number of next dispensing (set as NA for last dispensing)
# last - indicator for last dispensing in dataset (1 = yes)
#
# 5) other derived columns:
# recent_exp - length of recent exposure window (default 7 days)
# death - death indicator (1 if DeathDate exists, 0 otherwise)
# cens_date - censoring date (minimum of DeathDate and study EndDate)
# rec_num - row counter that increments when exposure status (es) changes
#
# Noting that this code aligns with SAS behaviour and includes NA values for:
# es=2 (recent) and es=3 (former) intervals, DateSupplied,
# Quantity, e_n, and episode_dispensing, because they represent gaps between dispensings


## Malcolm : issues / questions ----
# *fixed*
# 1) quantile calculation in R vs SAS alignment (type = 2)
# 2) SAS uses intnx for dates, think current implementation aligns with it

# *potential stuff to change*
# 3) do we need to allow selection of pop quantile different to P80? SAS explicitly uses this one only
# 4) grace period handling?


## Add grace period option? ----
#
#   grace period is not included in SAS code, at least did not spot it
#   if using grace with IDP, prob makes more sense to add grace days to current exposure periods,
#   so would delay the start of the next episode
#   if I were to add grace days to the input data, this may cause issues in the weighed formula or
#   introduce assumptions about the average dispensing cycle -- or is this ok?
#
# exposure_by_drug <- function(drug_code, macro_d, EndDate, pop_estimates,
#                              new_episode_threshold = 365, recent_exposure_window = 7,
#                              keep_tmp_variable = FALSE, percentile_col = "P_80",
#                              DeathDate = NULL, DateSupplied_index = NULL, output_name = NULL,
#                              grace_period = 0) {
#
#   # grace period could be introduced where I calculate patient exposure, i.e. when checking for new episode:
#   calculate_patient_exposure <- function(patient_data, grace_period) {
#
#     processed_data <- lapply(patients, function(p) {
#       patient_data <- dt[PPN == p][order(DateSupplied)]
#       calculate_patient_exposure(patient_data, grace_period)
#     })
#
#     results$recent_exp <- recent_exposure_window
#     results$grace_period <- grace_period
#
#     if (gap_days > (ceiling(prev_e_n) + recent_exposure_window + grace_period)) {
#       is_index_dispensing <- TRUE
#     }
#   }
#
#   # and also in create_exposure_intervals:
#   theoretical_end <- current_start + ceiling(row$e_n) + row$grace_period
#
#   recent_theoretical_end <- theoretical_end + row$recent_exp
#
# }
