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
#' @return get a df of percentile values for days per unit
#' (P_20, P_50, P_60, P_70, P_80, P_85, P_90, Group)
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
    # flag used to be as.numeric
    e_pope[, `:=`(dif = as.integer(DateSupplied - t_nm1),
                  e_pop_est = as.integer(DateSupplied - t_nm1) / q_nm1)]
    
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
      # included more pop estimates for further testing
      P_60 = quantile(e_pope$e_pop_est, 0.60, type=2, na.rm = TRUE),
      P_70 = quantile(e_pope$e_pop_est, 0.70, type=2, na.rm = TRUE),
      P_80 = quantile(e_pope$e_pop_est, 0.80, type=2, na.rm = TRUE),
      P_85 = quantile(e_pope$e_pop_est, 0.85, type=2, na.rm = TRUE),
      P_90 = quantile(e_pope$e_pop_est, 0.90, type=2, na.rm = TRUE),
      Group = drug_code)
    
    # checks
    print("percentiles:")
    print(percentiles)
    return(percentiles)
  } 
  else {
    warning("no continuous episodes found for this drug code")
    return(data.frame(N = 0, P_20 = NA, P_50 = NA, P_60 = NA, P_70 = NA,
                      P_80 = NA, P_85 = NA, P_90 = NA, Group = drug_code))
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
#' @param combined_item_code3 population estimates output from e_pop_estimate
#' @param new_episode_threshold days threshold for defining new episodes (default 365)
#' @param recent_exposure_window days for recent exposure window (default 7)
#' @param keep_tmp_variable keep temporary variables in final df (default FALSE)
#' OPTIONAL PARAMETERS:
#' @param DeathDate if not supplied, will create it and assign NA (default NULL)
#' @param DateSupplied_index index supply date for study period, if not supplied,
#' will assign earliest supply date of chosen drug_code for each person id 
#' @param output_name name for the output df in global environment (default NULL)

#' @return get a df with exposure periods
#' 
exposure_by_drug <- function(drug_code, macro_d, EndDate, combined_item_code3, 
                             new_episode_threshold = 365, recent_exposure_window = 7, 
                             keep_tmp_variable = FALSE, DeathDate = NULL, 
                             DateSupplied_index = NULL, output_name = NULL) {
  
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
  
  # get population estimate
  e_est_macro <- combined_item_code3 %>%
    filter(Group == drug_code) %>%
    select(Group, P_80)
  
  # check that there exactly one value
  if(nrow(e_est_macro) != 1) {
    stop(paste0("there is more than 1 value for population quantile estimate for drug code ", drug_code))
  }
  
  # add to dispensing data directly
  P_80_value <- e_est_macro$P_80[1]
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
    results$t_nm1 <- as.Date(rep(NA, n_rows))
    results$t_nm2 <- as.Date(rep(NA, n_rows))
    results$t_nm3 <- as.Date(rep(NA, n_rows))
    results$q_nm1 <- numeric(n_rows)
    results$q_nm2 <- numeric(n_rows)
    results$q_nm3 <- numeric(n_rows)
    results$e_n <- numeric(n_rows)
    results$last_time <- as.Date(rep(NA, n_rows))
    results$last_q <- numeric(n_rows)
    results$first_date <- as.Date(patient_data$DateSupplied[1])
    results$unique_ep_id <- integer(n_rows)
    results$recent_exp <- recent_exposure_window
    results$no_formerly_exposed <- 0
    
    # init tracking variables
    ep_num <- 0
    episode_dispensing <- 0
    t_nm1 <- as.Date(NA)
    t_nm2 <- as.Date(NA)
    t_nm3 <- as.Date(NA)
    q_nm1 <- NA_real_
    q_nm2 <- NA_real_
    q_nm3 <- NA_real_
    e_n <- 99999999
    last_time <- as.Date(NA)
    last_q <- NA_real_
    unique_ep_id <- 0
    
    # row wise process
    for (i in 1:n_rows) {
      # curent values
      current_time <- as.Date(patient_data$DateSupplied[i])
      current_q <- patient_data$Quantity[i]
      
      # update history variables
      if (episode_dispensing >= 3) {
        t_nm3 <- t_nm2
        q_nm3 <- q_nm2
      }
      if (episode_dispensing >= 2) {
        t_nm2 <- t_nm1
        q_nm2 <- q_nm1
      }
      if (episode_dispensing >= 1) {
        t_nm1 <- last_time
        q_nm1 <- last_q
      }
      
      # check for new episode
      # flag floor, used to be as.numeric
      if (ep_num == 0 || (!is.na(t_nm1) && current_time > (t_nm1 + as.integer(e_n) + recent_exposure_window))) {
        # start new episode
        ep_num <- ep_num + 1
        unique_ep_id <- unique_ep_id + 1
        episode_dispensing <- 1
        t_nm1 <- as.Date(NA)
        t_nm2 <- as.Date(NA)
        t_nm3 <- as.Date(NA)
        q_nm1 <- NA_real_
        q_nm2 <- NA_real_
        q_nm3 <- NA_real_
        e_n <- current_q * patient_data$P_80[i]
      } else {
        # continue episode
        # check if formerly exposed / using ceiling here too
        # flag floor, used to be as.numeric
        if (!is.na(t_nm1) && current_time > (t_nm1 + as.integer(e_n) + recent_exposure_window)) {
          ep_num <- ep_num + 1
          results$no_formerly_exposed[i] <- 1
        }
        episode_dispensing <- episode_dispensing + 1
        
        # helper fct include safe division
        safe_divide <- function(numerator, denominator, default_value) {
          ifelse(denominator > 0, numerator / denominator, default_value)
        }
        
        # calculate exposure days using weighted formula, using the helper
        # here keeping as.numeric as don't want to truncate too early -- trying to match SAS intnx
        if (is.na(t_nm1)) {
          e_n <- current_q * patient_data$P_80[i]
        } else if (is.na(t_nm2)) {
          # flag floor, used as.integer
          term1 <- safe_divide(as.integer(current_time - t_nm1), q_nm1, patient_data$P_80[i])
          e_n <- current_q * ((3/6) * term1 + 
                                (2/6) * patient_data$P_80[i] + 
                                (1/6) * patient_data$P_80[i])
        } else if (is.na(t_nm3)) {
          term1 <- safe_divide(as.integer(current_time - t_nm1), q_nm1, patient_data$P_80[i]) 
          term2 <- safe_divide(as.integer(t_nm1 - t_nm2), q_nm2, patient_data$P_80[i])
          e_n <- current_q * ((3/6) * term1 + 
                                (2/6) * term2 + 
                                (1/6) * patient_data$P_80[i])
        } else {
          term1 <- safe_divide(as.integer(current_time - t_nm1), q_nm1, patient_data$P_80[i]) 
          term2 <- safe_divide(as.integer(t_nm1 - t_nm2), q_nm2, patient_data$P_80[i])
          term3 <- safe_divide(as.integer(t_nm2 - t_nm3), q_nm3, patient_data$P_80[i])
          e_n <- current_q * ((3/6) * term1 + 
                                (2/6) * term2 + 
                                (1/6) * term3)
        }
      }
      
      # store results for this row
      results$ep_num[i] <- ep_num
      results$episode_dispensing[i] <- episode_dispensing
      results$t_nm1[i] <- t_nm1
      results$t_nm2[i] <- t_nm2
      results$t_nm3[i] <- t_nm3
      results$q_nm1[i] <- q_nm1
      results$q_nm2[i] <- q_nm2
      results$q_nm3[i] <- q_nm3
      results$e_n[i] <- e_n
      results$unique_ep_id[i] <- unique_ep_id
      
      # update for next iteration
      last_time <- current_time
      last_q <- current_q
      results$last_time[i] <- last_time
      results$last_q[i] <- last_q
    }
    
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
  macro_d1[, row_num := NULL] # remove tmp column
  
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
    if (row$DateSupplied > death_end_date) {
      return(NULL) # want to return nothing if after deceased date 
    }
    
    # calculate endpoints
    # !! trying to align to SAS intnx, as.integer / floor truncate to lower
    current_end <- min(as.Date(row$DateSupplied) + as.integer(row$e_n), 
                       as.Date(row$SEE1) - 1,
                       death_end_date,
                       na.rm = TRUE)
    
    recent_end <- min(as.Date(row$DateSupplied) + as.integer(row$e_n + row$recent_exp), 
                      as.Date(row$SEE1) - 1,
                      death_end_date,
                      na.rm = TRUE)
    
    former_end <- min(as.Date(row$SEE1) - 1,
                      death_end_date,
                      na.rm = TRUE)
    
    # generate intervals
    intervals <- list()
    
    # es=1 -- current exposure interval
    interval1 <- row
    interval1$es <- 1L
    interval1$start_date <- as.Date(row$DateSupplied) - 1
    interval1$end_date <- as.Date(current_end)
    intervals[[1]] <- interval1 
    
    # es=2 -- recent exposure interval (if remaining time true)
    if (current_end < min(as.Date(row$SEE1) - 1, death_end_date, na.rm = TRUE)) {
      interval2 <- row
      interval2$es <- 2L
      interval2$start_date <- as.Date(current_end) # + 1 caused issues in recent exposure calculation
      interval2$end_date <- as.Date(recent_end)
      intervals[[2]] <- interval2
      
      # es=3 -- former exposure  
      if (recent_end < min(as.Date(row$SEE1) - 1, death_end_date, na.rm = TRUE)) {
        interval3 <- row
        interval3$es <- 3L
        interval3$start_date <- as.Date(recent_end) # + 1 caused issues in former exposure calculation
        interval3$end_date <- as.Date(former_end)
        intervals[[3]] <- interval3
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
  # flag floor, used to be as.integer
  macro_episodes[, pdays := as.integer(end_date - start_date)] # floor causing issues
  
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


## Column outputs detail if needed for debug or testing ----
#
# historical dispensing dates columns:
# t_nm1 date of the most recent previous dispensing in this episode
# t_nm2 date of the second most recent previous dispensing in this episode
# t_nm3 date of the third most recent previous dispensing in this episode
# last_time is most recent dispensing date
#
# historical dispensing quantities columns:
# q_nm1 quantity dispensed at the most recent previous dispensing
# q_nm2 quantity dispensed at the second most recent previous dispensing
# q_nm3 quantity dispensed at the third most recent previous dispensing
# last_q quantity from most recent dispensing
#
# tracking episodes columns: 
# episode_dispensing is number of dispensings within current episode
# unique_ep_id is unique identifier for each medication episode
# first_date first dispensing date for this patient
# ep_st_date start date of current episode
# no_formerly_exposed indicator if ppn was formerly exposed before this dispensing
# 
# exposure interval columns:
# es is exposure status (1=current, 2=recent, 3=former)
# start_date start date of this exposure interval
# end_date end date of this exposure interval
# date_of_supply dispensing date (should be NA for es=2 & es=3)
# pdays person-days for this interval
# start_t days from index date to interval start
# end_t days from index date to interval end
# 
# next dispensing columns:
# PPN1 technically is ppn of next record, but not useful and have not populated it 
# SEE1 date of next dispensing
# EP1 episode number of next dispensing
# last is indicator for last dispensing (with 1 == yes)
# 
# state transition tracking columns:
# lag_es is exposure status of previous interval
# rec_num is record number, where it increments as the exposure state changes
#

## Malcolm : issues / questions ----
# 1) quantile calculation in R vs SAS alignment (potentially fixed with type = 2)
# 2) SAS uses intnx for dates and I am not sure how to fully align to it in R -- may cause difference in final output
# 3) do we need to allow selection of pop quantile different to P80? SAS explicitly uses this one only
# 4) grace period handling? -- see below 


## Add grace period option? ----
# exposure_by_drug <- function(drug_code, macro_d, EndDate, combined_item_code3, 
#                              output_name, new_episode_threshold = 365, recent_exposure_window = 7, 
#                              keep_tmp_variable = FALSE,
#                              grace_period = 0) {
# 
#   grace period is not included in SAS code, at least did not spot it 
#   if using grace with IDP, prob makes more sense to add grace days to current exposure periods,
#   so would delay the start of the next episode
#   if I were to add grace days to the input data, this may cause issues in the weighed formula or
#   introduce assumptions about the average dispensing cycle -- or is this ok?
#   !! TODO ask Malcolm about this
#
#   # grace period could be introduced where I calculate patient exposure (ceiling here), i.e. when checking for new episode:
#   if (ep_num == 0 || (!is.na(t_nm1) && current_time > (t_nm1 + as.integer(e_n) + recent_exposure_window + grace_period))) {
#     # start new episode code
#   }
#   
#   # and also would need it for when I calculate exposure end:
#   current_end <- min(as.Date(row$DateSupplied) + as.integer(row$e_n) + # explicitly use ceiling
#                        grace_period, 
#                      as.Date(row$SEE1) - 1,
#                      death_end_date,
#                      na.rm = TRUE)
#
#   recent_end <- min(as.Date(row$DateSupplied) + as.integer(row$e_n + row$recent_exp) + # explicitly use ceiling
#                        grace_period,
#                      as.Date(row$SEE1) - 1,
#                      death_end_date,
#                      na.rm = TRUE)
# }
