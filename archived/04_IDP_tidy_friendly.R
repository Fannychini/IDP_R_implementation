# ------------------------------------- #
# Tidyverse-friendly IDP implementation #
# ------------------------------------- #

## Libraries
library(tidyverse)
library(lubridate)

## Set seed for reproducibility
set.seed(239)

## this is not completed and is more complicated to implement than I thought :)


#' check required columns in a dataframe
#' 
#' @param data df to validate
#' @param required_cols vector of required column names
#' @return get TRUE if validation passes, otherwise stops with error
#' 
validate_columns <- function(data, required_cols) {
  missing_cols <- required_cols[!required_cols %in% names(data)]
  if (length(missing_cols) > 0) {
    stop(paste("missing required columns:", paste(missing_cols, collapse=", ")))
  }
  return(TRUE)
}

#' check date fields in a df
#' 
#' @param data df to check
#' @param date_cols vector of column names that should be dates
#' @return get input df
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


#' calculate population estimates for medication patterns
#' 
#' @param item_code the medication code to analyse
#' @param macro_d the input dispensing data with required columns: 
#' PPN, Date_of_Supply, q_D, group
#' @return get a df of percentile values for days per unit (P_20, P_50, P_80, P_90)
#' 
e_pop_estimate_tidy <- function(item_code, macro_d) {
  # check required columns
  required_cols <- c("PPN", "Date_of_Supply", "q_D", "group")
  validate_columns(macro_d, required_cols)
  
  # check dates
  date_cols <- c("Date_of_Supply", "DeathDate", "Date_of_Supply_index")
  macro_d <- validate_dates(macro_d, date_cols)
  
  # all dates are properly formatted
  macro_d <- macro_d %>%
    mutate(across(contains("date") | contains("Date"), as.Date))
  
  # filter and arrange relevant data
  filtered_data <- macro_d %>%
    filter(group == item_code) %>%
    arrange(PPN, Date_of_Supply)
  
  # Calculate intervals between dispensings and identify new episodes
  continuous_episodes <- filtered_data %>%
    group_by(PPN) %>%
    mutate(prev_date = lag(Date_of_Supply),
           prev_quantity = lag(q_D),
           days_since_prev = as.numeric(Date_of_Supply - prev_date),
           is_new_episode = is.na(prev_date) | days_since_prev >= 365) %>%
    filter(!is_new_episode) %>%
    # Calculate days per unit with division by zero protection
    mutate(e_pop_est = case_when(
      prev_quantity > 0 ~ days_since_prev / prev_quantity,
      TRUE ~ NA_real_),
      first_interval = row_number() == 1) %>%
    ungroup()
  
  # Output diagnostic information
  if (nrow(continuous_episodes) > 0) {
    message("First 50 rows of continuous episodes:")
    print(head(continuous_episodes, 50))
  } else {
    warning("No continuous episodes found for this item_code")
    return(data.frame(N = 0, P_20 = NA, P_50 = NA,
                      P_80 = NA, P_90 = NA, item_code = item_code))
  }
  
  # Calculate percentiles
  percentiles <- continuous_episodes %>%
    summarise(N = n(),
              P_20 = quantile(e_pope$e_pop_est, 0.20, type=2, na.rm = TRUE),
              P_50 = quantile(e_pope$e_pop_est, 0.50, type=2, na.rm = TRUE),
              # included more pop estimates for further testing
              P_60 = quantile(e_pope$e_pop_est, 0.60, type=2, na.rm = TRUE),
              P_70 = quantile(e_pope$e_pop_est, 0.70, type=2, na.rm = TRUE),
              P_80 = quantile(e_pope$e_pop_est, 0.80, type=2, na.rm = TRUE),
              P_85 = quantile(e_pope$e_pop_est, 0.85, type=2, na.rm = TRUE),
              P_90 = quantile(e_pope$e_pop_est, 0.90, type=2, na.rm = TRUE)) %>%
    mutate(item_code = item_code)
  
  message("percentiles:")
  print(percentiles)
  
  return(percentiles)
}

#' Calculate exposure periods for a medication
#' 
#' @param drug_number medication code
#' @param macro_d df with dispensing data
#' @param EndDate study end date
#' @param combined_item_code3 population estimates from e_pop_estimate
#' @param output_name name for the output df in global environment (default NULL)
#' @param new_episode_threshold days threshold for defining new episodes (default 365)
#' @param recent_exposure_window days for recent exposure window (default 7)
#' @param keep_tmp_variable keep temporary variables in final df (default FALSE)
#' @return get a df with exposure periods
#' 
exposure_by_drug_tidy <- function(drug_number, macro_d, EndDate, combined_item_code3, 
                                  output_name = NULL, new_episode_threshold = 365, 
                                  recent_exposure_window = 7, keep_tmp_variable = FALSE) {
  # Validate required columns
  required_cols <- c("PPN", "Date_of_Supply", "q_D", "group", "Date_of_Supply_index")
  validate_columns(macro_d, required_cols)
  
  # Validate dates
  date_cols <- c("Date_of_Supply", "DeathDate", "Date_of_Supply_index")
  macro_d <- validate_dates(macro_d, date_cols)
  
  # Ensure dates are properly formatted
  macro_d <- macro_d %>%
    mutate(across(contains("date") | contains("Date"), as.Date))
  
  # Validate EndDate
  EndDate <- as.Date(EndDate)
  
  # Filter dispensing data for the specific drug
  macro_d_temp <- macro_d %>%
    filter(group == drug_number & Date_of_Supply >= Date_of_Supply_index)
  
  # Get population estimate
  e_est_macro <- combined_item_code3 %>%
    filter(group == drug_number) %>%
    select(group, P_80)
  
  # IMPROVED: Check if we have exactly one value
  if(nrow(e_est_macro) != 1) {
    stop(paste0("Expected exactly one row in population estimates for drug number ", drug_number))
  }
  
  # Add population estimate to all rows
  P_80_value <- e_est_macro$P_80[1]
  macro_d_temp <- macro_d_temp %>%
    mutate(P_80 = P_80_value)
  
  # Process each patient to calculate exposure
  # This part is complex and requires tracking state, so we'll use nest/unnest
  processed_data <- macro_d_temp %>%
    group_by(PPN) %>%
    arrange(Date_of_Supply) %>%
    # Step 1: Pre-process each patient to create needed variables
    mutate(patient_index = row_number(),
           first_date = first(Date_of_Supply)) %>%
    # Step 2: Calculate exposure periods
    nest() %>%
    # Apply a function to process each patient's data
    mutate(processed = map(data, ~ process_patient_exposure(., 
                                                            recent_exp = recent_exposure_window,
                                                            new_episode_threshold = new_episode_threshold))) %>%
    unnest(processed) %>%
    ungroup()
  
  # Add death and censoring information
  processed_data <- processed_data %>%
    mutate(death = !is.na(DeathDate),
           cens_date = pmin(DeathDate, EndDate, na.rm = TRUE))
  
  # Calculate next dispensing date for each patient
  processed_with_next <- processed_data %>%
    group_by(PPN) %>%
    arrange(PPN, desc(Date_of_Supply)) %>%
    mutate(SEE1 = lead(Date_of_Supply), 
           EP1 = lead(ep_num)) %>%
    # Set values for last dispensing chronologically
    mutate(SEE1 = if_else(row_number() == 1, as.Date("9999-12-31"), SEE1),
           last = if_else(row_number() == 1, 1L, 0L),
           EP1 = if_else(row_number() == 1, NA_integer_, EP1)) %>%
    ungroup() %>%
    # Go back to chronological order
    arrange(PPN, Date_of_Supply)
  
  # Create exposure intervals (current, recent, former)
  exposure_intervals <- processed_with_next %>%
    # Use pmap to generate all intervals for each row
    mutate(intervals = pmap(list(Date_of_Supply, e_n, SEE1, DeathDate, recent_exp), 
                            ~ create_exposure_intervals(..1, ..2, ..3, ..4, EndDate, ..5))) %>%
    unnest(intervals) %>%
    # Clean up variables for non-active periods
    mutate(Date_of_Supply = if_else(es %in% c(2L, 3L), as.Date(NA), Date_of_Supply),
           q_D = if_else(es %in% c(2L, 3L), NA_real_, q_D),
           e_n = if_else(es %in% c(2L, 3L), NA_real_, e_n),
           episode_dispensing = if_else(es %in% c(2L, 3L), NA_integer_, episode_dispensing)) %>%
    # Calculate start_t and end_t
    mutate(start_t = as.numeric(start_date - (Date_of_Supply_index - 1)),
           end_t = as.numeric(end_date - (Date_of_Supply_index - 1))) %>%
    # Get episode start date
    group_by(PPN, ep_num) %>%
    mutate(ep_st_date = first(Date_of_Supply)) %>%
    ungroup() %>%
    # Generate person-days
    mutate(pdays = as.numeric(end_date - start_date)) %>%
    # Get lag_es for tracking state changes
    group_by(PPN) %>%
    arrange(PPN, start_date) %>%
    mutate(lag_es = lag(es)) %>%
    ungroup()
  
  # Create custom drug-specific columns
  ep_num_col <- paste0("ep_num", drug_number)
  first_date_col <- paste0("first_", drug_number, "_date")
  
  # Process episode numbering
  output_data <- exposure_intervals %>%
    group_by(PPN) %>%
    arrange(PPN, start_date) %>%
    mutate(rec_num = cumsum(es != lag(es, default = -1)),
           !!ep_num_col := cumsum(es == 1 & lag(es, default = 0) == 3) + 1L,
           !!first_date_col := first(Date_of_Supply)) %>%
    ungroup() %>%
    select(-ep_num) # Remove original ep_num to avoid confusion
  
  # Debugging info
  message("sample of final output:")
  print(head(output_data, 10))
  
  # Assign to global environment with provided output name
  assign(output_name, output_data, envir = .GlobalEnv)
  
  return(output_data)
}

#' Helper function to process a patient's exposure data
#' 
#' @param patient_data Dataframe with one patient's data
#' @param recent_exp Days for recent exposure
#' @param new_episode_threshold Days threshold for new episodes
#' @return Processed dataframe with exposure calculations
#' 
process_patient_exposure <- function(patient_data, recent_exp = 7, new_episode_threshold = 365) {
  # Initialise result dataframe
  results <- patient_data %>%
    mutate(ep_num = 0L,
           episode_dispensing = 0L,
           t_nm1 = as.Date(NA),
           t_nm2 = as.Date(NA),
           t_nm3 = as.Date(NA),
           q_nm1 = NA_real_,
           q_nm2 = NA_real_,
           q_nm3 = NA_real_,
           e_n = NA_real_,
           last_time = as.Date(NA),
           last_q = NA_real_,
           unique_ep_id = 0L,
           recent_exp = recent_exp,
           no_formerly_exposed = 0L)
  
  # If no rows, return empty dataframe
  if (nrow(results) == 0) return(results)
  
  # Process row by row (unfortunately this is hard to vectorize)
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
  
  for (i in 1:nrow(results)) {
    # Current values
    current_time <- results$Date_of_Supply[i]
    current_q <- results$q_D[i]
    
    # Update history variables
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
    
    # Check for new episode - using parameter
    if (ep_num == 0 || (!is.na(t_nm1) && current_time > (t_nm1 + e_n + recent_exp))) {
      # Start new episode
      ep_num <- ep_num + 1
      unique_ep_id <- unique_ep_id + 1
      episode_dispensing <- 1
      t_nm1 <- as.Date(NA)
      t_nm2 <- as.Date(NA)
      t_nm3 <- as.Date(NA)
      q_nm1 <- NA_real_
      q_nm2 <- NA_real_
      q_nm3 <- NA_real_
      e_n <- current_q * results$P_80[i]
    } else {
      # Continue episode
      # Check if formerly exposed
      if (!is.na(t_nm1) && current_time > (t_nm1 + e_n + recent_exp)) {
        ep_num <- ep_num + 1
        results$no_formerly_exposed[i] <- 1
      }
      episode_dispensing <- episode_dispensing + 1
      
      # Calculate exposure days using weighted formula with division by zero protection
      if (is.na(t_nm1)) {
        e_n <- current_q * results$P_80[i]
      } else if (is.na(t_nm2)) {
        if (q_nm1 > 0) {
          e_n <- current_q * ((3/6) * as.numeric(current_time - t_nm1) / q_nm1 + 
                                (2/6) * results$P_80[i] + 
                                (1/6) * results$P_80[i])
        } else {
          e_n <- current_q * results$P_80[i]  # Fallback if q_nm1 is zero
        }
      } else if (is.na(t_nm3)) {
        if (q_nm1 > 0 && q_nm2 > 0) {
          e_n <- current_q * ((3/6) * as.numeric(current_time - t_nm1) / q_nm1 + 
                                (2/6) * as.numeric(t_nm1 - t_nm2) / q_nm2 + 
                                (1/6) * results$P_80[i])
        } else if (q_nm1 > 0) {
          e_n <- current_q * ((3/6) * as.numeric(current_time - t_nm1) / q_nm1 + 
                                (3/6) * results$P_80[i])
        } else {
          e_n <- current_q * results$P_80[i]  # Fallback
        }
      } else {
        if (q_nm1 > 0 && q_nm2 > 0 && q_nm3 > 0) {
          e_n <- current_q * ((3/6) * as.numeric(current_time - t_nm1) / q_nm1 + 
                                (2/6) * as.numeric(t_nm1 - t_nm2) / q_nm2 + 
                                (1/6) * as.numeric(t_nm2 - t_nm3) / q_nm3)
        } else if (q_nm1 > 0 && q_nm2 > 0) {
          e_n <- current_q * ((3/6) * as.numeric(current_time - t_nm1) / q_nm1 + 
                                (3/6) * as.numeric(t_nm1 - t_nm2) / q_nm2)
        } else if (q_nm1 > 0) {
          e_n <- current_q * ((3/6) * as.numeric(current_time - t_nm1) / q_nm1 + 
                                (3/6) * results$P_80[i])
        } else {
          e_n <- current_q * results$P_80[i]  # Fallback
        }
      }
    }
    
    # Store results for this row
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
    
    # Update for next iteration
    last_time <- current_time
    last_q <- current_q
    results$last_time[i] <- last_time
    results$last_q[i] <- last_q
  }
  
  return(results)
}

#' Create exposure intervals from a dispensing row
#'
#' @param dispense_date Dispensing date
#' @param exposure_days Calculated exposure days
#' @param next_date Next dispensing date
#' @param death_date Death date (or NA)
#' @param end_date Study end date
#' @param recent_days Days for recent exposure window
#' @return A tibble with exposure intervals
create_exposure_intervals <- function(dispense_date, exposure_days, next_date, 
                                      death_date, end_date, recent_days = 7) {
  # Handle NA inputs
  next_date <- if(is.na(next_date)) as.Date("9999-12-31") else next_date
  
  # Improved NA handling
  death_end_date <- if(is.na(death_date)) end_date else min(death_date, end_date)
  
  # Check if dispensing is after death/end date
  if (dispense_date > death_end_date) {
    return(tibble()) # Return empty tibble if after deceased date
  }
  
  # Calculate end points
  current_end <- min(dispense_date + exposure_days,
                     next_date - 1,
                     death_end_date,
                     na.rm = TRUE)
  
  recent_end <- min(dispense_date + exposure_days + recent_days,
                    next_date - 1,
                    death_end_date,
                    na.rm = TRUE)
  
  former_end <- min(next_date - 1,
                    death_end_date, 
                    na.rm = TRUE)
  
  # Generate intervals - start with an empty list
  intervals <- tibble()
  
  # es=1 -- current exposure interval
  intervals <- bind_rows(intervals,
                         tibble(es = 1L,
                                start_date = dispense_date,
                                end_date = current_end))
  
  # es=2 -- recent exposure interval (if remaining time)
  if (current_end < min(next_date - 1, death_end_date, na.rm = TRUE)) {
    intervals <- bind_rows(intervals,
                           tibble(es = 2L,
                                  start_date = current_end + 1,
                                  end_date = recent_end))
    
    # es=3 -- former exposure
    if (recent_end < min(next_date - 1, death_end_date, na.rm = TRUE)) {
      intervals <- bind_rows(intervals,
                             tibble(es = 3L,
                                    start_date = recent_end + 1,
                                    end_date = former_end))
    }
  }
  
  return(intervals)
}

# Test function to create sample data
create_test_data <- function() {
  tibble(
    PPN = rep(c(1001, 1002, 1003, 1004), times = c(5, 3, 6, 4)),
    group = rep(c(101, 101, 101, 102), times = c(5, 3, 6, 4)),
    Date_of_Supply = as.Date(c("2020-01-15", "2020-02-15", "2020-03-15", "2020-04-15", "2020-05-15",
                               "2020-01-10", "2020-02-10", "2021-03-20",
                               "2020-01-05", "2020-02-10", "2020-04-20", "2020-07-01", "2020-08-15", "2021-09-20",
                               "2020-01-15", "2020-03-15", "2020-06-15", "2020-09-15")),
    q_D = c(30, 30, 30, 30, 30, # PPN 1001
            28, 28, 28, # PPN 1002
            14, 28, 30, 60, 30, 30, # PPN 1003
            60, 60, 60, 60), # PPN 1004
    Date_of_Supply_index = as.Date(c(rep("2020-01-01", 5), rep("2020-01-01", 3),
                                     rep("2020-01-01", 6), rep("2020-01-01", 4))),
    DeathDate = as.Date(c(rep(NA, 5), rep(NA, 2), "2021-06-30", # PPN 1002 deceased status
                          rep(NA, 6), rep(NA, 4))))
}

# Full test of tidyverse implementation
test_tidyverse_implementation <- function() {
  message("\nTesting tidyverse IDP implementation\n")
  
  # Create test data
  test_data <- create_test_data()
  
  # Calculate population estimates
  message("\nCalculating population estimates\n")
  percentiles_101 <- e_pop_estimate_tidy(101, test_data)
  percentiles_102 <- e_pop_estimate_tidy(102, test_data)
  
  # Combine percentiles
  combined_item_code <- bind_rows(percentiles_101, percentiles_102)
  combined_item_code3 <- combined_item_code %>%
    mutate(group = as.numeric(item_code))
  
  message("\nCombined percentiles:\n")
  print(combined_item_code3)
  
  # Define study end date
  EndDate <- as.Date("2021-12-31")
  
  # Calculate exposures
  message("\nCalculating exposures for group 101\n")
  exposure_101 <- exposure_by_drug_tidy(101, test_data, EndDate, combined_item_code3, "output_101_tidy")
  
  message("\nCalculating exposures for group 102\n")
  exposure_102 <- exposure_by_drug_tidy(102, test_data, EndDate, combined_item_code3, "output_102_tidy")
  
  # Verify results
  message("\nVerifying results\n")
  
  # Check exposure status
  message("\nExposure status for group 101\n")
  print(table(exposure_101$es))
  
  message("\nExposure status for group 102\n")
  print(table(exposure_102$es))
  
  # Check patient counts
  message("\nPPN count for group 101:", length(unique(exposure_101$PPN)))
  message("\nPPN count for group 102:", length(unique(exposure_102$PPN)))
  
  # Check specific patient
  message("\nExposure periods for PPN 1001 (group 101)")
  print(exposure_101 %>% 
          filter(PPN == 1001) %>% 
          select(PPN, start_date, end_date, es))
  
  # Check deceased patient
  message("\nExposure for PPN 1002 with death set as 2021-06-30")
  print(exposure_101 %>% 
          filter(PPN == 1002) %>% 
          select(PPN, start_date, end_date, es, DeathDate))
  
  # Check episode numbers
  ep_num_col <- paste0("ep_num", 101)
  message("\nEpisode counts for group 101")
  print(table(exposure_101[[ep_num_col]]))
  
  # Average exposure durations
  message("\nAverage duration of current exposure (es=1) for group 101:", 
          mean(exposure_101$pdays[exposure_101$es == 1]), "days")
  message("\nAverage duration of recent exposure (es=2) for group 101:", 
          mean(exposure_101$pdays[exposure_101$es == 2]), "days")
  message("\nAverage duration of former exposure (es=3) for group 101:", 
          mean(exposure_101$pdays[exposure_101$es == 3]), "days")
  
  # Return results
  list(group_101 = exposure_101, group_102 = exposure_102)
}

# Run the test
# test_tidyverse_implementation()
