# ------------------------------------- #
# R script for basic IDP implementation #
# ------------------------------------- #

## Library
library(tidyverse)
library(lubridate)

## Reproducibility
set.seed(239) 
sessionInfo()

# This is the basic R implementation that aims to reproduce the SAS 
# implementation as per https://github.com/c-bharat/IDP_exposure_method

## e_pop_estimate ----
e_pop_estimate <- function(item_code, macro_d) {
  
  # sort data by PPN and dispensing date
  medicine_data <- macro_d %>%
    arrange(PPN, Date_of_Supply) %>%
    # keep data for the specific medicine/class
    filter(group == item_code)
  
  # need empty df to collect continuous episodes and new episodes
  e_pope <- data.frame()
  new_episode <- data.frame()
  
  # process data row by row, sticking to SAS code 
  patients <- unique(medicine_data$PPN)
  for (patient in patients) {
    patient_data <- medicine_data %>% filter(PPN == patient)
    # initialise tracking variables
    last_time <- as.Date(NA)
    last_q <- NA
    t_nm1 <- as.Date(NA)
    q_nm1 <- NA
    new_episode_flag <- 0
    
    for (i in 1:nrow(patient_data)) {
      # update tracking variables
      t_nm1 <- last_time
      q_nm1 <- last_q
      last_q <- patient_data$q_D[i]
      last_time <- as.Date(patient_data$Date_of_Supply[i])
      
      # check if new episode (first dispensing or gap ≥ 365 days)
      if (is.na(t_nm1) || (patient_data$Date_of_Supply[i] - t_nm1) >= 365) {
        new_episode_flag <- 1
        # add row to new_episode dataset
        row <- patient_data[i, ]
        row$last_time <- last_time
        row$last_q <- last_q
        row$t_nm1 <- t_nm1
        row$q_nm1 <- q_nm1
        row$new_episode <- new_episode_flag
        
        new_episode <- rbind(new_episode, row)
        # reset flag
        new_episode_flag <- 0} 
      
      else {
        # add row to e_pope dataset (continuous episode)
        row <- patient_data[i, ]
        row$last_time <- last_time
        row$last_q <- last_q
        row$t_nm1 <- t_nm1
        row$q_nm1 <- q_nm1
        row$new_episode <- new_episode_flag
        
        e_pope <- rbind(e_pope, row)
      }
    }
  }
  
  # checks with first 50 rows of each dataset
  if (nrow(e_pope) > 0) {
    print("First 50 rows of e_pope:")
    print(head(e_pope, 50))}
  
  if (nrow(new_episode) > 0) {
    print("First 50 rows of new_episode:")
    print(head(new_episode, 50))}
  
  # calculate days per unit following a dispensing day
  if (nrow(e_pope) > 0) {
    e_pop_est <- e_pope %>%
      mutate(dif = as.numeric(Date_of_Supply - t_nm1),
             e_pop_est = dif / q_nm1) %>%
      group_by(PPN) %>%
      mutate(first_interval = row_number() == 1) %>%
      ungroup()
    
    # checks
    print("first 50 rows of e_pop_est:")
    print(head(e_pop_est, 50))
    
    # calculate percentiles (20th 50th 80th 90th)
    percentiles <- e_pop_est %>%
      summarise(N = n(),
                P_20 = quantile(e_pop_est, 0.20, na.rm = TRUE),
                P_50 = quantile(e_pop_est, 0.50, na.rm = TRUE),
                P_80 = quantile(e_pop_est, 0.80, na.rm = TRUE),
                P_90 = quantile(e_pop_est, 0.90, na.rm = TRUE)) 
    percentiles$item_code <- item_code # add parameter value
    
    # checks
    print("percentile:")
    print(percentiles)
    
    return(percentiles)
  } 
  
  else {
    warning("no continuous episodes found for this item_code")
    
    return(data.frame(N = 0, P_20 = NA, P_50 = NA,
                      P_80 = NA, P_90 = NA, item_code = item_code))
  }
}


## `exposure_by_drug` ----
# including necessary steps
exposure_by_drug <- function(drug_number, macro_d, EndDate, combined_item_code3, output_name) {
  # EndDate must be a date object
  EndDate <- as.Date(EndDate)
  
  # subset dispensing data for the specific drug
  macro_d_temp <- macro_d %>%
    filter(group == drug_number & Date_of_Supply >= Date_of_Supply_index) %>%
    mutate(all = 1) # using this for merging later
  
  # checks
  print("first 5 rows of filtered data:")
  print(head(macro_d_temp, 5))
  
  # generate population estimate = 80th percentile
  e_est_macro <- combined_item_code3 %>%
    filter(group == drug_number) %>%
    select(group, P_80) %>%
    mutate(all = 1) # using this for merging later
  
  # checks
  print("population estimate:")
  print(e_est_macro)
  
  # merge population estimate with dispensing data
  macro_d_temp <- merge(macro_d_temp, e_est_macro, by = "all")
  
  # initialise empty df for processed data
  macro_d1 <- data.frame()
  
  # will process data row by row by PPN for evaluation of exposure
  patients <- unique(macro_d_temp$PPN)
  
  # set counter for unique episode IDs
  unique_ep_id_counter <- 0  
  
  for (patient in patients) {
    patient_data <- macro_d_temp %>% 
      filter(PPN == patient) %>%
      arrange(Date_of_Supply) 
    
    # initialise patient-level tracking variables
    ep_num <- 0
    episode_dispensing <- 0
    recent_exp <- 7 # this is the recent exposure period length = 7d
    t_nm1 <- as.Date(NA)
    t_nm2 <- as.Date(NA)
    t_nm3 <- as.Date(NA)
    q_nm1 <- NA
    q_nm2 <- NA
    q_nm3 <- NA
    e_n <- 99999999
    last_time <- as.Date(NA)
    last_q <- NA
    first_date <- as.Date(NA)
    no_formerly_exposed <- 0
    
    patient_processed <- data.frame()
    
    for (i in 1:nrow(patient_data)) {
      row <- patient_data[i, ]
      
      # set 1st date if this is 1st dispensing
      if (i == 1) {
        first_date <- as.Date(row$Date_of_Supply)
      }
      
      # Save current values for later updating
      current_time <- as.Date(row$Date_of_Supply)
      current_q <- row$q_D
      
      # update variables from prior dispensings
      if (episode_dispensing >= 3) {
        t_nm3 <- as.Date(t_nm2)
        q_nm3 <- q_nm2
      }
      if (episode_dispensing >= 2) {
        t_nm2 <- as.Date(t_nm1)
        q_nm2 <- q_nm1
      }
      if (episode_dispensing >= 1) {
        t_nm1 <- as.Date(last_time)
        q_nm1 <- last_q
      }
      
      # check if this is a new episode or continuing
      if (ep_num == 0 || (current_time > (t_nm1 + e_n + recent_exp))) {
        # new this is new episode
        ep_num <- ep_num + 1
        unique_ep_id_counter <- unique_ep_id_counter + 1
        episode_dispensing <- 1
        t_nm1 <- as.Date(NA)
        t_nm2 <- as.Date(NA) 
        t_nm3 <- as.Date(NA)
        q_nm1 <- NA
        q_nm2 <- NA
        q_nm3 <- NA
        e_n <- current_q * row$P_80 # = initial estimate based on population
      } else {
        # else if this is continuing episode
        # check if formerly exposed
        if (current_time > (t_nm1 + e_n + recent_exp)) {
          ep_num <- ep_num + 1
          no_formerly_exposed <- 1
        }
        episode_dispensing <- episode_dispensing + 1
        
        # !! WEIGHTED FORMULA !!
        # check if this works now
        if (is.na(t_nm1)) {
          e_n <- current_q * row$P_80
        } else if (is.na(t_nm2)) {
          e_n <- current_q * ((3/6) * (as.numeric(current_time - t_nm1)) / q_nm1 + 
                                (2/6) * row$P_80 + 
                                (1/6) * row$P_80)
        } else if (is.na(t_nm3)) {
          e_n <- current_q * ((3/6) * (as.numeric(current_time - t_nm1)) / q_nm1 + 
                                (2/6) * (as.numeric(t_nm1 - t_nm2)) / q_nm2 + 
                                (1/6) * row$P_80)
        } else {
          e_n <- current_q * ((3/6) * (as.numeric(current_time - t_nm1)) / q_nm1 + 
                                (2/6) * (as.numeric(t_nm1 - t_nm2)) / q_nm2 + 
                                (1/6) * (as.numeric(t_nm2 - t_nm3)) / q_nm3)
        }
      }
      
      # update variables for next iteration
      last_time <- as.Date(current_time)
      last_q <- current_q
      
      # add processed row to df
      row$ep_num <- ep_num
      row$episode_dispensing <- episode_dispensing
      row$t_nm1 <- t_nm1
      row$t_nm2 <- t_nm2
      row$t_nm3 <- t_nm3
      row$q_nm1 <- q_nm1
      row$q_nm2 <- q_nm2
      row$q_nm3 <- q_nm3
      row$e_n <- e_n
      row$last_time <- last_time
      row$last_q <- last_q
      row$first_date <- first_date
      row$unique_ep_id <- unique_ep_id_counter
      row$recent_exp <- recent_exp
      row$no_formerly_exposed <- no_formerly_exposed
      
      patient_processed <- rbind(patient_processed, row)
    }
    
    macro_d1 <- rbind(macro_d1, patient_processed)
  }
  
  # need death and censoring information
  macro_d1 <- macro_d1 %>%
    mutate(death = ifelse(!is.na(DeathDate), 1, 0),
           cens_date = pmin(DeathDate, EndDate, na.rm = TRUE))
  
  # get date of next dispensing
  macro_d1 <- macro_d1 %>%
    arrange(PPN, desc(Date_of_Supply)) %>%
    group_by(PPN) %>%
    mutate(PPN1 = lag(PPN),
           SEE1 = lag(Date_of_Supply),
           EP1 = lag(ep_num)) %>%
    # here this is each PPN first row = last dispensing chronologically
    mutate(SEE1 = if_else(row_number() == 1, as.Date("9999-12-31"), as.Date(SEE1)),
           last = if_else(row_number() == 1, 1L, 0L),
           EP1 = if_else(row_number() == 1, NA_integer_, EP1)) %>%
    ungroup() %>%
    # go back to chronological order
    arrange(PPN, Date_of_Supply)
  
  # checs again
  print("Sample of macro_d1 with next dispensing dates:")
  print(head(macro_d1, 5))
  
  # create the exposure based on current / recent / former
  macro_episodes <- data.frame()
  post_death <- data.frame()
  
  for (i in 1:nrow(macro_d1)) {
    row <- macro_d1[i, ]
    
    # endpoints accounting for next dispensing, death, and study end
    current_end <- min(as.Date(row$Date_of_Supply) + row$e_n, 
                       as.Date(row$SEE1) - 1,
                       as.Date(row$DeathDate),
                       as.Date(EndDate),
                       na.rm = TRUE)
    
    recent_end <- min(as.Date(row$Date_of_Supply) + row$e_n + row$recent_exp,
                      as.Date(row$SEE1) - 1,
                      as.Date(row$DeathDate),
                      as.Date(EndDate),
                      na.rm = TRUE)
    
    former_end <- min(as.Date(row$SEE1) - 1,
                      as.Date(row$DeathDate),
                      as.Date(EndDate), 
                      na.rm = TRUE)
    
    # error with dates, check if dispensing is after death/end date
    if (row$Date_of_Supply > min(row$DeathDate, EndDate, na.rm = TRUE)) {
      # add to post_death dataset
      post_death <- rbind(post_death, row)
    } else {
      # es=1 -> current exposure
      interval <- row
      interval$es <- 1L
      interval$start_date <- as.Date(row$Date_of_Supply) - 1
      interval$end_date <- as.Date(current_end)
      macro_episodes <- rbind(macro_episodes, interval)
      
      # es=2 -> recent exposure + if time is remaining
      if (current_end < min(row$SEE1 - 1, row$DeathDate, EndDate, na.rm = TRUE)) {
        interval <- row
        interval$es <- 2L
        interval$start_date <- as.Date(current_end) + 1
        interval$end_date <- as.Date(recent_end)
        macro_episodes <- rbind(macro_episodes, interval)
        
        # es=3 -> former exposure 
        if (recent_end < min(row$SEE1 - 1, row$DeathDate, EndDate, na.rm = TRUE)) {
          interval <- row
          interval$es <- 3L
          interval$start_date <- as.Date(recent_end) + 1 
          interval$end_date <- as.Date(former_end)
          macro_episodes <- rbind(macro_episodes, interval)
        }
      }
    }
  }
  
  # sort episodes by PPN and start date
  macro_episodes <- macro_episodes %>%
    arrange(PPN, start_date)
  
  # cleaning and fixing issues:
  # need to clean variables for non active periods
  # fix needed as as.Date not preserved 
  macro_episodes2 <- macro_episodes %>%
    # make sure dates are dates
    mutate(date_of_supply = if_else(es == 2L | es == 3L, as.Date(NA), as.Date(Date_of_Supply)),
           q_D = if_else(es == 2L | es == 3L, NA_real_, q_D),
           e_n = if_else(es == 2L | es == 3L, NA_real_, e_n),
           episode_dispensing = if_else(es == 2L | es == 3L, NA_integer_, episode_dispensing),
           start_t = as.numeric(start_date - (Date_of_Supply_index - 1)),
           end_t = as.numeric(end_date - (Date_of_Supply_index - 1)),
           # extra check all date columns need to be as.Date
           t_nm1 = as.Date(t_nm1),
           t_nm2 = as.Date(t_nm2),
           t_nm3 = as.Date(t_nm3),
           last_time = as.Date(last_time),
           first_date = as.Date(first_date)) %>%
    group_by(PPN, ep_num) %>%
    mutate(ep_st_date = as.Date(first(Date_of_Supply))) %>%
    ungroup()
  
  # create final output df
  ep_num_col <- paste0("ep_num", drug_number)
  first_date_col <- paste0("first_", drug_number, "_date")
  
  # get initial output with PPN groups
  output_data <- macro_episodes2 %>%
    mutate(pdays = as.numeric(end_date - start_date)) %>%
    group_by(PPN) %>%
    mutate(lag_es = lag(es)) %>%
    ungroup()
  
  # adding this step as this caused issues with data types later on
  output_data[[ep_num_col]] <- NA_integer_
  output_data[[first_date_col]] <- as.Date(NA)
  
  # fill in episode numbers and dates correctly
  for (p in unique(output_data$PPN)) {
    patient_rows <- which(output_data$PPN == p)
    
    # set initial values for first row
    output_data[[ep_num_col]][patient_rows[1]] <- 1L
    output_data$rec_num[patient_rows[1]] <- 1L
    output_data[[first_date_col]][patient_rows[1]] <- as.Date(output_data$Date_of_Supply[patient_rows[1]])
    
    # other rows
    if (length(patient_rows) > 1) {
      for (i in 2:length(patient_rows)) {
        idx <- patient_rows[i]
        prev_idx <- patient_rows[i-1]
        
        # first date value and also must be a date
        output_data[[first_date_col]][idx] <- as.Date(output_data[[first_date_col]][prev_idx])
        
        # increase record number if exposure state changes
        if (output_data$es[idx] != output_data$es[prev_idx]) {
          output_data$rec_num[idx] <- output_data$rec_num[prev_idx] + 1L
        } else {
          output_data$rec_num[idx] <- output_data$rec_num[prev_idx]
        }
        
        # icrease episode number if returning to current exposure after former exposure
        if (output_data$es[idx] == 1 && output_data$lag_es[idx] == 3) {
          output_data[[ep_num_col]][idx] <- output_data[[ep_num_col]][prev_idx] + 1L
        } else {
          output_data[[ep_num_col]][idx] <- output_data[[ep_num_col]][prev_idx]
        }
      }
    }
  }
  
  # don't need original ep_num -- avoiding confusion
  output_data <- output_data %>%
    select(-ep_num)
  
  # checks
  print("sample of final output:")
  print(head(output_data, 10))
  
  # assign to global environment with provided output name
  assign(output_name, output_data, envir = .GlobalEnv)
  
  return(output_data)
}



## Testing ----

### Generate simple data ----
create_test_data <- function() {
  test_data <- data.frame(
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
  
  return(test_data)
}


### Implementation testing ----

test_e_pop_estimate <- function() {
  cat("\ne_pop_estimate testing\n")
  # generate test data
  test_data <- create_test_data()
  
  # run the function for group 101 and 102
  percentiles_101 <- e_pop_estimate(101, test_data)
  percentiles_102 <- e_pop_estimate(102, test_data)
  
  # show results
  cat("\npercentiles for group 101:\n")
  print(percentiles_101)
  cat("\npercentiles for group 102:\n")
  print(percentiles_102)
  
  # combine data
  combined_item_code <- rbind(percentiles_101, percentiles_102)
  combined_item_code3 <- combined_item_code %>%
    mutate(group = as.numeric(item_code))
  
  cat("\ncombined percentiles:\n")
  print(combined_item_code3)
  
  return(combined_item_code3)
}


test_exposure_by_drug <- function(combined_item_code3) {
  cat("\nexposure_by_drug testing\n")
  
  # generate test data
  test_data <- create_test_data()
  
  # define study end date
  EndDate <- as.Date("2021-12-31")
  
  # Run the function for group 101
  cat("\nexposure calculation for group 101\n")
  exposure_101 <- exposure_by_drug(101, test_data, EndDate, combined_item_code3, "output_101")
  
  # Run the function for group 102
  cat("\nexposure calculation for group 102\n")
  exposure_102 <- exposure_by_drug(102, test_data, EndDate, combined_item_code3, "output_102")
  
  # Return the results
  return(list(group_101 = exposure_101, group_102 = exposure_102))
}


verify_results <- function(exposure_results) {
  cat("\ncheck results\n")
  
  # the results are in exp_101 and exp_102
  exposure_101 <- exposure_results$group_101
  exposure_102 <- exposure_results$group_102
  
  # check exposure
  cat("\nExposure status for group 101\n")
  print(table(exposure_101$es))
  
  cat("\nExposure status for group 102\n")
  print(table(exposure_102$es))
  
  # PPM counts
  cat("\nPPN count for group 101", length(unique(exposure_101$PPN)), "\n")
  cat("\nPPN count for group 102", length(unique(exposure_102$PPN)), "\n")
  
  # check a specific PPN 
  cat("\nexposure periods for PPN 1001 (group 101)\n")
  print(exposure_101 %>% filter(PPN == 1001) %>% select(PPN, start_date, end_date, es))
  
  # check if PPN 1002 (deceased) has correct exposure truncation
  cat("\nexposure PPN 1002 with death set as 2021-06-30\n")
  print(exposure_101 %>% filter(PPN == 1002) %>% select(PPN, start_date, end_date, es, DeathDate))
  
  # check episode numbers are assigned correctly
  ep_num_col <- paste0("ep_num", 101)
  cat("\nepisode counts for group 101\n")
  print(table(exposure_101[[ep_num_col]]))
  
  # calculate average exposure duration
  cat("\naverage duration of current exposure (es=1) for group 101", 
      mean(exposure_101$pdays[exposure_101$es == 1]), "days\n")
  cat("\naverage duration of recent exposure (es=2) for group 101", 
      mean(exposure_101$pdays[exposure_101$es == 2]), "days\n")
  cat("\naverage duration of former exposure (es=3) for group 101", 
      mean(exposure_101$pdays[exposure_101$es == 3]), "days\n")
}



run_complete_test <- function() {
  # e_pop_estimate
  combined_item_code3 <- test_e_pop_estimate()
  
  # exposure_by_drug
  exposure_results <- test_exposure_by_drug(combined_item_code3)
  
  # check results
  verify_results(exposure_results)
}

run_complete_test()


## Results ----

### e_pop_estimate expected outputs ----

# group 101 percentiles
# n: 7
# P_20: should be ~ 0.9-1.1 
# P_50: should be ~ 1.0-1.2
# P_80: should be ~ 1.3-1.5
# P_90: should be ~ 1.4-1.7

# group 102 percentiles
# n: 3
# P_20: should be ~ 0.9-1.1
# P_50: should be ~ 1.0-1.2
# P_80: should be ~ 1.3-1.5
# P_90: should be ~ 1.4-1.7


### exposure_by_drug expected outputs ----

## expected summary for patient 1001 set for regular monthly dispensing:
# currently exposed periods should be 5(es=1), ~ 30 days each
# recently exposed periods (es=2) ~  7 days each
# formerly exposed periods (es=3) should be of varying lengths
# should be 1 episode (ep_num101 = 1) for all periods

## expected summary for patient 1002 (with death date):
# 2 normal dispensings, then a large gap -- need to generate visualisation for that
# the 3rd dispensing is expected to start a new episode (ep_num101 = 2)
# should not see any exposure periods after death date (2021-06-30) 

## expected summary for patient 1003 (set with varied quantities):
# period 1 expected to use pop P_80 estimate
# period 2 use weighted formula with 50% individual data
# period 3 use weighted formula with 83% individual data
# - Fourth+ periods should use 100% individual data from previous 3 dispensings
# - Varying dispensing quantities should result in different exposure duration

## expected summary for patient 1004 (set with different medication):
# shoudl fall within group 102, exposure periods should be longer for 60days


### overall expected outputs ----

# Total exposure periods:
# - group 101 should be ~ 30-35 periods
# - group 102 should be ~ 10-12 periods

# Exposure distribution:
# es1, es2 and es3 should have each approx 33% 

# Episode counts:
# - PPN 1001: 1 episode
# - PPN 1002: 2 episodes, becayse of 365+ day gap
# - PPN 1003: 2 episodes, because of gap in dispensing
# - PPN 1004: 1 episode
