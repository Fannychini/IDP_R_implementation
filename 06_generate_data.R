# ----------------------------------------- #
# Generate fake dispensing data for testing #
# ----------------------------------------- #

## Libraries
library(tidyverse)
library(lubridate)
library(data.table)

## Set seed for reproducibility
set.seed(239)

## Generate disp_data df for debugging and testing
# The idea is to try and generate 'realistic' dispensing intervals for different medication types
# Will generate fake medications based on some common dispensing interval patterns

generate_diverse_data <- function(n_patients_per_med, n_dispenses_per_patient) {
  # medication types
  medications <- c(
    "Levodrax",       # Regular ~30 day intervals
    "Cyclobine",      # ~21 day intervals (14 days on, 7 off)
    "Biforalin",      # ~14 day intervals
    "Quartifuse",     # ~90 day intervals
    "Flexitol")       # Irregular intervals
  
  result <- data.frame()
  
  patient_id <- 1
  for (med in medications) {
    for (i in 1:n_patients_per_med) {
      start_date <- as.Date("2023-01-01") + days(sample(0:10, 1))
      
      if (med == "Levodrax") { # = daily oral
        # regular monthly dispensing with small variations
        intervals <- c(NA, round(rnorm(n_dispenses_per_patient - 1, mean = 30, sd = 2)))
        quantity <- rep(30, n_dispenses_per_patient)
        strength_mg <- rep(50, n_dispenses_per_patient)
      } 
      else if (med == "Cyclobine") { # = cyclic oral
        # 21-day cycle medication
        intervals <- c(NA, round(rnorm(n_dispenses_per_patient - 1, mean = 21, sd = 1.5)))
        quantity <- rep(42, n_dispenses_per_patient)  # 2 pills per day for 21 days
        strength_mg <- rep(25, n_dispenses_per_patient)
      }
      else if (med == "Biforalin") { # = biweekly cyclic
        # biweekly dispensing
        intervals <- c(NA, round(rnorm(n_dispenses_per_patient - 1, mean = 14, sd = 1)))
        quantity <- rep(14, n_dispenses_per_patient)
        strength_mg <- rep(100, n_dispenses_per_patient)
      }
      else if (med == "Quartifuse") { # = quarterly IV
        # quarterly dispensing with more variation
        intervals <- c(NA, round(rnorm(n_dispenses_per_patient - 1, mean = 90, sd = 7)))
        quantity <- rep(1, n_dispenses_per_patient)  # Single infusion
        strength_mg <- rep(200, n_dispenses_per_patient)
      }
      else if (med == "Flexitol") { # = variable dosing
        # variable dispensing with some missing doses and early refills
        base_intervals <- c(NA, rep(30, n_dispenses_per_patient - 1))
        # include variability - some early refills, some late
        variation <- sample(c(-10, -5, 0, 5, 15, 25), n_dispenses_per_patient - 1, replace = TRUE, 
                            prob = c(0.1, 0.2, 0.4, 0.1, 0.1, 0.1))
        intervals <- c(NA, base_intervals[-1] + variation)
        quantity <- rep(30, n_dispenses_per_patient)
        strength_mg <- rep(75, n_dispenses_per_patient)
      }
      
      # calculate actual dates from intervals
      dates <- vector("list", length(intervals))
      dates[[1]] <- start_date
      
      for (j in 2:length(intervals)) {
        dates[[j]] <- dates[[j-1]] + days(intervals[j])
      }
      
      dates <- do.call(c, dates)
      
      # create patient data
      patient_data <- data.frame(
        patient_id = patient_id,
        medication = med,
        dispense_date = dates,
        quantity_dispensed = quantity,
        strength_mg = strength_mg,
        age = sample(30:75, 1))
      
      result <- rbind(result, patient_data)
      patient_id <- patient_id + 1
    }
  }
  
  # age groups
  result <- result %>%
    mutate(age_group = case_when(
      age < 45 ~ "Under 45",
      age >= 45 & age < 65 ~ "45-64",
      age >= 65 ~ "65 and above"))
  
  return(result)
}
