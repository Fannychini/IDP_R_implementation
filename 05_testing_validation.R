# ---------------------------------------------- #
# Testing and validation of IDP-R implementation #
# ---------------------------------------------- #

## Libraries
library(tidyverse)
library(lubridate)
library(data.table)

## Set seed for reproducibility
set.seed(239)

## Source implementation and generated data
source("03_IDP_optimised.R")
#source("04_IDP_tidy_friendly.R") 
source("06_generate_data.R")

## Generate data ----
disp_data <- generate_diverse_data(n_patients_per_med = 200, n_dispenses_per_patient = 75)
head(disp_data)

## Prep data ----
# rename columns to match function requirements
df <- disp_data %>%
  rename(PPN = patient_id,               
         group = medication,             
         Date_of_Supply = dispense_date, 
         q_D = quantity_dispensed) %>%
  # also need Date_of_Supply_index (starting date for each patient)
  # will use min Date here
  group_by(PPN) %>%
  mutate(Date_of_Supply_index = min(Date_of_Supply)) %>%
  ungroup() %>%
  # also need status and will set all NA 
  mutate(DeathDate = as.Date(NA))

## Get population estimates ----
unique_drugs <- unique(df$group)
combined_item_code <- data.frame()

for (drug in unique_drugs) {
  drug_percentiles <- e_pop_estimate(drug, df)
  combined_item_code <- rbind(combined_item_code, drug_percentiles)
}

combined_item_code3 <- combined_item_code %>%
  mutate(group = as.character(item_code))
#          N       P_20 P_50       P_80      P_90  item_code      group
# 20%  14800  0.9333333  1.0  1.0666667  1.100000   Levodrax   Levodrax
# 20%1 14800  0.4761905  0.5  0.5238095  0.547619  Cyclobine  Cyclobine
# 20%2 14800  0.9285714  1.0  1.0714286  1.071429  Biforalin  Biforalin
# 20%3 14800 84.0000000 90.0 96.0000000 99.000000 Quartifuse Quartifuse
# 20%4 14800  0.8333333  1.0  1.5000000  1.533333   Flexitol   Flexitol

# Knowing that I wanted to set the medications as:
# "Levodrax" = regular ~30 day intervals
# "Cyclobine" = ~21 day intervals (14 days on, 7 off)
# "Biforalin" = ~14 day intervals
# "Quartifuse" = ~90 day intervals
# "Flexitol" = irregular intervals

## Calculate exposures ----
# study end date
end_date <- as.Date("2023-12-31")  

# process each drug
for (drug in unique_drugs) {
  # define output variable name
  output_name <- paste0("exposure_", drug)
  
  # calculate exposures using exposure_by_drug function
  exposure_result <- exposure_by_drug(drug_number = drug,
                                      macro_d = df,
                                      EndDate = end_date,
                                      combined_item_code3 = combined_item_code3,
                                      output_name = output_name,
                                      new_episode_threshold = 365,  
                                      recent_exposure_window = 7)
}

# may need to fix warning at some point
# Warning message:
  # In `[.data.table`(macro_episodes, , `:=`(ep_st_date, first(Date_of_Supply)),  :
  # Invalid .internal.selfref detected and fixed by taking a (shallow) copy of the data.table 
  # so that := can add this new column by reference. At an earlier point, this data.table 
  # has been copied by R (or was created manually using structure() or similar).

## Analyse exposure for one drug ----
drug_example <- unique_drugs[1] # Levodrax
exposure_data <- get(paste0("exposure_", drug_example))

# summary by exposure status (1=current, 2=recent, 3=former)
exposure_summary <- exposure_data %>%
  group_by(es) %>%
  summarise(patients = n_distinct(PPN),
            total_days = sum(pdays),
            avg_duration = mean(pdays))

print(exposure_summary)
#      es patients total_days avg_duration
# 1     1      200    70932.        28.5  
# 2     2      200      393.         0.599
# 3     3        1       -0.5       -0.5  

# problem with calculations for es=3
problematic_record <- exposure_data %>%
  filter(es == 3 & pdays < 0)
# PPN 143, start_date and end_date are the same 2023-06-01
# calculated pdays is -0.5 so there is an issue in the function
# this is because of date handling with as.numeric, need to change to as.integer to match SAS

# Rerun with corrected as.integer in function:
#      es patients total_days avg_duration
# 1     1      200      70590        28.4 
# 2     2      200        714         1.05
# 3     3        1          0         0   

check_143 <- exposure_data %>% filter(PPN %in% c('143')) %>% arrange(Date_of_Supply)

# PPN-level analysis 
pat_example <- exposure_data$PPN[1]
pat_timeline <- exposure_data %>%
  filter(PPN == pat_example) %>%
  select(PPN, start_date, end_date, es, pdays) %>%
  arrange(start_date)

print(pat_timeline) # seems ok

## Exposure distribution viuals ----
exposure_data %>%
  mutate(exposure_type = case_when(
    es == 1 ~ "current",
    es == 2 ~ "recent",
    es == 3 ~ "former")) %>%
  ggplot(aes(x = pdays, fill = exposure_type)) +
  geom_histogram(bins = 30) +
  facet_wrap(~exposure_type, scales = "free_y") +
  labs(title = paste("Distribution of exposure days for", drug_example),
       x = "person-days",
       y = "frequency",
       fill = "exposure type") +
  theme_minimal()
