# ---------------------------------------------- #
# Testing and validation of IDP-R implementation #
# ---------------------------------------------- #

## Libraries
library(tidyverse)
library(lubridate)
library(data.table)
library(patchwork)
library(beepr)

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
  # will use min date here
  group_by(PPN) %>%
  mutate(Date_of_Supply_index = min(Date_of_Supply)) %>%
  ungroup() %>%
  # also need status and will set all NA 
  mutate(DeathDate = as.Date(NA))

# Knowing that I wanted to set the medications as:
# "Levodrax" = regular ~30 day intervals
# "Cyclobine" = ~21 day intervals (14 days on, 7 off)
# "Biforalin" = ~14 day intervals
# "Quartifuse" = ~90 day intervals
# "Flexitol" = irregular intervals


## Single drug analysis ----
# Here we look at Flexitol

### e_pop_estimate ----
pop_output <- e_pop_estimate("Flexitol", df) %>%
  mutate(group = as.character(item_code))

### exposure_by_drug ----
exposure_result <- exposure_by_drug(drug_number = "Flexitol",
                                    macro_d = df,
                                    EndDate = as.Date("2023-12-31"), 
                                    # or max(df$Date_of_Supply) 
                                    combined_item_code3 = pop_output,
                                    new_episode_threshold = 365,  
                                    recent_exposure_window = 7) %>%
  mutate(PPN = factor(PPN)) %>%
  mutate(es = factor(es, levels=1:3, labels=c("Current", "Recent", "Former")))

### checks ----
# summary by exposure status (1=current, 2=recent, 3=former)
exposure_summary <- exposure_result %>%
  group_by(es) %>%
  summarise(patients = n_distinct(PPN),
            total_days = sum(pdays),
            avg_duration = mean(pdays))
#   es      patients total_days avg_duration
# 1 Current      200      65424        28.2 
# 2 Recent       200       3381         5.50
# 3 Former       175       3213         9.51

# visuals for 5 PPNs
exposure_5ppn <- exposure_result %>%
  distinct(PPN) %>%
  # sample 5 PPNs
  slice_sample(n = 5) %>%
  inner_join(exposure_result, by = "PPN")
  
ggplot(exposure_5ppn, aes(y=PPN, color=es)) +
  geom_segment(aes(x=start_date+1, xend=end_date, yend=PPN)) +
  geom_point(aes(x=Date_of_Supply)) +
  theme_minimal()


## SAS vs R ----
compare_columns <- c("PPN", "group", "EP1", "first_date", "ep_st_date", 
                     "start_date", "end_date", "episode_dispensing", "e_n", 
                     "unique_ep_id", "SEE1", "last", "es", "pdays", "rec_num")

### Run e_pop_estimate and exposure_by_drug ----
#### e_pop_estimate ----
pop_r <- e_pop_estimate("item_code", validation_data) %>%
  mutate(group = as.character(item_code))

#### exposure_by_drug ----
exposure_r <- exposure_by_drug(drug_number = "item_code",
                                    macro_d = validation_data,
                                    EndDate = as.Date("2023-12-31"), 
                                    # or max(validation_data$Date_of_Supply) 
                                    combined_item_code3 = pop_output,
                                    new_episode_threshold = 365,  
                                    recent_exposure_window = 7) %>%
  mutate(PPN = factor(PPN)) %>%
  mutate(es = factor(es, levels=1:3, labels=c("Current", "Recent", "Former")))

r_results <- exposure_r %>% arrange(PPN, start_date)

### SAS results ----
# load data
# arrange for comparison
sas_results <- exposure_sas %>% arrange(PPN, start_date)

### Check differences in both implementations ----
differences <- anti_join(r_results[, compare_columns], 
                         sas_results[, compare_columns])



## All codes at once analysis ----

### Get population estimates ----
unique_drugs <- unique(df$group)
combined_item_code <- data.frame()

for (drug in unique_drugs) {
  drug_percentiles <- e_pop_estimate(drug, df)
  combined_item_code <- rbind(combined_item_code, drug_percentiles)
}

combined_item_code3 <- combined_item_code %>%
  mutate(group = as.character(item_code))
#          N       P_20 P_50      P_60       P_70       P_80      P_85      P_90  item_code      group
# 20%  14800  0.9333333  1.0  1.033333  1.0333333  1.0666667  1.066667  1.100000   Levodrax   Levodrax
# 20%1 14800  0.4761905  0.5  0.500000  0.5238095  0.5238095  0.547619  0.547619  Cyclobine  Cyclobine
# 20%2 14800  0.9285714  1.0  1.000000  1.0714286  1.0714286  1.071429  1.071429  Biforalin  Biforalin
# 20%3 14800 84.0000000 90.0 92.000000 94.0000000 96.0000000 97.000000 99.000000 Quartifuse Quartifuse
# 20%4 14800  0.8333333  1.0  1.000000  1.0000000  1.5000000  1.500000  1.666667   Flexitol   Flexitol


### Calculate exposures ----
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
                                      recent_exposure_window = 7) %>%
    mutate(PPN = factor(PPN)) %>%
    mutate(es = factor(es, levels=1:3, labels=c("Current", "Recent", "Former")))
}

# this should generate several files in global environment, as per all drugs listed
# in this case, exposure_Biforalin, exposure_Cyclobine, etc.

# may need to fix warning at some point
# Warning message:
  # In `[.data.table`(macro_episodes, , `:=`(ep_st_date, first(Date_of_Supply)),  :
  # Invalid .internal.selfref detected and fixed by taking a (shallow) copy of the data.table 
  # so that := can add this new column by reference. At an earlier point, this data.table 
  # has been copied by R (or was created manually using structure() or similar).

### Analyse exposure for one drug ----
# or access through exposure_Levodrax
drug_example <- unique_drugs[1] # Levodrax
exposure_data <- get(paste0("exposure_", drug_example)) %>%
  mutate(PPN = factor(PPN)) %>%
  mutate(es = factor(es, levels=1:3, labels=c("Current", "Recent", "Former")))

# summary by exposure status (1=current, 2=recent, 3=former)
exposure_summary <- exposure_data %>%
  group_by(es) %>%
  summarise(patients = n_distinct(PPN),
            total_days = sum(pdays),
            avg_duration = mean(pdays))

print(exposure_summary)
#      es patients total_days avg_duration
# 1     1      200      71200        28.6 
# 2     2      193        782         1.84

check_143 <- exposure_data %>% filter(PPN %in% c('143')) %>% arrange(Date_of_Supply)

# check first PPN
pat_example <- exposure_data$PPN[1]
pat_timeline <- exposure_data %>%
  filter(PPN == pat_example) %>%
  arrange(start_date) %>% print()
# seems ok

ggplot(pat_timeline, aes(y=PPN, color=es)) +
  geom_segment(aes(x=start_date+1, xend=end_date, yend=PPN)) +
  geom_point(aes(x=Date_of_Supply)) +
  theme_minimal()

### Exposure distribution visuals ----
exposure_data %>%
  ggplot(aes(x = pdays, fill = es)) +
  geom_histogram(bins = 30) +
  facet_wrap(~es, scales = "free_y") +
  labs(title = paste("Distribution of exposure days for", drug_example),
       x = "person-days",
       y = "frequency",
       fill = "exposure type") +
  theme_minimal()







