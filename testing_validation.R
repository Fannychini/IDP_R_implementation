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
source("IDP_optimised.R")
source("generate_data.R")

## Generate data ----
disp_data <- generate_diverse_data(n_patients_per_med = 50, n_dispenses_per_patient = 5)
head(disp_data)
# write.csv(disp_data, file = "data_april.csv")

## Prep data ----
# Data requirements:
# individual-level daily dispensing data, where dispensings on the same day for the same medicine(s) are combined/summed together (daily total dispensed)

# Col names if want to match function requirements:
# rename columns
df <- disp_data %>%
  rename(PPN = patient_id,
         Group = medication,
         DateSupplied = dispense_date,
         Quantity = quantity_dispensed)


## Single drug analysis ----
# Here we look at Levodrax

### e_pop_estimate ----
pop_output <- e_pop_estimate("Levodrax", df)

### exposure_by_drug ----
exposure_result <- exposure_by_drug(drug_code = "Levodrax",
                                    macro_d = df,
                                    EndDate = max(df$DateSupplied[df$Group == 'Levodrax']),
                                    pop_estimates = pop_output,
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
# 1 Current       50       6944        27.8
# 2 Recent        37        134         2.73
# 3 Former         3         18         6

# # visuals for 5 PPNs
# exposure_5ppn <- exposure_result %>%
#   distinct(PPN) %>%
#   # sample 5 PPNs
#   slice_sample(n = 5) %>%
#   inner_join(exposure_result, by = "PPN")
#
# ggplot(exposure_5ppn, aes(y=PPN, color=es)) +
#   geom_segment(aes(x=start_date+1, xend=end_date, yend=PPN)) +
#   geom_point(aes(x=DateSupplied)) +
#   theme_minimal()

ggplot(exposure_result, aes(y=PPN, color=es)) +
  geom_segment(aes(x=start_date+1, xend=end_date, yend=PPN)) +
  geom_point(aes(x=DateSupplied)) +
  theme_minimal()


## SAS vs R ----
drug_code <- "code"
compare_columns <- c("PPN", "Group", paste0("ep_num", drug_code),
                     paste0("first_", drug_code, "_date"), "ep_st_date",
                     "start_date", "end_date", "episode_dispensing", "e_n",
                     "unique_ep_id", "SEE1", "last", "es", "pdays", "rec_num")

### Run e_pop_estimate and exposure_by_drug ----
#### e_pop_estimate ----
pop_r <- e_pop_estimate("code", validation_data) %>%
  mutate(group = as.character(item_code))

#### exposure_by_drug ----
exposure_r <- exposure_by_drug(drug_code = "code",
                               macro_d = validation_data,
                               EndDate = max(validation_data$Date_of_Supply),
                               pop_estimates = pop_r,
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

combined_item_code3 <- combined_item_code
#          N       P_20 P_50      P_60       P_70       P_80      P_85      P_90      group
# 20%  14800  0.9333333  1.0  1.033333  1.0333333  1.0666667  1.066667  1.100000   Levodrax
# 20%1 14800  0.4761905  0.5  0.500000  0.5238095  0.5238095  0.547619  0.547619  Cyclobine
# 20%2 14800  0.9285714  1.0  1.000000  1.0714286  1.0714286  1.071429  1.071429  Biforalin
# 20%3 14800 84.0000000 90.0 92.000000 94.0000000 96.0000000 97.000000 99.000000 Quartifuse
# 20%4 14800  0.8333333  1.0  1.000000  1.0000000  1.5000000  1.500000  1.666667   Flexitol


### Calculate exposures ----
# study end date
end_date <- as.Date("2023-12-31")

# process each drug
for (drug in unique_drugs) {
  # define output variable name
  output_name <- paste0("exposure_", drug)

  # calculate exposures using exposure_by_drug function
  exposure_result <- exposure_by_drug(drug_code = drug,
                                      macro_d = df,
                                      EndDate = end_date,
                                      pop_estimates = pop_estimates,
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

check_143 <- exposure_data %>% filter(PPN %in% c('143')) %>% arrange(DateSupplied)

# check first PPN
pat_example <- exposure_data$PPN[1]
pat_timeline <- exposure_data %>%
  filter(PPN == pat_example) %>%
  arrange(start_date) %>% print()
# seems ok

ggplot(pat_timeline, aes(y=PPN, color=es)) +
  geom_segment(aes(x=start_date+1, xend=end_date, yend=PPN)) +
  geom_point(aes(x=DateSupplied)) +
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







