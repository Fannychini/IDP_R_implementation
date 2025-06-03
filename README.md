## R implementation of Individualised Dispensing Patterns (IDP) methodology

This repository contains the R implementation of the Individualised Dispensing Patterns (IDP) methodology for estimating medication exposure in pharmaceutical claims data. 

!! The finalised implementation is located in `03_IDP_optimised.R`. !!

### Overview

The IDP methodology provides a data-driven approach to estimate medication exposure periods based on individual persons' dispensing
patterns. 

This R implementation is based on the methodology described in: 

> Bharat, C., Degenhardt, L., Pearson, S-A., et al. A data-informed approach using individualised dispensing patterns to estimate medicine exposure periods and dose from pharmaceutical claims data. Pharmacoepidemiol Drug Saf. 2023; 32(3): 352-365. <doi:10.1002/pds.5567>

### Required formatting

The input data must have these columns:

- `PPN` person identifier
- `Group` dedication code/identifier
- `DateSupplied` dispensing date
- `Quantity` quantity dispensed
- `DateSupplied_index` starting date for analysis
    - ! if not supplied, calculated as per min `DateSupplied` for each PPN
- `DeathDate` date of death (can be NA)
    - ! if not supplied, generated with NA as date

### Key functions

#### `e_pop_estimate()`

Calculates population-level estimates of medication usage patterns.

```{r}
e_pop_estimate(drug_code, macro_d)
```

- `drug_code` medication code/identifier
- `macro_d` dataframe with dispensing data

--\> **Returns** a dataframe with various percentile values for days per unit

#### `exposure_by_drug()`

Calculates individual person exposure periods based on dispensing data.

```{r}
exposure_by_drug(drug_code, macro_d, EndDate, combined_item_code3, 
                 output_name = NULL, new_episode_threshold = 365, 
                 recent_exposure_window = 7, keep_tmp_variable = FALSE)
```

- `drug_code` medication code/identifier
- `macro_d` dataframe with dispensing data
- `EndDate` study end date
- `combined_item_code3` population estimates from `e_pop_estimate()`
- `output_name` name for the output dataframe in global environment (optional)
- `new_episode_threshold` days threshold for new episodes (default: 365)
- `recent_exposure_window` days for recent exposure window (default: 7)
- `keep_tmp_variable` whether to keep intermediate variables (default: FALSE)

--\> **Returns** a dataframe with exposure periods for each person


### Use of the R implementation

```{r}
# Load libraries
library(tidyverse)
library(data.table)
library(lubridate)

# Source R implementation
source("03_IDP_optimised.R")

# Load data
data <- read_csv("data.csv")

# Align with col name requirements
# data <- data %>%
#   rename(PPN = x1,               
#          Group = x2,             
#          DateSupplied = x3, 
#          Quantity = x4,
#          # optional
#          DeathDate = x5,
#          DateSupplied_index = x6)

# Calculate population estimates
pop_r <- e_pop_estimate("code", data) 

# Calculate exposure periods
exposure_r <- exposure_by_drug(drug_code = "code",
                               macro_d = data,
                               EndDate = as.Date("2023-12-31"), # or max(data$DateSupplied[df$Group == 'code'])
                               combined_item_code3 = pop_r, # output from e_pop_estimate
                               new_episode_threshold = 365,  
                               recent_exposure_window = 7) %>%
  mutate(PPN = factor(PPN)) %>%
  mutate(es = factor(es, levels=1:3, labels=c("Current", "Recent", "Former")))
```

### Output details

The output contains exposure periods with these key variables:

- `PPN` person identifier
- `es` exposure status (1=current, 2=recent, 3=former)
- `start_date` start date of the exposure interval
- `end_date` end date of the exposure interval
- `pdays` person-days of exposure
- `ep_num1` episode number for analysed drug (episodes are periods of continuous treatment)
- `first_1_date` first dispensing date for analysed drug
- `ep_st_date` episode start date
- `episode_dispensing` number of dispensings in this episode
- `e_n` calculated exposure duration in days
- `unique_ep_id` unique episode identifier
- `SEE1` next dispensing date
- `last` flag for last dispensing (1=yes)
- `rec_num` record number (increments with exposure state changes)

### Validation

Please note that *this R implementation has not been validated yet*.

Proposed flow for validation and comparison between SAS and R implementation:

```{r}
compare_columns <- c("PPN", "group", "ep_num1", "first_1_date", "ep_st_date", 
                     "start_date", "end_date", "episode_dispensing", "e_n", 
                     "unique_ep_id", "SEE1", "last", "es", "pdays", "rec_num")

## Run R functions
# e_pop_estimate
pop_r <- e_pop_estimate("item_code", validation_data) 

# exposure_by_drug
exposure_r <- exposure_by_drug(drug_code = "item_code",
                               macro_d = validation_data,
                               EndDate = as.Date("2023-12-31"), 
                               combined_item_code3 = pop_r,
                               new_episode_threshold = 365,  
                               recent_exposure_window = 7) %>%
  mutate(PPN = factor(PPN)) %>%
  mutate(es = factor(es, levels=1:3, labels=c("Current", "Recent", "Former")))

r_results <- exposure_r %>% arrange(PPN, start_date)

## SAS code
# load data

# arrange for comparison
sas_results <- exposure_sas %>% arrange(PPN, start_date)

## Check differences
differences <- anti_join(r_results[, compare_columns], 
                         sas_results[, compare_columns])
```

### Citation

If you use this implementation in your research, please cite both the original paper and this repository:

> Bharat, C., Degenhardt, L., Pearson, S-A., et al. A data-informed
> approach using individualised dispensing patterns to estimate medicine
> exposure periods and dose from pharmaceutical claims data.
> Pharmacoepidemiol Drug Saf. 2023; 32(3): 352-365.
> <doi:10.1002/pds.5567>

> Franchini, F. (2025). R implementation of the individualised
> dispensing patterns methodology. GitHub repository
> <https://github.com/fannychini/IDP_R_implementation>

### License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
This package is released under the GNU General Public License v3.0 (GPL-3.0).

This means: 

- You are free to use, modify, and distribute this code 
- If you distribute modified versions, you must also distribute them under GPL-3.0 
- Any software that incorporates this code must also be released under GPL-3.0 
- Full license details can be found in the LICENSE file or at <https://www.gnu.org/licenses/gpl-3.0.html>
