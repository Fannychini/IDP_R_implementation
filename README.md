## R implementation of Individualised Dispensing Patterns (IDP) methodology

This repository contains the R implementation of the Individualised Dispensing Patterns (IDP) methodology for estimating medication exposure in pharmaceutical claims data. 

The finalised implementation is located in `IDP_optimised.R` !!

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
exposure_by_drug(drug_code, macro_d, EndDate, pop_estimates, 
                 new_episode_threshold = 365, recent_exposure_window = 7, 
                 keep_tmp_variable = FALSE, percentile_col = "P_80",
                 DeathDate = NULL, DateSupplied_index = NULL, output_name = NULL)
```

- `drug_code` medication code/identifier
- `macro_d` dataframe with dispensing data
- `EndDate` study end date
- `pop_estimates` population estimates from `e_pop_estimate()`
- `new_episode_threshold` days threshold for new episodes (default: 365)
- `recent_exposure_window` days for recent exposure window (default: 7)

Optional:

- `keep_tmp_variable` whether to keep intermediate variables (default: FALSE)
- `percentile_col` other population estimate to use from pop_estimates (default: P_80)
- `DeathDate` death date, created if not supplied
- `DateSupplied_index` index date, uses earliest supply date if not supplied 
- `output_name` name for the output dataframe in global environment 


--\> **Returns** a dataframe with exposure periods for each person


### Use of the R implementation

```{r}
# Load libraries
library(tidyverse)
library(data.table)
library(lubridate)

# Source IDP functions
source("IDP_optimised.R")

# Load data
data <- read_csv("data.csv")

# Align with col name requirements
# data <- data %>%
#   rename(PPN = x1,               
#          Group = x2,             
#          DateSupplied = x3, 
#          Quantity = x4) 

# Calculate population estimates
pop_r <- e_pop_estimate("code", data) 

# Calculate exposure periods
exposure_r <- exposure_by_drug(drug_code = "code",
                               macro_d = data,
                               EndDate = max(data$DateSupplied[df$Group == 'code']), # or as.Date("yyyy-mm-dd"),
                               pop_estimates = pop_r, # output from e_pop_estimate
                               new_episode_threshold = 365,  
                               recent_exposure_window = 7) %>%
  mutate(PPN = factor(PPN)) %>%
  mutate(es = factor(es, levels=1:3, labels=c("Current", "Recent", "Former")))
```

### Output details

The output contains exposure periods with these key variables:

- `PPN` person identifier
- `DateSupplied` original dispensing date (kept for all rows)
- `date_of_supply` dispensing date (NA for recent/former periods)
- `Quantity` original quantity dispensed (NA for recent/former periods)
- `Group` medication/drug code
- `es` exposure status (1=current, 2=recent, 3=former)
- `start_date` start date of the exposure interval
- `end_date` end date of the exposure interval
- `pdays` person-days in interval 
- `ep_num{drug_code}` episode number for analysed drug (i.e. ep_num123)
- `first_{drug_code}_date` first dispensing date for analysed drug (i.e. first_123_date)
- `ep_st_date` episode start date
- `episode_dispensing` dispensing count in episode (NA for recent/former)
- `e_n` estimated exposure duration in days (NA for recent/former)
- `unique_ep_id` unique episode identifier
- `SEE1` next dispensing date (set as 9999-12-31 for last dispensing)
- `last` flag for last dispensing (1=yes)
- `rec_num` record number (increments with exposure state changes)

### Validation

*This R implementation has not been validated yet*.

Proposed flow for validation and comparison between SAS and R implementation:

```{r}
drug_code <- "item_code"
compare_columns <- c("PPN", "Group", paste0("ep_num", drug_code), 
                     paste0("first_", drug_code, "_date"), "ep_st_date", 
                     "start_date", "end_date", "episode_dispensing", "e_n", 
                     "unique_ep_id", "SEE1", "last", "es", "pdays", "rec_num")
                     
## Run R functions
# e_pop_estimate
pop_r <- e_pop_estimate("item_code", validation_data) 

# exposure_by_drug
exposure_r <- exposure_by_drug(drug_code = "item_code",
                               macro_d = validation_data,
                               EndDate = max(validation_data$DateSupplied[validation_data$Group == 'item_code']), 
                               pop_estimates = pop_r,
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

If you use this implementation in your research, please cite both the original paper and this R implementation:

   ```
   Bharat, C., Degenhardt, L., Pearson, S-A., et al. A data-informed approach using individualised dispensing patterns to estimate medicine exposure periods and dose from pharmaceutical claims data.
   Pharmacoepidemiol Drug Saf. 2023; 32(3): 352-365. doi:10.1002/pds.5567
   ```

   ```
   Franchini F. (2025). R implementation of the individualised dispensing patterns methodology.
   GitHub repository, https://github.com/fannychini/IDP_R_implementation
   ```

### License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
This package is released under the GNU General Public License v3.0 (GPL-3.0).

This means: 

- You are free to use, modify, and distribute this code 
- If you distribute modified versions, you must also distribute them under GPL-3.0 
- Any software that incorporates this code must also be released under GPL-3.0 
- Full license details can be found in the LICENSE file or at <https://www.gnu.org/licenses/gpl-3.0.html>
