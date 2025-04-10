# ---------------------------
# Purpose of script: Check the output of gcplyr against rolling regression.
#
# What this script does:
# 1. Reads in a Run file
# 2. Calculates all the growth metrics from gcplyr
# 3. Calculates growth rate from rolling regression two ways
# 4. Checks how close they are
#
# Author: Dr. Daniel Padfield
#
# Date Created:  2025-03-13
#
# Copyright (c) Daniel Padfield, 2025
# This code is licensed under a modified MIT non-AI license. The code and any modifications made to it may not be used for the purpose of training or improving machine learning algorithms, including but not limited to artificial intelligence, natural language processing, or data mining. This condition applies to any derivatives, modifications, or updates based on the Software code. Any usage of the Software in an AI-training dataset is considered a breach of this License.
# The full license can be found here: https://github.com/padpadpadpad/non-ai-licenses/blob/main/NON-AI-MIT
#
# ---------------------------
#
# Notes:
#
# ---------------------------

# if librarian is not installed, install it
if (!requireNamespace("librarian", quietly = TRUE)){
  install.packages("librarian")
}

# if BiocManager is not installed, install it
if (!requireNamespace("BiocManager", quietly = TRUE)){
install.packages("BiocManager")
}

# if Biobase is not installed, install it from Bioconductor
if (!requireNamespace("Biobase", quietly = TRUE)){
BiocManager::install("Biobase")
}

# load packages
librarian::shelf(tidyverse, gcplyr, growthrates, GGally)

## ---------------------------

# load in a dataset

# name of files to read in
file_names <- list.files('data/processed', full.names = TRUE, pattern = 'od_data')

# read in data
d <- file_names[1] %>% map_df(read.csv)

# check the number of unique combinations there should be
select(d, run, well, temp, id, serial_no) %>%
  distinct() %>%
  nrow() 

# set up your grouping variables to get each individual well for each plate
# this bit is interactive
groupings <- c("run", "temp", "well", "serial_no", "id")

#--------------------#
# calculate blank ####
#--------------------#

# method 1: Calculate well specific blank from the first 3 times points (12 minutes)

# take the first three time points 
# 8 minutes
# average per well per temperature
d_blank <- d %>%
  group_by(across(all_of(groupings))) %>%
  filter(measurement_time_hr %in% sort(measurement_time_hr)[1:3]) %>%
  summarise(ave_blank = mean(raw_absorbance),
            sd_blank = sd(raw_absorbance),
            .groups = 'drop')

hist(d_blank$ave_blank)

# method 2: Calculate blank from average of all blank wells in each well

outside_wells <- c('A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12', 'B1', 'B12', 'C1', 'C12', 'D1', 'D12', 'E1', 'E12', 'F1', 'F12', 'G1', 'G12', 'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12')

# use the median
d_blank <- d %>%
  filter(! well %in% outside_wells) %>%
  filter(substr(id, 1, 1) == "X") %>%
  group_by(across(all_of(groupings[!groupings %in% c('well', 'id')]))) %>%
  summarise(ave_blank = median(raw_absorbance),
            sd_blank = sd(raw_absorbance),
            .groups = 'drop')

hist(d_blank$ave_blank)
# we will use this method as it is more commonly used, especially when the well is not completely empty of bacteria

#-----------------#
# wrangle data ####
#-----------------#

# add the blank to the data
d <- left_join(d, d_blank)

# make corrected absorbance
d <- mutate(d, od_cor = raw_absorbance - ave_blank)

# set negative values to 0
d <- mutate(d, od_cor = ifelse(od_cor < 0, 0, od_cor))

# remove all wells that are blanks - need to edit this code maybe
d_filt <- filter(d, substr(id, 1, 1) != "X")

# do a quick plot
d_filt %>%
  filter(well == sample(well, 1))  %>%
  ggplot(aes(x = measurement_time_hr, y = log(od_cor), group = interaction(well, run))) +
  geom_point(alpha = 0.3) +
  facet_wrap(~temp) +
  theme_bw()

# calculate the expected number of rows in the output - 1 per group
n_expected <- select(d_filt, all_of(groupings)) %>%
  distinct() %>%
  nrow()
n_expected

#----------------------------------------#
# calculate model-free growth metrics ####
#----------------------------------------#

# this uses gcplyr: https://mikeblazanin.github.io/gcplyr/index.html

# calculate per-capita derivative - this is the slope across a window of time
# a key parameter here is window_width_n, which is the number of time points to use in the calculation

# this group by is the important step here, needs to contain all the variables to get a single well
d_filt <- d_filt %>%
  filter(od_cor > 0.001) %>%
  group_by(across(all_of(groupings))) %>%
  mutate(deriv_percap = calc_deriv(x = measurement_time_hr, 
                                   y = od_cor,
                                   percapita = TRUE, 
                                   blank=0,
                                   window_width_n = 7,
                                   trans_y = "log")) %>%
  ungroup()

# calculate growth metrics
d_sum <- d_filt %>%
  group_by(across(all_of(groupings))) %>%
  summarise(lag_time = lag_time(y = od_cor, 
                                x = measurement_time_hr,
                                deriv = deriv_percap,
                                blank=0),
            max_percap = max_gc(deriv_percap, na.rm = TRUE),
            max_percap_time = extr_val(measurement_time_hr, which_max_gc(deriv_percap)),
            max_percap_dens = extr_val(od_cor, which_max_gc(deriv_percap)),
            min_dens = min_gc(od_cor),
            max_dens = max_gc(od_cor, na.rm = TRUE),
            max_time = extr_val(measurement_time_hr, which_max_gc(od_cor)),
            doub_time = doubling_time(y = max_percap),
            auc = auc(y = od_cor, x = measurement_time_hr),
            .groups = 'drop')

# run a check to see if the number of rows is as expected
if (nrow(d_sum) == n_expected) {
  message("✅ Dataframe has the expected number of rows: ", n_expected)
} else {
  message("⚠️ Warning: Expected ", n_expected, " rows, but found ", nrow(d_sum))
}


# do rolling regression using custom function
roll_regress <- function(x){
  temp <- data.frame(x)
  mod <- lm(temp, na.action = na.omit)
  temp <- data.frame(slope = coef(mod)[[2]],
                     slope_lwr = confint(mod)[2, ][[1]],
                     slope_upr = confint(mod)[2, ][[2]],
                     intercept = coef(mod)[[1]],
                     rsq = summary(mod)$r.squared, stringsAsFactors = FALSE)
  return(temp)
}

models <- d_filt %>%
  mutate(ln_od = log(od_cor),
         ln_od = ifelse(is.infinite(ln_od), NA, ln_od)) %>%
  filter(od_cor > 0.001) %>%
  group_by(across(all_of(groupings))) %>%
  do(cbind(model = select(., ln_od, measurement_time_hr) %>% 
             zoo::rollapplyr(width = 7, roll_regress, by.column = FALSE, fill = NA, align = 'center'),
           measurement_time_hr = select(., measurement_time_hr),
           ln_od = select(., ln_od))) %>%
  rename_all(., gsub, pattern = 'model.', replacement = '')

growth_rates <- models %>%
  filter(slope == max(slope, na.rm = TRUE)) %>%
  ungroup()

# do rolling regression using growthrates

# get_par_df
growthrates_get_coef <- function(x){
  return(coef(x) %>% t() %>%  as_tibble())
}

models <- d_filt %>%
  filter(od_cor > 0.001) %>%
  group_by(across(all_of(groupings))) %>%
  nest() %>%
  mutate(models = map(data, ~growthrates::fit_easylinear(.x$measurement_time_hr, .x$od_cor, h = 7, quota = 1)))

models <- models %>%
  mutate(output = map(models, growthrates_get_coef)) %>%
  unnest(output)

#--------------------#
# compare methods ####
#--------------------#

# lets compare growth rates from each of the methods
d_compare <- left_join(select(d_sum, all_of(groupings), gcplyr_gr = max_percap, gcplyr_lag = lag_time),
                       select(growth_rates, all_of(groupings), rollregress_gr = slope)) %>%
  left_join(., select(models, all_of(groupings), growthrates_gr = mumax, growthrates_lag = lag))

# plot results

# compare growth rate
select(d_compare, ends_with('gr')) %>%
  ggpairs()

# compare lag time
ggplot(d_compare, aes(x = gcplyr_lag, y = growthrates_lag)) +
  geom_point() +
  geom_abline()

# pretty good agreement between the three methods