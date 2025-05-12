# ---------------------------
# Purpose of script: Calculate growth rate metrics from processed data
#
# What this script does:
# 1. Read in processed data from 03_process_od.R
# 2. Remove noise at low ODs
# 3. Calculate derivative for each well
# 4. Calculate growth metrics and save them out
# 5. Visualise growth metrics on growth curves
#
# Author: Dr. Daniel Padfield
#
# Date Created:  2025-02-20
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
if (!requireNamespace("librarian", quietly = TRUE)) {
  install.packages("librarian")
}

# if BiocManager is not installed, install it
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# if Biobase is not installed, install it from Bioconductor
if (!requireNamespace("Biobase", quietly = TRUE)) {
  BiocManager::install("Biobase")
}

# load packages
# fmt: skip
librarian::shelf(tidyverse, growthrates, gcplyr)

## ---------------------------

#---------------------#
# things to change ####
#---------------------#

# name of files to read in
file_names <- list.files('data/processed', full.names = TRUE)

# read in data
d <- file_names %>% map_df(read.csv)

# split file identifier if needed

# check the number of unique combinations there should be
select(d, file, well, id, serial_no) %>%
  distinct() %>%
  nrow()

# set up your grouping variables to get each individual well for each plate
# this bit is interactive
groupings <- c("file", "well", "serial_no", "id")

# set blank method: options are: well_specific or blank_median
blank_method <- "blank_median"

# output name
output <- 'all_growth_metrics.csv'

# window_width - set to calculate a rolling regression over the course of an hour
measurements_every_X_min = 4
window_width = 60 / measurements_every_X_min
rm(measurements_every_X_min)

#---------------------#

#---------------------#
# calculate blanks ####
#---------------------#

# method 1: create a well-specific blank
if (blank_method == "well_specific") {
  # take the first three time points
  # average per grouping
  d_blank <- d %>%
    group_by(across(all_of(groupings))) %>%
    filter(measurement_time_hr %in% sort(measurement_time_hr)[1:3]) %>%
    summarise(
      ave_blank = mean(raw_absorbance),
      sd_blank = sd(raw_absorbance),
      .groups = 'drop'
    )
}

# method 2 Calculate blank from average of all blank wells in each plate
# use the median
# assumes all blanks start with X
if (blank_method == "blank_median") {
  d_blank <- d %>%
    filter(!well %in% outside_wells) %>%
    filter(substr(id, 1, 1) == "X") %>%
    group_by(across(all_of(groupings[!groupings %in% c('well', 'id')]))) %>%
    summarise(
      ave_blank = median(raw_absorbance),
      sd_blank = sd(raw_absorbance),
      .groups = 'drop'
    )
}

# visualise average blank
hist(d_blank$ave_blank)

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

# remove any od values below certain value - remove some noise
d_filt2 <- filter(d_filt, od_cor > 0.005)

# quick plot of data ####
d_filt %>%
  ggplot(aes(
    x = measurement_time_hr,
    y = od_cor,
    group = interaction(well, run)
  )) +
  geom_line(alpha = 0.3) +
  facet_wrap(~temp) +
  theme_bw()

#--------------------------------------------------------------#
# remove samples because of biologically implausible signal ####
#--------------------------------------------------------------#

# from the plots in first_look_plots, you can remove some wells that have biologically implausible signal (e.g. massive spikes) that are not going to be easy to model and will give poor estimates

# create unique ID of file and well
d_filt2 <- d_filt2 %>%
  mutate(plate_well_id = paste(file, well, sep = '_'))

to_remove <- c()

# remove these wells
d_filt2 <- d_filt2 %>%
  filter(!plate_well_id %in% to_remove)

# calculate the expected number of rows in the output - 1 per group
n_expected <- select(d_filt2, all_of(groupings)) %>%
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
d_filt2 <- d_filt2 %>%
  group_by(across(all_of(groupings))) %>%
  mutate(
    deriv_percap = calc_deriv(
      x = measurement_time_hr,
      y = od_cor,
      percapita = TRUE,
      blank = 0,
      window_width_n = window_width,
      trans_y = "log"
    )
  ) %>%
  ungroup()

# calculate growth metrics
d_sum <- d_filt2 %>%
  group_by(across(all_of(groupings))) %>%
  summarise(
    min_dens = first_minima(od_cor, return = 'y'),
    lag_time = lag_time(
      y = od_cor,
      x = measurement_time_hr,
      deriv = deriv_percap,
      blank = 0,
      y0 = min_dens
    ),
    gr = max_gc(deriv_percap, na.rm = TRUE),
    gr_time = extr_val(measurement_time_hr, which_max_gc(deriv_percap)),
    gr_dens = extr_val(od_cor, which_max_gc(deriv_percap)),
    max_dens = max_gc(od_cor, na.rm = TRUE),
    max_time = extr_val(measurement_time_hr, which_max_gc(od_cor)),
    doub_time = doubling_time(y = gr),
    auc = auc(y = od_cor, x = measurement_time_hr),
    .groups = 'drop'
  )

# run a check to see if the number of rows is as expected
if (nrow(d_sum) == n_expected) {
  message("✅ Dataframe has the expected number of rows: ", n_expected)
} else {
  message("⚠️ Warning: Expected ", n_expected, " rows, but found ", nrow(d_sum))
}

# save this out ####
write.csv(d_sum, file.path('data/metrics', output), row.names = FALSE)

#----------------------#
# Visualise metrics ####
#----------------------#

# make a plot to visualise every plate
runs <- unique(d_filt$run)

# open a pdf
pdf(
  file.path('plots', paste('check_gcplyr_log', '.pdf', sep = '')),
  width = 10,
  height = 6.5
)

for (i in 1:length(runs)) {
  temp_od <- filter(d_filt, run == runs[i]) |>
    mutate(column = str_extract(well, "[A-Z]+"), row = parse_number(well))

  temp_temps <- unique(temp_od$temp)

  # create a plot for each temperature
  for (j in 1:length(temp_temps)) {
    temp_od2 <- temp_od |>
      filter(temp == temp_temps[j])

    temp_params <- filter(d_sum, run == runs[i] & temp == temp_temps[j]) %>%
      mutate(column = str_extract(well, "[A-Z]+"), row = parse_number(well))

    temp_id <- select(temp_od2, column, row, id) %>%
      distinct() %>%
      mutate(x = 0, y = 1)

    # fmt: skip
    temp_plot <- temp_od2 %>%
        ggplot(aes(x = measurement_time_hr, y = log(od_cor))) +
        geom_line() +
        geom_vline(aes(xintercept = lag_time), temp_params, col = 'red') +
        geom_point(aes(gr_time, log(gr_dens)), temp_params, col = 'red') +
        geom_point(aes(max_time, log(max_dens)), temp_params, col = 'blue') +
        geom_label(aes(x = 0, y = 1.25, label = id), temp_id, hjust = 0, vjust = 0.8, label.size = NA) +
        facet_grid(column~row, switch = 'y') +
        theme_bw() +
        labs(title = paste('Run:', runs[i], ';Temp:', temp_temps[j], 'ºC', sep = ' '),
             x = 'Time (hr)',
             y = 'log Absorbance (OD600)') +
        ylim(c(min(log(temp_od2$od_cor)[is.finite(log(temp_od2$od_cor))]), 1.25)) +
        NULL

    print(temp_plot)
  }
}

dev.off()

# do the same plot but unlogged
# open a pdf
pdf(
  file.path('plots', paste('check_gcplyr_nolog', '.pdf', sep = '')),
  width = 10,
  height = 6.5
)

for (i in 1:length(runs)) {
  temp_od <- filter(d_filt, run == runs[i]) |>
    mutate(column = str_extract(well, "[A-Z]+"), row = parse_number(well))

  temp_temps <- unique(temp_od$temp)

  # create a plot for each temperature
  for (j in 1:length(temp_temps)) {
    temp_od2 <- temp_od |>
      filter(temp == temp_temps[j])

    temp_params <- filter(d_sum, run == runs[i] & temp == temp_temps[j]) %>%
      mutate(column = str_extract(well, "[A-Z]+"), row = parse_number(well))

    temp_id <- select(temp_od2, column, row, id) %>%
      distinct() %>%
      mutate(x = 0, y = 1)

    # fmt: skip
    temp_plot <- temp_od2 %>%
      ggplot(aes(x = measurement_time_hr, y = od_cor)) +
      geom_line() +
      geom_vline(aes(xintercept = lag_time), temp_params, col = 'red') +
      geom_point(aes(gr_time, gr_dens), temp_params, col = 'red') +
      geom_point(aes(max_time, max_dens), temp_params, col = 'blue') +
      geom_label(aes(x = 0, y = 1.25, label = id), temp_id, hjust = 0, vjust = 0.8, label.size = NA) +
      facet_grid(column~row, switch = 'y') +
      theme_bw() +
      labs(title = paste('Run:', runs[i], ';Temp:', temp_temps[j], 'ºC', sep = ' '),
           x = 'Time (hr)',
           y = 'Absorbance (OD600)') +
      ylim(c(min(temp_od2$od_cor[is.finite(temp_od2$od_cor)]), 1.25)) +
      NULL

    print(temp_plot)
  }
}

dev.off()
