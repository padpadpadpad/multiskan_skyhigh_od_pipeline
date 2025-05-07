# ---------------------------
# Purpose of script: process OD data
#
# What this script does:
# 1. Read in all data
# 2. Correct for blanks
# 3. Remove blanks
# 4. Filter out biologically implausible data
# 5. Save out file
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
librarian::shelf(tidyverse)

## ---------------------------

#---------------------#
# things to change ####
#---------------------#

# name of files to read in
file_names <- list.files('data/processed', full.names = TRUE)

# read in data
d <- file_names %>% map_df(read.csv)

# check the number of unique combinations there should be
select(d, file, well, id, serial_no) %>%
  distinct() %>%
  nrow()

# split file identifier
d <- separate(d, file, into = c('run', 'temp'), sep = '_', remove = FALSE)

# set up your grouping variables to get each individual well for each plate
# this bit is interactive
groupings <- c("file", "run", "temp", "well", "serial_no", "id")

# set blank method: options are: well_specific or blank_median
blank_method <- "blank_median"

# output name
output <- 'all_growth_processed.csv'

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

# name outside wells
# fmt: skip
outside_wells <- c('A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12',
                   'B1', 'B12', 'C1', 'C12', 'D1', 'D12', 'E1', 'E12',
                   'F1', 'F12', 'G1', 'G12', 'H1', 'H2', 'H3', 'H4',
                   'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11')

# do blank subtraction
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

#--------------------------------------------------------------#
# remove samples because of biologically implausible signal ####
#--------------------------------------------------------------#

# from the plots in first_look_plots, you can remove some wells that have biologically implausible signal (e.g. massive spikes) that are not going to be easy to model and will give poor estimates

# create unique ID of file and well
d_filt <- d_filt %>%
  mutate(plate_well_id = paste(file, well, sep = '_'))

to_remove <- c()

# remove these wells
d_filt <- d_filt %>%
  filter(!plate_well_id %in% to_remove)

#------------------#
# save out data ####
#------------------#

write.csv(d_filt, file.path('data/processed', output))
