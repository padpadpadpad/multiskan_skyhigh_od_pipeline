# ---------------------------
# Purpose of script: Process files from Multiskan Skyhigh
#
# What this script does:
# 1. Lists files ready to be merged, this could be changed to improve file selection (e.g. only do "run1")
# 2. Calculates measurement time to add to the second files
# 3. Removes outer wells
# 4. Adds in correct master plate
# 5. Save out
#
# Author: Dr. Daniel Padfield
#
# Date Created:  2024-12-16
#
# Copyright (c) Daniel Padfield, 2024
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
librarian::shelf(tidyverse, readxl)

## ---------------------------

# load in extra functions
source('scripts/00_functions.R')

#---------------------#
# things to change ####
#---------------------#

# list files in
# change this to only be the correct files
files <- list.files('data/raw', full.names = TRUE, pattern = 'xlsx')

# set output name
output <- 'output.csv'

#-----------------#
# read in data ####
#-----------------#

# remove file if there is a $ in the name - means they are temporary
files <- files[!grepl('\\$', files)]

# bind together all the OD data using get_od() - can see how it works in 00_functions.R
d_od <- map_df(files, get_od) %>%
  mutate(well = gsub(" ", "", well))

# bind all info data using get_info() - can see how it works in 00_functions.R
d_info <- map_df(files, get_info)

# bind all the runlog data - using get_runlog() - can see how it works in 00_functions.R
d_runlog <- map_df(files, possibly(get_runlog))
# check which files are not in d_runlog
no_runlog <- files[
  !tools::file_path_sans_ext(basename(files)) %in% unique(d_runlog$file)
]

# read in master plates
# d_master_plates
d_layout <- map_df(files, get_plate_layout)

#--------------------#
# wrangle od data ####
#--------------------#

# bind well data with od data
d_od <- left_join(d_od, d_layout)

# grab serial number from the info file
d_serial <- filter(d_info, str_detect(param, 'Serial number')) %>%
  select(file, serial_no = value)

# bind with od
d_od <- left_join(d_od, d_serial)

# split the file name of d_od to its useful parts
d_od <- d_od %>%
  select(
    file,
    serial_no,
    measurement_time_s,
    wavelength_nm,
    well,
    raw_absorbance,
    id
  )

# convert time to minutes and hours
d_od <- mutate(
  d_od,
  measurement_time_min = measurement_time_s / 60,
  measurement_time_hr = measurement_time_min / 60
)

# save out data
write.csv(d_od, file.path('data/processed', output), row.names = FALSE)
