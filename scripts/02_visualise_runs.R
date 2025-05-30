# ---------------------------
# Purpose of script: Create a visualisations for raw OD readings from the Multiskan Skyhigh
#
# What this script does:
# 1. Loads in data
# 2. Creates a plot for each file, with all wells combined
#
# Author: Dr. Daniel Padfield
#
# Date Created:  2025-01-06
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
librarian::shelf(tidyverse)

## ---------------------------

#---------------------#
# things to change ####
#---------------------#

# set run, this will be used to label the plots
# should be the same as output from 01_process_od.R
input <- 'output.csv'
input <- paste('01_', input, sep = '')

input_no_ext <- tools::file_path_sans_ext(input)

# change this as needed
d_od <- read.csv(file.path("data/processed", input))

#---------------------#

# plot each well across each combination of factors ####

# create an ID column for each plate
d_od <- d_od %>%
  mutate(plate_id = paste(file, serial_no, sep = '_'))

plate_ids <- unique(d_od$plate_id)

# open a pdf
pdf(
  file.path(
    'plots/first_look_plots',
    paste(input_no_ext, '.pdf', sep = '')
  ),
  width = 10,
  height = 6.5
)

for (i in 1:length(plate_ids)) {
  temp_od <- d_od |>
    filter(plate_id == plate_ids[i]) |>
    mutate(column = str_extract(well, "[A-Z]+"), row = parse_number(well))

  temp_id <- select(temp_od, column, row, id) %>%
    distinct() %>%
    mutate(x = 0, y = 1)

  temp_plot <- temp_od %>%
    ggplot(aes(x = measurement_time_hr, y = raw_absorbance)) +
    geom_line() +
    geom_text(aes(x = 0, y = 1, label = id), temp_id, hjust = 0) +
    facet_grid(column ~ row, switch = 'y') +
    theme_bw() +
    labs(
      title = paste('Plate ID:', plate_ids[i]),
      x = 'Time (hr)',
      y = 'Absorbance (OD600)'
    ) +
    ylim(c(0, 1.25))

  print(temp_plot)
}

dev.off()

# remove temporary items in plot
rm(list = c('temp_od', 'temp_plot'))
