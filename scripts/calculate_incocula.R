# ---------------------------
# Purpose of script: To calculate volumes of media and inoculum to pipette into wells for an inoculum plate.
#
# What this script does:
# 1. Takes in readings of OD.
# 2. Calculates volume of media and inoculum needed for each well
# 3. Saves this out as a csv
#
# Author: Dr. Daniel Padfield
#
# Date Created:  2024-12-05
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
librarian::shelf(tidyverse, openxlsx)

## ---------------------------

#---------------------#
# custom functions ####
#---------------------#

# write function to try all of the culture volumes and return the answer which has the closest total volume to the target
# it will pick the volumes that lie within the low and high volume targets, and if not return the one closest to the mean
get_culture_volume <- function(od, target_od, culture_volume, target_total_volume_low, target_total_volume_high){
  sapply(od, function(current_od) {
    
    # Calculate total volume for each culture volume
    total_volume <- round((current_od * culture_volume) / target_od, -1)
    
    # Find the culture volumes that are within the target total volume range
    to_keep <- which(total_volume >= target_total_volume_low & total_volume <= target_total_volume_high)
    
    # find the value that is closest to the high target
    closest <- to_keep[which.min(abs(total_volume[to_keep] - target_total_volume_high))]
    
    # if none are within range, find the volume closest to the middle of the low and high target
    if(length(to_keep) == 0){
      closest <- which.min(abs(total_volume - mean(c(target_total_volume_low, target_total_volume_high))))
      }
    
    # Return the corresponding culture volume
    return(culture_volume[closest])
  })
}

# 1. Input your readings of OD ####
d_OD <- rnorm(15, mean = 1.3, sd = 0.2)

# 2. give clones a name, probably what they are on your freezer stock, or the label for the plate
d_clone <- paste('label', 1:length(d_OD), sep = '_')

# 3. combine these into a dataframe
d <- data.frame(clone = d_clone,
                od = d_OD)

d

# 4. calculate concentrations ####

# set target OD
target_od <- 0.1

# set amount of culture to transfer
# this can be changed as you get a better feel for your overnight ODs
# we will calculate multiple and pick the one closest to the desired volume
culture_volumes <- seq(0, 280, by = 5)
culture_volumes

# are any of your OD readings below the target od?
dplyr::filter(d, od < target_od)
d_low <- dplyr::filter(d, od < target_od)
# option is to resuspend these in a lower volume to concentrate, but not ideal

d <- dplyr::filter(d, !clone %in% d_low$clone)

# use the following equation for calculating stock solution volumes:
# original OD volume = (target OD x total volume of target OD) / original OD concentration
# rearrange this as we will always put in 100µl of the original OD volume
# total volume of target OD = (original OD concentration x original OD volume) / target OD
# you can read more about the equation here: https://www.thoughtco.com/dilutions-from-stock-solutions-606085
d2 <- dplyr::mutate(d,
                   # calculate total volume of target OD
                   culture_volume = get_culture_volume(od, target_od, culture_volumes, target_total_volume_low = 240, target_total_volume_high = 290),
                   # can calculate amount of media needed
                   total_volume = round((od * culture_volume)/target_od, -1),
                   # if any volumes are below 240, make them 240, if any are above 300, make them 300
                   total_volume2 = case_when(total_volume < 240 ~ 240,
                                             total_volume > 300 ~ 300,
                                             TRUE ~ total_volume),
                   media_volume = total_volume2 - culture_volume,
                   target_od = target_od)

d2

# look at any times where total_volume != total_volume2
dplyr::filter(d2, total_volume != total_volume2)

# look at the unique number of different culture volumes
unique(d2$culture_volume)

d2 <- select(d2, -total_volume)

# if we're putting in 20 µl into 5 plates for each well, we need at least 120µl in each well, so need total volume to be at least 240µL.

# save this out as a csv somewhere
folder <- 'FOLDER PATH HERE'
write.xlsx(d2, file.path('~/Desktop', 'inoculum_plate1.xlsx'), rowNames = FALSE)