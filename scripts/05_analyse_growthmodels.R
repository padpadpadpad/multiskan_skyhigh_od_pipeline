# ---------------------------
# Purpose of script: fit growth models to processed data
#
# What this script does:
# 1. Read in all data
# 2. Correct for blanks
# 3. Remove empty wells
# 4. Fits growth models
# 5. Visualises growth models
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
librarian::shelf(tidyverse, nls.multstart, foreach, parallel, progress, doSNOW)

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

# set up your grouping variables to get each individual well for each plate
# this bit is interactive
groupings <- c("file", "well", "serial_no", "id")

# output name
output <- paste('04_', output, sep = '')

# window_width - set to calculate a rolling regression over the course of an hour
measurements_every_X_min = 4
window_width = 60 / measurements_every_X_min
rm(measurements_every_X_min)

#---------------------#

#-----------------#
# wrangle data ####
#-----------------#

# remove any od values below certain value - remove some noise
# generally leave this in for the logistic growth models
d_filt <- filter(d, od_cor > 0.005)

# quick plot of data ####
d_filt %>%
  ggplot(aes(
    x = measurement_time_hr,
    y = od_cor,
    group = interaction(well)
  )) +
  geom_line(alpha = 0.3) +
  facet_wrap(~file) +
  theme_bw()

# calculate the expected number of rows in the output - 1 per group
n_expected <- select(d_filt, all_of(groupings)) %>%
  distinct() %>%
  nrow()
n_expected


#----------------------------#
# calculate model metrics ####
#----------------------------#

# fit a bunch of parametric non-linear models using equations from growthrates, but in nls.multstart
# info on fitting models can be found here: https://tpetzoldt.github.io/growthrates/doc/Introduction.html
# GitHub link is here: https://github.com/tpetzoldt/growthrates

# write functions for the models

# baranyi
# ?grow_baranyi for explanation of parameters
baranyi <- function(time, n0, k, mumax, h0) {
  temp <- time +
    1 / mumax * log(exp(-mumax * time) + exp(-h0) - exp(-mumax * time - h0))
  return(exp(
    log(n0) +
      mumax * temp -
      log(1 + (exp(mumax * temp) - 1) / exp(log(k) - log(n0)))
  ))
}

# logistic
# ?grow_logistic for explanation of parameters
logistic <- function(time, n0, k, mumax) {
  return((k * n0) / (n0 + (k - n0) * exp(-mumax * time)))
}

# gompertz with no lag
# ?grow_gompertz2 for explanation of parameters
gompertz_nolag <- function(time, n0, k, mumax) {
  return(n0 * (k / n0)^(exp(-exp((-exp(1) * mumax * time) / log(k / n0) + 1))))
}

# gompertz with lag
# ?grow_gompertz2 for explanation of parameters
gompertz_lag <- function(time, n0, k, mumax, lag) {
  return(
    n0 * (k / n0)^(exp(-exp((exp(1) * mumax * (lag - time)) / log(k / n0) + 1)))
  )
}

# huang model
# ?grow_huang for explanation of parameters
huang <- function(time, n0, k, mumax, lag, alpha) {
  temp <- time +
    1 / alpha * log((1 + exp(-alpha * (time - lag))) / (1 + exp(alpha * lag)))
  log_y <- log(n0) + log(k) - log(n0 + (k - n0) * exp(-mumax * temp))
  return(exp(log_y))
}

# richards model
# ?grow_richard for explanation of parameters
richards <- function(time, n0, mumax, k, beta) {
  return(k * (1 - exp(-beta * mumax * time) * (1 - (n0 / k)^-beta))^(-1 / beta))
}

# add unique id to loop through for models
d_filt3 <- unite(d_filt, id, all_of(groupings), sep = ':', remove = FALSE) %>%
  select(
    id,
    file,
    run,
    temp,
    well,
    serial_no,
    measurement_time_hr,
    od_cor,
    plate_well_id
  )

ids <- unique(d_filt3$id)

# fit models using foreach
# DO NOT FIT GOMPERTZ NO LAG

# set timer running
runstart <- lubridate::now()

# set up cluster
cl <- parallel::makeCluster(parallel::detectCores() - 1)
#doParallel::registerDoParallel(cl)
registerDoSNOW(cl)

# set up progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent | :elapsed | eta: :eta",
  total = length(ids),
  width = 80
)

progress <- function(n) {
  pb$tick(tokens = list(id = ids[n]))
}

opts <- list(progress = progress)

# run foreach model
fits <- foreach(
  i = 1:length(ids),
  .options.snow = opts,
  .combine = bind_rows,
  .packages = c("dplyr", "nls.multstart")
) %dopar%
  {
    temp_id <- ids[i]

    # filter for specific dataset
    temp <- dplyr::filter(d_filt3, id == temp_id) %>%
      dplyr::arrange(measurement_time_hr)

    # set start vals and fit logistic model
    start_params <- c(n0 = min(temp$od_cor), k = max(temp$od_cor), mumax = 0.5)
    lower <- c(n0 = 0, k = 0, mumax = 0)
    upper <- c(
      n0 = min(temp$od_cor[temp$od_cor > 0]),
      k = max(temp$od_cor) * 10,
      mumax = Inf
    )

    mod1 <- tryCatch(
      {
        nls_multstart(
          od_cor ~ logistic(time = measurement_time_hr, n0, k, mumax),
          data = temp,
          start_lower = start_params * 0.5,
          start_upper = start_params * 1.5,
          lower = lower,
          upper = upper,
          iter = 1000,
          supp_errors = 'Y',
          lhstype = 'random',
          na.action = na.omit
        )
      },
      error = function(e) NULL
    )

    # set start values and fit baranyi
    start_params <- c(
      n0 = min(temp$od_cor),
      k = max(temp$od_cor),
      mumax = 0.5,
      h = 0.5
    )
    lower <- c(n0 = 0, k = 0, mumax = 0, h = 0)
    upper <- c(
      n0 = min(temp$od_cor[temp$od_cor > 0]),
      k = max(temp$od_cor) * 10,
      mumax = Inf,
      h = Inf
    )

    mod2 <- tryCatch(
      {
        nls_multstart(
          od_cor ~ baranyi(time = measurement_time_hr, n0, k, mumax, h0),
          data = temp,
          start_lower = start_params * 0.5,
          start_upper = start_params * 1.5,
          lower = lower,
          upper = upper,
          iter = 1000,
          supp_errors = 'Y',
          lhstype = 'random',
          na.action = na.omit
        )
      },
      error = function(e) NULL
    )

    # set start values and fit gompertz lag
    start_params <- c(
      n0 = min(temp$od_cor),
      k = max(temp$od_cor),
      mumax = 0.5,
      lag = 0.5
    )
    lower <- c(n0 = 0, k = 0, mumax = 0, lag = 0)
    upper <- c(
      n0 = min(temp$od_cor[temp$od_cor > 0]),
      k = max(temp$od_cor) * 10,
      mumax = Inf,
      lag = Inf
    )

    mod3 <- nls_multstart(
      od_cor ~ gompertz_lag(time = measurement_time_hr, n0, k, mumax, lag),
      data = temp,
      start_lower = start_params * 0.5,
      start_upper = start_params * 1.5,
      lower = lower,
      upper = upper,
      iter = 1000,
      supp_errors = 'Y',
      lhstype = 'random',
      na.action = na.omit
    )

    # set start values and fit huang
    # start values
    # fmt: skip
    start_params <- c(n0 = min(temp$od_cor), k = max(temp$od_cor), mumax = 0.5, lag = 0.5,alpha = 0.5)
    lower <- c(n0 = 0, k = 0, mumax = 0, lag = 0, alpha = -Inf)
    # fmt: skip
    upper <- c(n0 = min(temp$od_cor[temp$od_cor > 0]), k = max(temp$od_cor) * 10, mumax = Inf, lag = Inf, alpha = Inf)

    mod4 <- nls_multstart(
      od_cor ~ huang(time = measurement_time_hr, n0, k, mumax, lag, alpha),
      data = temp,
      start_lower = start_params * 0.5,
      start_upper = start_params * 1.5,
      lower = lower,
      upper = upper,
      iter = 1000,
      supp_errors = 'y',
      lhstype = 'random',
      na.action = na.omit
    )

    # fit and set start params richards
    # fmt: skip
    start_params <- c(n0 = min(temp$od_cor), k = max(temp$od_cor), mumax = 0.5, beta = 0.5)
    lower <- c(n0 = 0, k = 0, mumax = 0, beta = -Inf)
    # fmt: skip
    upper <- c(n0 = min(temp$od_cor[temp$od_cor > 0]), k = max(temp$od_cor) * 10, mumax = Inf, beta = Inf)

    mod5 <- nls_multstart(
      od_cor ~ richards(time = measurement_time_hr, n0, mumax, k, beta),
      data = temp,
      start_lower = start_params * 0.5,
      start_upper = start_params * 1.5,
      lower = lower,
      upper = upper,
      iter = 1000,
      supp_errors = 'Y',
      lhstype = 'random',
      na.action = na.omit
    )

    # set up for loop for all the fits
    tibble(
      id = temp_id,
      fit_logistic = list(mod1),
      fit_baranyi = list(mod2),
      fit_gompertz = list(mod3),
      fit_huang = list(mod4),
      fit_richards = list(mod5)
    )
  }

stopCluster(cl)

time1 <- lubridate::as.duration(lubridate::now() - runstart)
time1

# save this out
saveRDS(fits, file.path('data/metrics', output))

#---------------------#
# Visualise models ####
#---------------------#

# make a plot to visualise every plate
runs <- unique(d_filt$run)

fits <- separate(fits, id, into = groupings, sep = ':', remove = FALSE)

# open a pdf
pdf(
  file.path(
    'plots/metric_plots',
    paste(tools::file_path_sans_ext(output), '.pdf', sep = '')
  ),
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

    temp_models <- filter(fits, run == runs[i] & temp == temp_temps[j]) %>%
      mutate(column = str_extract(well, "[A-Z]+"), row = parse_number(well)) %>%
      pivot_longer(
        cols = c(
          fit_logistic,
          fit_baranyi,
          fit_gompertz,
          fit_huang,
          fit_richards
        ),
        names_to = 'model',
        values_to = 'fit'
      ) %>%
      mutate(preds = map(fit, broom::augment)) %>%
      unnest(preds)

    temp_id <- select(temp_od2, column, row, id) %>%
      distinct() %>%
      mutate(x = 0, y = 1)

    # fmt: skip
    temp_plot <- temp_od2 %>%
      ggplot(aes(x = measurement_time_hr, y = od_cor)) +
      geom_line(linewidth = 2, col = 'grey') +
      geom_line(aes(x = measurement_time_hr, y = .fitted, col = model), temp_models, linewidth = 0.75) +
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
