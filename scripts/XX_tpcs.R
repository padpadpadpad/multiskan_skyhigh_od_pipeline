# ---------------------------
# Purpose of script: Analyse the thermal performance curves
#
# What this script does:
# 1.
# 2.
# 3.
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
librarian::shelf(tidyverse, rTPC, nls.multstart, broom, ggbeeswarm, lme4, brms, padpadpadpad/dataViewer)

## ---------------------------

# read in growth metrics
d_growth <- read.csv("data/metric_data/all_growth_metrics.csv")

# read in data about the different strains
d_isolates <- read.csv('data/well_id.csv')

d_isolates %>% group_by(biome, habitat) %>%
  tally()

# bind together
d_growth <- left_join(d_growth, d_isolates, by = c('id' = 'id'))

# remove any curves that have less than 7 temperatures
d_growth <- group_by(d_growth, id) %>%
  filter(length(unique(temp)) > 7) %>%
  ungroup()

unique(d_growth$id)

# time to see whether growth rate and auc are temperature dependent

# facet by ID
ggplot(d_growth, aes(temp, auc)) +
  geom_point() +
  facet_wrap(~id, scales = 'free_y') +
  theme_bw()

ggplot(d_growth, aes(temp, max_percap)) +
  geom_point(aes(text = run)) +
  facet_wrap(~id, scales = 'free_y') +
  theme_bw()

# growth rate is quite messy

#-----------------------------------#
# fit thermal performance curves ####
#-----------------------------------#

get_stan_code <- "
  real sharpeschoolhigh_1981(real temp, real rtref, real e, real eh, real th, real tref) {
    
    real pred;
    real tref2;
    real k;
    real boltzmann;
    real inactivation;
    
    tref2 = tref + 273.15;
    k = 8.62e-05;
    
    boltzmann = rtref * exp(e/k * (1/tref2 - 1/(temp + 
        273.15)));
    
    inactivation = 1/(1 + exp(eh/k * (1/(th + 273.15) - 
        1/(temp + 273.15))));
        
    pred = boltzmann*inactivation;
        
    return(pred);
    }
"

get_stan_code <- function(x_var = 'temp', model_name = 'sharpeschoolhigh_1981'){
  # define sharpe schoolfield model
  if(model_name == 'sharpeschoolhigh_1981'){
    return(gsub('temp', x_var, "
    real sharpeschoolhigh_1981(real temp, real rtref, real e, real eh, real th, real tref) {
    
    real pred;
    real tref2;
    real k;
    real boltzmann;
    real inactivation;
    
    tref2 = tref + 273.15;
    k = 8.62* 10^-5;
    
    boltzmann = rtref * exp(e/k * (1/tref2 - 1/(temp + 
        273.15)));
    
    inactivation = 1/(1 + exp(eh/k * (1/(th + 273.15) - 
        1/(temp + 273.15))));
        
    pred = boltzmann*inactivation;
        
    return(pred);
    }
"))
  }
}

# cheeky bayesian model attempt
# fit bacteria model using brms
brms_form <- bf(auc ~ rtref * exp(e/(8.62* 10^-5) * (1/(273.15+18) - 1/(temp + 273.15))) * 1/(1 + exp(eh/(8.62 * 10^-5) * (1/(th + 273.15) - 1/(temp + 273.15)))),
                rtref ~ 0 + biome + (1|habitat/id), 
                e ~ 0 + biome + (1|habitat/id), 
                eh ~ 0 + biome + (1|habitat/id), 
                th ~ 0 + biome + (1|habitat/id), 
                nl = TRUE)

d_auc <- filter(d_growth, !is.na(auc))

# see available priors
priors <- get_prior(brms_form, d_auc)

priors

priors <- filter(priors, class == 'b' & coef == '')
priors

# use rTPC helper functions
start_vals <- get_start_vals(d_auc$temp, d_auc$auc, model_name = 'sharpeschoolhigh_1981')
lower_lims <- get_lower_lims(d_auc$temp, d_auc$auc, model_name = 'sharpeschoolhigh_1981')
upper_lims <- get_upper_lims(d_auc$temp, d_auc$auc, model_name = 'sharpeschoolhigh_1981')

# set priors
priors$lb <- lower_lims[sort(names(lower_lims))]
priors$ub <- upper_lims[sort(names(upper_lims))]
priors$prior <- paste('normal(', start_vals[sort(names(start_vals))], ',', 10, ')', sep = '')

priors

fit_brms <- brm(brms_form,
                data = d_auc, 
                prior = priors,
                family = gaussian(),
                chains = 1,
                cores = 1,
                iter = 4000,
                control = list(adapt_delta = 0.99, max_treedepth = 30))

fit_bayes

# fit sharpe-schoolfield to both max_percap (growth rate) and auc (biomass accumulation)
# the possibly argument returns NA (where fits have not worked) to occur
# need to remove NAs
d_fits_auc <- d_growth %>%
  filter(!is.na(auc)) %>%
  group_by(id, habitat, biome) %>%
  nest() %>%
  mutate(auc_flex = map(data, possibly(~nls_multstart(auc ~ flextpc_2024(temp, tmin, tmax, rmax, topt, beta),
                                                      data = .,
                                                      iter = 1000,
                                                      start_lower = get_start_vals(.x$temp, .x$auc, model_name = 'flextpc_2024') - 10,
                                                      start_upper = get_start_vals(.x$temp, .x$auc, model_name = 'flextpc_2024') + 10,
                                                      lower = get_lower_lims(.x$temp, .x$auc, model_name = 'flextpc_2024'),
                                                      upper = get_upper_lims(.x$temp, .x$auc, model_name = 'flextpc_2024'),
                                                      supp_errors = 'Y',
                                                      convergence_count = FALSE,
                                                      lhstype = 'random'), NA)),
    auc_sharpeschool = map(data, possibly(~nls_multstart(auc ~ sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 18),
                                                         data = .,
                                                         iter = 1000,
                                                         start_lower = get_start_vals(.x$temp, .x$auc, model_name = 'sharpeschoolhigh_1981') - 10,
                                                         start_upper = get_start_vals(.x$temp, .x$auc, model_name = 'sharpeschoolhigh_1981') + 10,
                                                         lower = get_lower_lims(.x$temp, .x$auc, model_name = 'sharpeschoolhigh_1981'),
                                                         upper = get_upper_lims(.x$temp, .x$auc, model_name = 'sharpeschoolhigh_1981'),
                                                         supp_errors = 'Y',
                                                         convergence_count = FALSE,
                                                         lhstype = 'random'), NA)),
    auc_gaussian = map(data, possibly(~nls_multstart(auc ~ gaussian_1987(temp, rmax, topt, a),
                                                    data = .,
                                                    iter = 1000,
                                                    start_lower = get_start_vals(.x$temp, .x$auc, model_name = 'gaussian_1987') * 0.5,
                                                    start_upper = get_start_vals(.x$temp, .x$auc, model_name = 'gaussian_1987') * 1.5,
                                                    lower = get_lower_lims(.x$temp, .x$auc, model_name = 'gaussian_1987'),
                                                    upper = get_upper_lims(.x$temp, .x$auc, model_name = 'gaussian_1987'),
                                                    supp_errors = 'Y',
                                                    convergence_count = FALSE,
                                                    lhstype = 'random'), NA)),
    auc_thomas = map(data, possibly(~nls_multstart(auc ~ thomas_2017(temp, a,b,c,d,e),
                                                   data = .,
                                                   iter = 1000,
                                                   start_lower = get_start_vals(.x$temp, .x$auc, model_name = 'thomas_2017') * 0.5,
                                                   start_upper = get_start_vals(.x$temp, .x$auc, model_name = 'thomas_2017') * 1.5,
                                                   lower = get_lower_lims(.x$temp, .x$auc, model_name = 'thomas_2017'),
                                                   upper = get_upper_lims(.x$temp, .x$auc, model_name = 'thomas_2017'),
                                                   supp_errors = 'Y',
                                                   convergence_count = FALSE,
                                                   lhstype = 'random'), NA)))

# lets look at the predicted values of the fits
d_preds_auc <- 
  # stack models into a single column, with an id column for model_name
  pivot_longer(d_fits_auc, names_to = 'model_name', values_to = 'fit', contains('auc')) %>%
  filter(., !is.na(fit)) %>%
  mutate(new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100)))) %>%
  # get rid of original data column
  select(., -data) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(id, habitat, biome, preds, model_name) %>%
  # unlist the preds list column
  unnest(preds) %>%
  # filter out very low .fitted values
  filter(.fitted > -10)

# make plot
ggplot(d_growth, aes(temp, auc)) +
  geom_point() +
  facet_wrap(~id, scales = 'free_y') +
  theme_bw() +
  geom_line(aes(temp, .fitted, col = model_name), d_preds_auc) +
  facet_wrap(~id, scales = 'free_y') +
  theme_bw()

# calculate AIC scores, then rank for each ID
d_aic_auc <- pivot_longer(d_fits_auc, names_to = 'model_name', values_to = 'fit', contains('auc')) %>%
  # remove NA fits
  filter(., !is.na(fit)) %>%
  # calculate AIC
  mutate(bic = map(fit, ~MuMIn::AICc(.x))) %>%
  select(id, habitat, biome, model_name, aic) %>%
  unnest(aic) %>%
  # for each ID grab only the lowest AIC
  group_by(id) %>%
  mutate(weight = MuMIn::Weights(aic)) %>%
  slice_min(n = 1, order_by = aic) %>%
  ungroup()

# see which model is best
group_by(d_aic_auc, model_name) %>%
  tally()
# model to use is sharpeschool

# calculate the parameters to test
params <- pivot_longer(d_fits_auc, names_to = 'model_name', values_to = 'fit', contains('auc')) %>%
  filter(., !is.na(fit)) %>%
  # calculate the R2 while we're at it
  mutate(r2 = map_dbl(fit, function(x){soilphysics::Rsq(x)$pseudo.R.squared}),
        params = map(fit, calc_params)) %>%
  select(id, habitat, biome, model_name, r2, params) %>%
  unnest(params) %>%
  select(-c(e, eh)) %>%
  # filter just for the sharpeschool
  filter(model_name == 'auc_flex')

# grab the estimated parameters specifically from the model - gives a better estimate of e and eh and r_tref
params2 <- pivot_longer(d_fits_auc, names_to = 'model_name', values_to = 'fit', contains('auc')) %>%
  filter(., !is.na(fit)) %>%
  filter(model_name == 'auc_flex') %>%
  mutate(params = map(fit, tidy)) %>%
  select(id, habitat, biome, model_name, params) %>%
  unnest(params) %>%
  select(., -c(std.error:p.value)) %>%
  pivot_wider(names_from = term, values_from = estimate)

params <- filter(params, model_name == 'auc_flex') %>%
  left_join(., params2)

# have an initial visualisation. This is where you need to think about which parameters you think might vary!
select(params, id, habitat, biome, model_name, rmax, topt, ctmax, thermal_safety_margin, thermal_tolerance) %>%
  pivot_longer(., names_to = 'param', values_to = 'value', rmax:thermal_tolerance) %>%
  filter(!is.infinite(value)) %>%
  ggplot(aes(biome, value, col = habitat)) +
  geom_quasirandom(size = 1, width = 0.2) +
  facet_wrap(~param, scales = 'free_y')

# simple plot to visualisa a parameter and nestedness of the data
ggplot(params, aes(biome, rmax)) +
  geom_point(aes(col = habitat), size = 1, position = position_dodge(width = 0.2)) +
  stat_summary(aes(group = habitat, col = habitat),
    fun = "mean",        #argument updated in new version.
    geom = "point",
    col = "black",
    size = 3,
    position = position_dodge(width = 0.2)
  )

ggplot(params, aes(r_tref, rmax)) +
  geom_point() +
  facet_wrap(~biome)

# you might want to save out the data here and start a new script to do the stats

#--------------------------#
# time to do some stats ####
#--------------------------#

# probably just do linear models and ignore the nested nature of the data.

# for rmax
mod1 <- lmer(rmax ~ biome + (1|habitat), params)
mod2 <- lmer(rmax ~ 1 + (1|habitat), params)
anova(mod1, mod2) # signif

# check model performance 

# 1. Look at intraclass-correlation coefficient - the proportion of the variance explained by the grouping structure in the population. Is between 0 - 1. If it is close to 0, means the observation within groups are no more similar than observations from different groups, and setting it as a random factor might not be necessary
performance::icc(mod1)
performance::check_model(mod1)
summary(mod1)

# for topt
mod1 <- lmer(topt ~ biome + (1|habitat), params)
mod2 <- lmer(topt ~ 1 + (1|habitat), params)
anova(mod1, mod2) # non signif
performance::icc(mod1)
performance::check_model(mod1)

# for ctmax
mod1 <- lmer(ctmax ~ biome + (1|habitat), params)
mod2 <- lmer(ctmax ~ 1 + (1|habitat), params)
anova(mod1, mod2)

# for thermal_safty_margin
mod1 <- lmer(thermal_safety_margin ~ biome + (1|habitat), params)
mod2 <- lmer(thermal_safety_margin ~ 1 + (1|habitat), params)
anova(mod1, mod2)

# for e
mod1 <- lmer(e ~ biome + (1|habitat), params)
mod2 <- lmer(e ~ 1 + (1|habitat), params)
anova(mod1, mod2)

# for r_tref
mod1 <- lmer(r_tref ~ biome + (1|habitat), params)
mod2 <- lmer(r_tref ~ 1 + (1|habitat), params)
anova(mod1, mod2)

# test for local adaptation
local_adaptation <- filter(d_growth_filt, temp %in% c(15,21)) %>%
  group_by(id, biome, temp) %>%
  summarise(auc = mean(auc),
            .groups = 'drop')

ggplot(local_adaptation, aes(temp,auc)) +
  geom_quasirandom(width = 0.75) +
  facet_wrap(~biome)

mod1 <- lmer(auc ~ as.character(temp) * biome + (1|id), local_adaptation)
mod2 <- lmer(auc ~ as.character(temp) + biome + (1|id), local_adaptation)
anova(mod1, mod2)
summary(mod2)

# do the same but with growth rate
d_fits_gr <- d_growth %>%
  filter(!is.na(max_percap)) %>%
  group_by(id, habitat, biome) %>%
  nest() %>%
  mutate(
    gr_sharpeschool = map(data, possibly(~nls_multstart(max_percap ~ sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 18),
                                                         data = .,
                                                         iter = 1000,
                                                         start_lower = get_start_vals(.x$temp, .x$max_percap, model_name = 'sharpeschoolhigh_1981') - 10,
                                                         start_upper = get_start_vals(.x$temp, .x$max_percap, model_name = 'sharpeschoolhigh_1981') + 10,
                                                         lower = get_lower_lims(.x$temp, .x$max_percap, model_name = 'sharpeschoolhigh_1981'),
                                                         upper = get_upper_lims(.x$temp, .x$max_percap, model_name = 'sharpeschoolhigh_1981'),
                                                         supp_errors = 'Y',
                                                         convergence_count = FALSE,
                                                         lhstype = 'random'), NA)),
    gr_gaussian = map(data, possibly(~nls_multstart(max_percap ~ gaussian_1987(temp, rmax, topt, a),
                                                     data = .,
                                                     iter = 1000,
                                                     start_lower = get_start_vals(.x$temp, .x$max_percap, model_name = 'gaussian_1987') * 0.5,
                                                     start_upper = get_start_vals(.x$temp, .x$max_percap, model_name = 'gaussian_1987') * 1.5,
                                                     lower = get_lower_lims(.x$temp, .x$max_percap, model_name = 'gaussian_1987'),
                                                     upper = get_upper_lims(.x$temp, .x$max_percap, model_name = 'gaussian_1987'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE,
                                                     lhstype = 'random'), NA)),
    gr_thomas = map(data, possibly(~nls_multstart(max_percap ~ thomas_2017(temp, a,b,c,d,e),
                                                   data = .,
                                                   iter = 1000,
                                                   start_lower = get_start_vals(.x$temp, .x$max_percap, model_name = 'thomas_2017') * 0.5,
                                                   start_upper = get_start_vals(.x$temp, .x$max_percap, model_name = 'thomas_2017') * 1.5,
                                                   lower = get_lower_lims(.x$temp, .x$max_percap, model_name = 'thomas_2017'),
                                                   upper = get_upper_lims(.x$temp, .x$max_percap, model_name = 'thomas_2017'),
                                                   supp_errors = 'Y',
                                                   convergence_count = FALSE,
                                                   lhstype = 'random'), NA)))

# lets look at the fits
d_preds_gr <- 
  # stack models into a single column, with an id column for model_name
  pivot_longer(d_fits_gr, names_to = 'model_name', values_to = 'fit', contains('gr')) %>%
  filter(., !is.na(fit)) %>%
  mutate(new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100)))) %>%
  # get rid of original data column
  select(., -data) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(id, habitat, biome, preds, model_name) %>%
  # unlist the preds list column
  unnest(preds) %>%
  filter(.fitted > -1)

ggplot(d_growth, aes(temp, max_percap)) +
  geom_point() +
  facet_wrap(~id, scales = 'free_y') +
  theme_bw() +
  geom_line(aes(temp, .fitted, col = model_name), d_preds_gr) +
  facet_wrap(~id, scales = 'free_y') +
  theme_bw()

# calculate quasi Rsq for each curve
d_fits_gr <- pivot_longer(d_fits_gr, names_to = 'model_name', values_to = 'fit', contains('gr')) %>%
  filter(., !is.na(fit)) %>%
  mutate(r2 = map_dbl(fit, function(x){soilphysics::Rsq(x)$pseudo.R.squared}))

# lets use sharpe schoolfield
params <- filter(d_fits_gr, r2 > 0.5) %>%
  mutate(params = map(fit, calc_params)) %>%
  select(id, habitat, biome, model_name, params) %>%
  unnest(params)

params2 <- filter(params, model_name == 'gr_sharpeschool') %>%
  select(id, habitat, biome, model_name, rmax, topt, ctmax, thermal_safety_margin, thermal_tolerance) %>%
  pivot_longer(., names_to = 'param', values_to = 'value', rmax:thermal_tolerance) %>%
  filter(!is.infinite(value))

ggplot(params2, aes(biome, value)) +
  geom_quasirandom(size = 1, width = 0.2) +
  facet_wrap(~param, scales = 'free_y')


d <- d_growth %>%
  filter(id == 'AUS7')

start_vals <- get_start_vals(d$temp, d$auc, model_name = 'sharpeschoolhigh_1981')

mod <- nls.multstart::nls_multstart(auc~sharpeschoolhigh_1981(temp = temp, r_tref, e, eh, th, tref = 18),
                                    data = d,
                                    iter = c(5,5,5,5),
                                    start_lower = start_vals - 10,
                                    start_upper = start_vals + 10,
                                    #lower = get_lower_lims(d$temp, d$auc, model_name = 'sharpeschoolhigh_1981'),
                                    #upper = get_upper_lims(d$temp, d$auc, model_name = 'sharpeschoolhigh_1981'),
                                    supp_errors = 'Y',
                                    convergence_count = FALSE)

# lets look at the fits
d_preds_gr <- filter(d_fits, !is.na(gr_sharpeschool)) %>%
  mutate(new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100)))) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', contains('gr')) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(id, preds, model_name) %>%
  # unlist the preds list column
  unnest(preds)

ggplot(d_growth, aes(temp, max_percap)) +
  geom_point() +
  geom_line(aes(temp, .fitted, col = model_name), d_preds_gr) +
  facet_wrap(~id, scales = 'free_y') +
  theme_bw()


 
# CODE STOPS ####

# create means and standard errors
d_growth_sum <- d_growth %>%
  group_by(id, temp) %>%
  summarise(mean_percap = mean(max_percap, na.rm = TRUE),
            sd_percap = sd(max_percap, na.rm = TRUE),
            mean_auc = mean(auc, na.rm = TRUE),
            sd_auc = sd(auc, na.rm = TRUE),
            .groups = 'drop') %>%
  mutate(across(starts_with('sd'), ~ifelse(is.na(.x), 0, .x)))

# identify potential outliers
dataViewer(d_growth, x = 'temp', y = 'auc', id_col = 'id')
to_check <- .Last.value
# growth data to check
select(to_check, run, temp, well, id) %>%
  write.csv('data/processed/auc_data_to_check.csv', row.names = FALSE)

# filter out the outliers
d_growth_filt <- get_unclicked(d_growth, to_check)




